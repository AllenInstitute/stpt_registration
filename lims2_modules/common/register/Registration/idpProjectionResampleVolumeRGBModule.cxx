/*=========================================================================

  idpProjectionResampleVolumeRGBModule.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/


#include "idpProjectionAlignmentDispatcher.h"
#include "itkMultiThreader.h"
#include "itkTransformFileWriter.h"
#include "itkFlipImageFilter.h"
#include "itkFixedArray.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkCastImageFilter.h"
#include <itksys/CommandLineArguments.hxx>
#include "itkCompose3DVectorImageFilter.h"
#include "itkVectorImage.h"

/**
 * Overview: the purpose of this module is to resample the image using the affine transform and nonrigid deformation field
 *
 * Inputs: the path to the folder where there are the outputs from the local alignment.
 *
 * Outputs: the path to the output folder
 * 
 */

int main(int argc, char *argv[]) {

    // Parse command line arguments using itksys tool
    std::string inputXMLFile = "";
    std::string transformXMLFile = "";
    std::string deformationFieldFile = "";
    std::string outputDirectory = "";
    std::string modelDirectory = "/projects/0378/vol1/informatics/model";
    //std::string addonDfmfld = ""; //the addon deformation field
    //std::string addonDfmfld = modelDirectory+"/P56/ARA2TC_dfmfld/dfmfld.nii.gz"; //the addon deformation field


    //std::string outputXMLFile = "";
    int debugLevel = 2;
    int threads = 2;
    int downsampleFactor = 64;
    //int chooseParameters = 1;

    
    {
        typedef itksys::CommandLineArguments argT;

        argT arg;
        arg.Initialize(argc, argv);

        arg.AddArgument("--modelDirectory", argT::SPACE_ARGUMENT, &modelDirectory, "path to model directory");
        arg.AddArgument("--threads", argT::SPACE_ARGUMENT, &threads, "maximum number of threads");
        arg.AddArgument("--downsampleFactor", argT::SPACE_ARGUMENT, &downsampleFactor, "factor of the input image downsampling");
		//arg.AddArgument("--addonDeformationField", argT::SPACE_ARGUMENT, &addonDfmfld, "deformation field to be composed");
		arg.AddArgument("--debugLevel", argT::SPACE_ARGUMENT, &debugLevel, "debugging level");
		

        arg.StoreUnusedArguments(true);

        if (!arg.Parse()) {
            std::cerr << "Problem parsing command line arguments" << std::endl;
            return EXIT_FAILURE;
        }

        char ** lastArgv;
        int lastArgc;
        arg.GetUnusedArguments(&lastArgc, &lastArgv);

        if (lastArgc < 5) {
            std::cerr << "Error: missing input arguments" << std::endl;
            std::cerr << "Usage: " << argv[0] << " ";
            std::cerr << "[--options] inputXMLFile transformXMLFile deformationFieldFile outputDirectory";
            std::cerr << std::endl;
            std::cerr << arg.GetHelp() << std::endl;
            arg.DeleteRemainingArguments(lastArgc, &lastArgv);
            return EXIT_FAILURE;
        }

		inputXMLFile = lastArgv[ lastArgc - 4 ];
        transformXMLFile = lastArgv[ lastArgc - 3 ];
        deformationFieldFile = lastArgv[ lastArgc - 2 ];
        outputDirectory = lastArgv[ lastArgc - 1 ];

		itk::idp::ChangePaths(inputXMLFile);
        itk::idp::ChangePaths(transformXMLFile);
        itk::idp::ChangePaths(deformationFieldFile);
        itk::idp::ChangePaths(outputDirectory);
        itk::idp::ChangePaths(modelDirectory);

		std::cout << "inputXMLFile: " << inputXMLFile << std::endl;
		std::cout << "transformXMLFile: " << transformXMLFile << std::endl;
        std::cout << "deformationFieldFile: " << deformationFieldFile << std::endl;
        std::cout << "outputDirectory: " << outputDirectory << std::endl;
        std::cout << "modelDirectory: " << modelDirectory << std::endl;
        std::cout << "debugLevel: " << debugLevel << std::endl;
        std::cout << "threads: " << threads << std::endl;
        std::cout << "downsampelFactor: " << downsampleFactor << std::endl;
        arg.DeleteRemainingArguments(lastArgc, &lastArgv);
    }
	
	//std::string inputXMLFile = inputDirectory+"/image_input.xml";
	//std::string transformXMLFile = inputDirectory+"/transform_input.xml";	

    // Note: image data is 16-bit
    typedef unsigned short PixelType;
	typedef itk::Image<float, 3>  RegistrationImageType; 		
    typedef itk::idp::ProjectionAlignmentDispatcher<PixelType> DispatcherType;
    DispatcherType::Pointer dispatcher = DispatcherType::New();

    typedef DispatcherType::VolumeType VolumeType;
    std::string fname;


    try {
        // Limit the number of threads used
        itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads);

        // Initialization
        dispatcher->SetVerbose(true);
        dispatcher->LoadImageSeriesFromXML(inputXMLFile.c_str());
        dispatcher->SetModelDirectory(modelDirectory.c_str());
		
        // Volume creation settings:
        // * downsample factor = 64
        // * no mask images were created for this type of data
        // * volume size is the minimal size to house stacked images aligned at (0,0) pixel
        // * volume background is zero
        // * intensity inversion off
        // * use red channel for registration
        dispatcher->GetUtility()->SetDownsampleFactor(downsampleFactor);
        dispatcher->GetUtility()->SetGenerateMask(false);
        dispatcher->GetUtility()->SetUseStandardSize(false);
        dispatcher->SetImageBackground(0);
        dispatcher->SetInvertIntensity(0);
        dispatcher->SetRedScale(1);
        dispatcher->SetGreenScale(0);
        dispatcher->SetBlueScale(0);

        // Initialize all the internal data and setup volume meta-information and sizes.
        // Initially, all transforms are set to identity.
        // Hence the initial reconstruction is a simple stacking aligned at (0,0) pixel
        dispatcher->Initialize();
		
		// Read in transforms
		dispatcher->LoadTransformsFromXML( transformXMLFile.c_str() );
		
		dispatcher->PopulateVolume();
		DispatcherType::Affine3DTransformType::Pointer vtoa = dispatcher->GetAtlasToVolumeTransform();
		vtoa->Print(std::cout);
		
		// get the deformation field
		std::cout<<"Reading the deformation field "<<deformationFieldFile;
		typedef itk::Vector<float, 3>    VectorType;
		typedef itk::Image<VectorType, 3> FieldType;
		typedef itk::ImageFileReader<FieldType> FldReaderType;
		FldReaderType::Pointer fldReader = FldReaderType::New();
		fldReader->SetFileName(deformationFieldFile.c_str());
		fldReader->Update();	
		FieldType::Pointer field = fldReader->GetOutput();		
		std::cout<<" ... done. "<<std::endl;
		
		// get the resolution of deformation field
		//fldReader->GetOutput()->GetLargestPossibleRegion().GetSpacing();
		std::cout<<" The resolution of input deformation field is "<< field->GetSpacing() << std::endl;		
                        
		// read sagittal atlas
		RegistrationImageType::Pointer sagittalAtlas;
                std::ostringstream pixelSize ;
                pixelSize << field->GetSpacing()[0];                
		fname = modelDirectory;
   		fname = fname+ "/P56/atlasVolume/average_template_"+pixelSize.str()+".nrrd";
		itk::idp::ReadImage<RegistrationImageType>( fname.c_str(), sagittalAtlas );
		std::cout<<"finish reading the atlas"<<std::endl;			
		
		// resample to the sagittal atlas affinely	
		typedef itk::ResampleImageFilter< VolumeType, VolumeType >    ResampleFilterType;
		ResampleFilterType::Pointer resampler = ResampleFilterType::New();
		resampler->SetInput( dispatcher->GetVolume() ); 
		
		resampler->SetTransform( vtoa.GetPointer() ); 
		resampler->SetSize( sagittalAtlas->GetLargestPossibleRegion().GetSize() );
		resampler->SetOutputOrigin(  sagittalAtlas->GetOrigin() );
		resampler->SetOutputSpacing( sagittalAtlas->GetSpacing() );
		resampler->SetOutputDirection( sagittalAtlas->GetDirection() );
		resampler->SetDefaultPixelValue( 0 );
		resampler->Update();
				
		// resample to the sagittal atlas using deformation field
		VolumeType::Pointer warpedImg;
		VolumeType::Pointer affineResampled = resampler->GetOutput();
		itk::idp::ResampleImageByDfmfld< VolumeType,RegistrationImageType,FieldType >(affineResampled,sagittalAtlas,field,warpedImg,"Linear" );
		
		// output the affine warped red channel image 
		fname = outputDirectory;
		fname += "/resampled_red.mhd";
		itk::idp::WriteImage< VolumeType >( fname.c_str(), warpedImg );  
		
		// warp and output the green channel image
		dispatcher->SetRedScale(0);
                dispatcher->SetGreenScale(1);
                dispatcher->SetBlueScale(0);
		dispatcher->PopulateVolume();
		resampler->SetInput( dispatcher->GetVolume() ); 
		resampler->Update();
		affineResampled = resampler->GetOutput();
		itk::idp::ResampleImageByDfmfld< VolumeType,RegistrationImageType,FieldType >(affineResampled,sagittalAtlas,field,warpedImg,"Linear" );
		fname = outputDirectory;
		fname += "/resampled_green.mhd";
		itk::idp::WriteImage< VolumeType >( fname.c_str(), warpedImg );  
		
		// warp and output the blue channel image
		dispatcher->SetRedScale(0);
                dispatcher->SetGreenScale(0);
                dispatcher->SetBlueScale(1);
		dispatcher->PopulateVolume();
		resampler->SetInput( dispatcher->GetVolume() ); 
		resampler->Update();
		affineResampled = resampler->GetOutput();
		itk::idp::ResampleImageByDfmfld< VolumeType,RegistrationImageType,FieldType >(affineResampled,sagittalAtlas,field,warpedImg,"Linear" );
		fname = outputDirectory;
		fname += "/resampled_blue.mhd";
		itk::idp::WriteImage< VolumeType >( fname.c_str(), warpedImg );
		

    } catch (itk::ExceptionObject & err) {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}
