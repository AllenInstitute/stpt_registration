/*=========================================================================

  idpProjectionAlignmentModule.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/


#include "idpProjectionAlignmentDispatcher.h"
#include "idpProjectionLocalAlignmentDispatcher.h"
#include "itkMultiThreader.h"
#include "itkTransformFileWriter.h"
#include "itkFlipImageFilter.h"
#include "itkFixedArray.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkCastImageFilter.h"
#include <itksys/CommandLineArguments.hxx>

#include "itkDisplacementFieldCompositionFilter.h"
#include "itkCompose3DVectorImageFilter.h"
#include "itkVectorImage.h"

/**
 * Overview: the purpose of this module is register each tissuecyte projection dataset
 * to the same canonical 3D space as the allen reference atlas.
 *
 * Inputs: An xml file containing meta-informatio about the specimen, image-series and images.
 *
 * Outputs: A warped image as well as the control points that can laterly reproduce the deformation field
 * 
 */

int main(int argc, char *argv[]) {

    // Parse command line arguments using itksys tool
    std::string inputXMLFile = "";
	std::string transformXMLFile = "";
    std::string outputDirectory = "";
    std::string modelDirectory = "/projects/0378/vol1/informatics/model";
	//std::string addonDfmfld = ""; //the addon deformation field
	    
	std::string outputXMLFile = "";
	int debugLevel = 2;
    int threads = 2;
    int downsampleFactor = 64;
	int channelR = 1;
	int channelG = 0;
	int channelB = 0;
	int addonDeformationFieldOn = 0;
	
    //int chooseParameters = 1;    
    
        typedef itksys::CommandLineArguments argT;

        argT arg;
        arg.Initialize(argc, argv);

        arg.AddArgument("--modelDirectory", argT::SPACE_ARGUMENT, &modelDirectory, "path to model directory");
        arg.AddArgument("--threads", argT::SPACE_ARGUMENT, &threads, "maximum number of threads [2]");
        arg.AddArgument("--downsampleFactor", argT::SPACE_ARGUMENT, &downsampleFactor, "factor of the input image downsampling [64]");
		arg.AddArgument("--addonDeformationFieldOn", argT::SPACE_ARGUMENT, &addonDeformationFieldOn, "set to '1' to compose with the addon deformation field [0]");
		arg.AddArgument("--debugLevel", argT::SPACE_ARGUMENT, &debugLevel, "debugging level [2]");
		arg.AddArgument("--channelR", argT::SPACE_ARGUMENT, &channelR, "set 1 to use the red channel [1] ");
		arg.AddArgument("--channelG", argT::SPACE_ARGUMENT, &channelG, "set 1 to use the green channel [0]");
		arg.AddArgument("--channelB", argT::SPACE_ARGUMENT, &channelB, "set 1 to use the blue channel [0]");

        arg.StoreUnusedArguments(true);

        if (!arg.Parse()) {
            std::cerr << "Problem parsing command line arguments" << std::endl;
            return EXIT_FAILURE;
        }

        char ** lastArgv;
        int lastArgc;
        arg.GetUnusedArguments(&lastArgc, &lastArgv);

        if (lastArgc < 4) {
            std::cerr << "Error: missing input arguments" << std::endl;
            std::cerr << "Usage: " << argv[0] << " ";
            std::cerr << "[--options] inputXMLFile transformXMLFile outputDirectory";
            std::cerr << std::endl;
            std::cerr << arg.GetHelp() << std::endl;
            arg.DeleteRemainingArguments(lastArgc, &lastArgv);
            return EXIT_FAILURE;
        }
		
		std::string addonDfmfld;
		
		if (addonDeformationFieldOn==1) {
			addonDfmfld = modelDirectory+"/P56/ARA2TC_dfmfld/dfmfld.nii.gz"; //the addon deformation field
		}

        inputXMLFile = lastArgv[ lastArgc - 3 ];
        transformXMLFile = lastArgv[ lastArgc - 2 ];
        outputDirectory = lastArgv[ lastArgc - 1 ];

        itk::idp::ChangePaths(inputXMLFile);
        itk::idp::ChangePaths(transformXMLFile);
        itk::idp::ChangePaths(outputDirectory);
        itk::idp::ChangePaths(modelDirectory);

        std::cout << "inputXMLFile: " << inputXMLFile << std::endl;
		std::cout << "transformXMLFile: " << transformXMLFile << std::endl;
        std::cout << "outputDirectory: " << outputDirectory << std::endl;
        std::cout << "modelDirectory: " << modelDirectory << std::endl;
        std::cout << "debugLevel: " << debugLevel << std::endl;
        std::cout << "threads: " << threads << std::endl;
        std::cout << "downsampelFactor: " << downsampleFactor << std::endl;
        arg.DeleteRemainingArguments(lastArgc, &lastArgv);
    

    // Note: image data is 16-bit
    typedef unsigned short PixelType;
	typedef itk::Image<float, 3>  RegistrationImageType; 		
    typedef itk::idp::ProjectionLocalAlignmentDispatcher<PixelType> DispatcherType;
    DispatcherType::Pointer dispatcher = DispatcherType::New();

    typedef DispatcherType::VolumeType VolumeType;
    typedef DispatcherType::IOImageType IOImageType;
	typedef DispatcherType::RegistrationImageType RegistrationImageType;
	
	
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
        
		dispatcher->SetRedScale(channelR);
		dispatcher->SetGreenScale(channelG);
		dispatcher->SetBlueScale(channelB);
		

        // Initialize all the internal data and setup volume meta-information and sizes.
        // Initially, all transforms are set to identity.
        // Hence the initial reconstruction is a simple stacking aligned at (0,0) pixel
        dispatcher->Initialize();
		
		// Read in transforms
		dispatcher->LoadTransformsFromXML( transformXMLFile.c_str() );
		
		      
		// do registration
		dispatcher->deformableReg(modelDirectory, outputDirectory, debugLevel);        
		
		// output the dfmfld
		if (debugLevel > 1) {
			typedef itk::Vector<float, 3>    VectorType;
			typedef itk::Image<VectorType, 3> FieldType;
			typedef itk::ImageFileReader<FieldType> FldReaderType;
			typedef itk::ImageFileWriter<FieldType> FldWriterType;
 				
			FldWriterType::Pointer fldWriter = FldWriterType::New();
			std::string fldOutput = outputDirectory+"/dfmfld.mhd";
			fldWriter->SetFileName(fldOutput.c_str());
		
			RegistrationImageType::Pointer fldX = RegistrationImageType::New();
			RegistrationImageType::Pointer fldY = RegistrationImageType::New();
			RegistrationImageType::Pointer fldZ = RegistrationImageType::New();
			
			char* gridFile=new char[200];		
			strcpy(gridFile,(outputDirectory+"/grids.raw").c_str());
			dispatcher->ResampleDeformationFld(gridFile, modelDirectory, fldX, fldY, fldZ);
		
			typedef itk::Compose3DVectorImageFilter<RegistrationImageType,FieldType> ImageToVectorImageFilterType;
			ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
			imageToVectorImageFilter->SetInput(0, fldX);
			imageToVectorImageFilter->SetInput(1, fldY);
			imageToVectorImageFilter->SetInput(2, fldZ);
			imageToVectorImageFilter->Update();
			
			// compose with the provided dfmfld
			if (!addonDfmfld.empty()) {
				std::cout<<"Compose with the addon deformation field "<<addonDfmfld<<" ..."<<std::endl;
			
				//read in the addon
				FldReaderType::Pointer fldReader = FldReaderType::New();
				fldReader->SetFileName(addonDfmfld.c_str()); 
  			
				typedef itk::DisplacementFieldCompositionFilter<FieldType, FieldType> ComposerType;
				ComposerType::Pointer composer = ComposerType::New();
				composer->SetInput( 1, fldReader->GetOutput() );
				composer->SetInput( 0, imageToVectorImageFilter->GetOutput() );
				//composer->Print( std::cout );
				composer->Update();
			
				fldWriter->SetInput(composer->GetOutput());
				fldWriter->Update();
			}
			else {
			
				fldWriter->SetInput(imageToVectorImageFilter->GetOutput());
				fldWriter->Update();
			}
			std::cout<<"Write out the deformation field ... done."<<std::endl;
		}
				
        // Write out a sagittal slice visual inspection
        if (debugLevel > 0) {
		
			// get the warped image			
			std::string warpedFile =outputDirectory+"/deformed.nii.gz";
			IOImageType::Pointer warpedImage = IOImageType::New();
			itk::idp::ReadImage< IOImageType >( warpedFile.c_str(), warpedImage );		
			
			double shift = 0.0;
			double scale = 0.2;
			typedef DispatcherType::ImageType SliceType;
			typedef itk::Image<unsigned char, 2 > OutputSliceType;
            
			SliceType::Pointer slice;
			OutputSliceType::Pointer sliceOutput;
            itk::idp::ExtractSlice< IOImageType,SliceType >( warpedImage,slice,250,2 );
			itk::idp::ShiftScale< SliceType, OutputSliceType > (slice, sliceOutput, shift, scale);

			fname = outputDirectory;
            fname += "/sliceSagittal.png";
            itk::idp::WriteImage< OutputSliceType > (fname.c_str(), sliceOutput);       

			if (debugLevel < 3) {
				std::remove(warpedFile.c_str());
			}
        }


    } catch (itk::ExceptionObject & err) {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}

