/*=========================================================================

  idpProjectionAlignmentModule.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/


#include "idpProjectionAlignmentDispatcher.h"
#include "itkMultiThreader.h"
#include "itkTransformFileWriter.h"
#include <itksys/CommandLineArguments.hxx>

/**
 * Overview: the purpose of this module is register each tissuecyte projection dataset
 * to the same canonical 3D space as the allen reference atlas.
 *
 * Inputs: An xml file containing meta-informatio about the specimen, image-series and images.
 *
 * Outputs: An xml file containing a 2D transform for each image to the reconstructed volume
 * space of the specimen and a 3D transform for each image-series mapping the reconstructed
 * volume to the canonical 3D space.
 *
 */

int main(int argc, char *argv[]) {

    // Parse command line arguments using itksys tool
    std::string inputXMLFile = "";
    std::string outputXMLFile = "";
    std::string outputDirectory = "";
    std::string modelDirectory = "/projects/0378/vol1/informatics/model";
    int debugLevel = 1;
    int threads = 1;
    int downsampleFactor = 64;
    double shift = 0.0;
    double scale = 0.2;
    int chooseParameters = 1;
    std::string rigidRegParaFile = "rigidVersor3DParameters-1.xml";
    std::string affineRegParaFile = "centeredAffineParameters.xml";
	
	int channelR = 1;
	int channelG = 0;
	int channelB = 0;
    
    {
        typedef itksys::CommandLineArguments argT;

        argT arg;
        arg.Initialize(argc, argv);

        arg.AddArgument("--modelDirectory", argT::SPACE_ARGUMENT, &modelDirectory, "path to model directory");
        arg.AddArgument("--debugLevel", argT::SPACE_ARGUMENT, &debugLevel, "debugging level");
        arg.AddArgument("--threads", argT::SPACE_ARGUMENT, &threads, "maximum number of threads");
        arg.AddArgument("--downsampleFactor", argT::SPACE_ARGUMENT, &downsampleFactor, "factor of the input image downsampling");
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
            std::cerr << "[--options] inputXMLFile outputXMLFile outputDirectory";
            std::cerr << std::endl;
            std::cerr << arg.GetHelp() << std::endl;
            arg.DeleteRemainingArguments(lastArgc, &lastArgv);
            return EXIT_FAILURE;
        }

        inputXMLFile = lastArgv[ lastArgc - 3 ];
        outputXMLFile = lastArgv[ lastArgc - 2 ];
        outputDirectory = lastArgv[ lastArgc - 1 ];

        itk::idp::ChangePaths(inputXMLFile);
        itk::idp::ChangePaths(outputXMLFile);
        itk::idp::ChangePaths(outputDirectory);
        itk::idp::ChangePaths(modelDirectory);

        std::cout << "inputXMLFile: " << inputXMLFile << std::endl;
        std::cout << "outputXMLFile: " << outputXMLFile << std::endl;
        std::cout << "outputDirectory: " << outputDirectory << std::endl;
        std::cout << "modelDirectory: " << modelDirectory << std::endl;
        std::cout << "debugLevel: " << debugLevel << std::endl;
        std::cout << "threads: " << threads << std::endl;
        std::cout << "downsampelFactor: " << downsampleFactor << std::endl;
        arg.DeleteRemainingArguments(lastArgc, &lastArgv);
    }

    // Note: image data is 16-bit
    typedef unsigned short PixelType;
    typedef itk::idp::ProjectionAlignmentDispatcher<PixelType> DispatcherType;
    DispatcherType::Pointer dispatcher = DispatcherType::New();

    typedef DispatcherType::ImageType ImageType;
    typedef DispatcherType::VolumeType VolumeType;
    std::string fname;

    // choose between the parameter files. xxx1.xml: optimal; xxx2.xml: suboptimal
    choosePara:

    try {
        // Limit the number of threads used
        itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads);

        // Initialization
        dispatcher->SetVerbose(true);
        dispatcher->LoadImageSeriesFromXML(inputXMLFile.c_str());
        dispatcher->SetModelDirectory(modelDirectory.c_str());
        dispatcher->SetGlobalAlignmentModelDirectory();

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
        //dispatcher->SetRedScale(1);
        //dispatcher->SetGreenScale(0);
        //dispatcher->SetBlueScale(0);
		
		dispatcher->SetRedScale(channelR);
		dispatcher->SetGreenScale(channelG);
		dispatcher->SetBlueScale(channelB);

        // Initialize all the internal data and setup volume meta-information and sizes.
        // Initially, all transforms are set to identity.
        // Hence the initial reconstruction is a simple stacking aligned at (0,0) pixel
        dispatcher->Initialize();

        // This section is for preparing data form the model directory
        // To-do: move this somewhere else
        if (debugLevel > 10) {
            typedef itk::Image<unsigned char, 2 > OutputImageType;
            // Create slice view and write out
            dispatcher->ReadAtlasData();
            ImageType::Pointer slice;
            VolumeType::Pointer atlas = dispatcher->GetAtlasVolume();
            dispatcher->ExtractViewSlice(atlas, 1, slice);
            OutputImageType::Pointer output;
            itk::idp::ShiftScale<ImageType, OutputImageType > (slice, output, shift, scale);

            fname = outputDirectory;
            fname += "/slice.png";
            itk::idp::WriteImage<OutputImageType > (fname.c_str(), output);

            // write out z projection image
            VolumeType::Pointer proj;
            itk::idp::AverageProjection<VolumeType > (atlas, proj, 2);
            fname = outputDirectory;
            fname += "/avgProj.png";
            itk::idp::WriteImage<VolumeType > (fname.c_str(), proj);

            exit(0);
        }

        // Function PopulateVolume creates the reconstructed volume based on
        // the current transform parameters
        dispatcher->PopulateVolume();

        // Function ReadAtlasData reads in atlas volume from the model directory
        dispatcher->ReadAtlasData();

        // Resample the atlas using initial transform and compute metric
           dispatcher->ResampleAtlas();
           dispatcher->ComputeMetric();
           dispatcher->Checkpoint();

        // Use command-line option (--debugLevel 2) to dump out intermediate volumes
        // For production use default debugLevel = 1
        if (debugLevel > 1) {
            dispatcher->WriteVolume(outputDirectory.c_str(), true);
            dispatcher->WriteResampledAtlas(outputDirectory.c_str(), true);
        }

        // set verbose before actrual registration
        if (debugLevel < 2) {
            dispatcher->SetVerbose(false);
        }   
		
		// detect the software pause of the scanner which usually causes the image shift between sections
		//dispatcher->WriteVolume(outputDirectory.c_str(), true); 
		dispatcher->DetectSoftwarePause();
		dispatcher->PopulateVolume();
		//dispatcher->WriteVolume(outputDirectory.c_str(), true);
		//return EXIT_SUCCESS; // exis early only for debug, remove this when done.
        
        // first do rigid versor as initialization
        dispatcher->RigidVersor3DRegistration(rigidRegParaFile.c_str());
        // then do the 12-parameter affine
        dispatcher->CenteredAffineVolumeRegistration(affineRegParaFile.c_str());
        
        // set verbose after actrual registration
        if (debugLevel > 0) {
            dispatcher->SetVerbose(true);
        }

        dispatcher->ResampleAtlas();
        dispatcher->ResampleVolume();
        
        // check if the alignment improves significantly, otherwise restore the transform at last checkpoint
        dispatcher->ComputeMetric();
        if ((dispatcher->Rollback())&&(chooseParameters==1)) //try sub-optimal parameters instead
	{
            dispatcher->ResampleAtlas();
            chooseParameters=2;
            rigidRegParaFile = "rigidVersor3DParameters-2.xml";
            affineRegParaFile = "centeredAffineParameters2.xml";
            std::cout<<"changing the parameter file to sub-optimal one ......"<<std::endl;
            goto choosePara;
        }
	else if ((dispatcher->Rollback())&&(chooseParameters==2))
	{
	    dispatcher->ResampleAtlas();
	    std::cout<<"Can't converge on sub-optimal parameters either...revert the transform"<<std::endl;
	}

        // write out the volume and atlas after registration
        if (debugLevel > 1) {
            dispatcher->WriteResampledAtlas(outputDirectory.c_str(), true);
            dispatcher->WriteResampledVolume(outputDirectory.c_str(), true);
        }

        // Write out a horizontal slice and z-projection for debugging purposes
        if (debugLevel > 0) {
            typedef itk::Image<unsigned char, 2 > OutputImageType;

            // Create slice view and write out
            dispatcher->ResampleVolume();
            ImageType::Pointer slice;
            dispatcher->ExtractViewSlice(2, slice);
            OutputImageType::Pointer output;
            itk::idp::ShiftScale<ImageType, OutputImageType > (slice, output, shift, scale);
            fname = outputDirectory;
            fname += "/sliceY.png";
            itk::idp::WriteImage<OutputImageType > (fname.c_str(), output);

            // write out z projection image
            ImageType::Pointer proj = dispatcher->GetZProjection();
            itk::idp::ShiftScale<ImageType, OutputImageType > (proj, output, shift, scale);
            fname = outputDirectory;
            fname += "/projectionZ.png";
            itk::idp::WriteImage<OutputImageType > (fname.c_str(), output);
        }

		// Write out the transform in itk format
        if (debugLevel > 1) {
			itk::TransformFileWriter::Pointer tranWriter = itk::TransformFileWriter::New();
			tranWriter->SetInput(dispatcher->GetVolumeToAtlasTransform());
			fname = outputDirectory;
            		fname += "/itkTran.txt";
			tranWriter->SetFileName(fname.c_str());
			tranWriter->Update(); 

			tranWriter->SetInput(dispatcher->GetAtlasToVolumeTransform());
			fname = outputDirectory;
            		fname += "/itkTranInv.txt";
			tranWriter->SetFileName(fname.c_str());
			tranWriter->Update(); 
		}

        // Write transform to XML file
        // First compose the registration with the bridge-to-canonical_space transform.
        // These functions take care of adding in the image resolution scale at the sub-image level.
        dispatcher->ComposeWithAtlasToGridTransform();
        dispatcher->WriteTransformsToXML(outputXMLFile.c_str());

    } catch (itk::ExceptionObject & err) {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;

}

