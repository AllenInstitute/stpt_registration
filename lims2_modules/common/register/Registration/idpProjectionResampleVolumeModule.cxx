/*=========================================================================

  idpProjectionResampleVolumeModule.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#include "idpImageSeries.h"
#include "idpImageSeriesUtilities.h"
#include "idpRegistrationDispatcher.h"
#include "idpRegistrationUtilities.h"

#include "itkMultiThreader.h"

#include <itksys/CommandLineArguments.hxx>
#include <sstream>


int main( int argc, char *argv[] )
{

  // Parse command line arguments using itksys tool
  std::string inputXMLFile = ""; 
  std::string transformXMLFile = "";
  std::string outputDirectory = ""; 
  std::string modelDirectory = "/projects/0378/vol1/informatics/model";
  int threads = 1;
  double redScale = 0;
  double greenScale = 1;
  double blueScale = 0;
  int downsampleFactor = 128;

  {
  typedef itksys::CommandLineArguments argT;

  argT arg;
  arg.Initialize( argc, argv );

  arg.AddArgument( "--modelDirectory", argT::SPACE_ARGUMENT, &modelDirectory, "model directory" );
  arg.AddArgument( "--threads", argT::SPACE_ARGUMENT, &threads, "maximum number of threads" );
  arg.AddArgument( "--redScale", argT::SPACE_ARGUMENT, &redScale, "red scale" );
  arg.AddArgument( "--greenScale", argT::SPACE_ARGUMENT, &greenScale, "green scale" );
  arg.AddArgument( "--blueScale", argT::SPACE_ARGUMENT, &blueScale, "blue scale" );
  arg.StoreUnusedArguments( true );

  if ( !arg.Parse() )
    {
    std::cerr << "Problem parsing command line arguments" << std::endl;
    return EXIT_FAILURE;
    }

  char ** lastArgv;
  int lastArgc;
  arg.GetUnusedArguments( &lastArgc, &lastArgv );
  
  if ( lastArgc < 4 )
    {
    std::cerr << "Error: missing input arguments" << std::endl;
    std::cerr << "Usage: " << argv[0] << " ";
    std::cerr << "[--options] inputXMLFile transformXMLFile outputDirectory";
    std::cerr << std::endl;
    std::cerr << arg.GetHelp() << std::endl;
    arg.DeleteRemainingArguments( lastArgc, &lastArgv );
    return EXIT_FAILURE;
    }

  inputXMLFile = lastArgv[ lastArgc - 3 ];
  transformXMLFile = lastArgv[ lastArgc - 2 ];
  outputDirectory = lastArgv[ lastArgc - 1 ];

  itk::idp::ChangePaths( inputXMLFile );
  itk::idp::ChangePaths( transformXMLFile );
  itk::idp::ChangePaths( outputDirectory );
  itk::idp::ChangePaths( modelDirectory );

  std::cout << "inputXMLFile: " << inputXMLFile << std::endl;
  std::cout << "transformXMLFile: " << transformXMLFile << std::endl;
  std::cout << "outputDirectory: " << outputDirectory << std::endl;
  std::cout << "modelDirectory: " << modelDirectory << std::endl;
  std::cout << "threads: " << threads << std::endl;
  std::cout << "redScale: " << redScale << std::endl;
  std::cout << "greenScale: " << greenScale << std::endl;
  std::cout << "blueScale: " << blueScale << std::endl;
  arg.DeleteRemainingArguments( lastArgc, &lastArgv );
  }

  typedef unsigned short PixelType;
  typedef itk::idp::RegistrationDispatcher<PixelType> DispatcherType;
  typedef DispatcherType::Affine3DTransformType Affine3DTransformType;
  typedef DispatcherType::VolumeType VolumeType;

  DispatcherType::Pointer dispatcher = DispatcherType::New();

  try
    {
    // Limit the number of threads used
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads( threads );

    // Initialization
    dispatcher->SetVerbose( true );
    dispatcher->LoadImageSeriesFromXML( inputXMLFile.c_str() );
    dispatcher->SetModelDirectory( modelDirectory.c_str() );

    dispatcher->GetUtility()->SetDownsampleFactor(downsampleFactor);
    dispatcher->GetUtility()->SetGenerateMask(false);
    dispatcher->GetUtility()->SetUseStandardSize(false);
    dispatcher->SetImageBackground(0);
    dispatcher->SetInvertIntensity(0);
    dispatcher->SetRedScale(redScale);
    dispatcher->SetGreenScale(greenScale);
    dispatcher->SetBlueScale(blueScale);

    dispatcher->Initialize();

    // Read in transforms
    dispatcher->LoadTransformsFromXML( transformXMLFile.c_str() );

    // populate volume
    dispatcher->PopulateVolume();
   // dispatcher->WriteVolume( outputDirectory.c_str() );

    // resample to grid50
    std::string mDir = dispatcher->GetModelDirectory();
    std::string path = mDir;
    path += "/grid50.mhd";
    VolumeType::Pointer grid50;

    Affine3DTransformType::Pointer atov = dispatcher->GetAtlasToVolumeTransform();
    itk::idp::ReadImage<VolumeType>( path.c_str(), grid50 );

    VolumeType::Pointer output;
    VolumeType::Pointer volume = dispatcher->GetVolume();
    itk::idp::ResampleImage<VolumeType>(  volume,
                                          grid50, atov, output, 0, "Linear" );
    path = outputDirectory.c_str();
    path += "/resampledToGrid50.mhd";
    itk::idp::WriteImage<VolumeType>( path.c_str(), output );

    // resample to grid
    path = mDir;
    path += "/grid.mhd";
    VolumeType::Pointer grid;
    itk::idp::ReadImage<VolumeType>( path.c_str(), grid );
    itk::idp::ResampleImage<VolumeType>( volume, grid, atov, output, 0, "Linear" );
    path = outputDirectory.c_str();
    path += "/resampledToGrid.mhd";
    itk::idp::WriteImage<VolumeType>( path.c_str(), output );
    

    return EXIT_SUCCESS;

    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;

}

