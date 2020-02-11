/*=========================================================================

  AffineCorrelation2DRegistration.cxx
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"
#include "idpAffineCorrelationVolumeRegistration.h"
#include <string>
#include <iostream>


int main( int argc, char *argv[] )
{
  if (argc < 7 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile targetFile parameterFile resampleFile transformFile metricFile";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFile = argv[1];
  std::string targetFile = argv[2];
  std::string parameterFile = argv[3];
  std::string resampleFile = argv[4];
  std::string transformFile = argv[5];
  std::string metricFile = argv[6];

  typedef itk::Image<unsigned char, 2> ImageType;
  typedef itk::idp::AffineCorrelationVolumeRegistration<ImageType> RegistrationType;
  typedef RegistrationType::TransformType TransformType;

  try
   {

   ImageType::Pointer input;
   ImageType::Pointer target;
   ImageType::Pointer resample;
   TransformType::Pointer transform;
   ImageType::Pointer inputMask;
   ImageType::Pointer targetMask;

   itk::idp::ReadImage<ImageType>( inputFile.c_str(), input );
   itk::idp::ReadImage<ImageType>( targetFile.c_str(), target );
   itk::idp::MakeImage<ImageType,ImageType>( input, inputMask, 255 );
   itk::idp::MakeImage<ImageType,ImageType>( target, targetMask, 255 );

   RegistrationType::Pointer reg = RegistrationType::New();

   reg->SetFixedVolume( target );
   reg->SetMovingVolume( input );
   reg->SetFixedMask( inputMask );
   reg->SetMovingMask( targetMask );

   reg->LoadParametersFromXML( parameterFile.c_str() );

   transform = TransformType::New();
   transform->SetIdentity();
   reg->SetInputTransform( transform );

   reg->Compute();

   itk::idp::ResampleImage<ImageType>( input, target, transform, resample, 255, "Linear" );
    
   itk::idp::WriteImage<ImageType>( resampleFile.c_str(), resample );
   itk::idp::WriteTransform<TransformType>( transformFile.c_str(), transform );

   double correlation;
   itk::idp::ComputeCorrelation<ImageType>( resample, target, inputMask, inputMask, correlation );

   std::ofstream output( metricFile.c_str(), std::ios::out );
    if ( output.fail() )
      {
      std::cout << "Error while opening output file " << metricFile << std::endl;
      return EXIT_FAILURE;
      }
   output << correlation << std::endl;
   output.close();


   }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown error" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


