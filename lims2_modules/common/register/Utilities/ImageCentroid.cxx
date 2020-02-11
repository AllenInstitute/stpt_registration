/*=========================================================================

  ImageCentroid.cxx 

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"
#include "itkImageMomentsCalculator.h"
#include "itksys/String.hxx"
#include <string>
#include "itkImageRegionIterator.h"
#include <fstream>
#include "itkLinearInterpolateImageFunction.h"


template <typename ImageType > 
int Pipeline( 
 const char * imageFile, 
 const char * outputFile )
{
  typename ImageType::Pointer image;

  itk::idp::ReadImage<ImageType>( imageFile, image );
  itk::idp::BinaryThreshold<ImageType>( image, image, 0, 0, 0, 1 );

  typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( image );
  
  typedef itk::ImageMomentsCalculator<ImageType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( image );

  double totalMass = 0.0;
  CalculatorType::VectorType cg;
  cg.Fill( 0.0 );

  double inMask = 0;

  try
    {
    calculator->Compute();

    std::cout << calculator->GetTotalMass() << std::endl;
    std::cout << calculator->GetCenterOfGravity() << std::endl; 

    totalMass = calculator->GetTotalMass();
    for( unsigned int j = 0; j < ImageType::ImageDimension; j++ )
      {
      totalMass *= image->GetSpacing()[j];
      }
    cg = calculator->GetCenterOfGravity();

    typename ImageType::PointType point;
    for( unsigned int j = 0; j < ImageType::ImageDimension; j++ )
      {
      point[j] = cg[j];
      }

    inMask = interpolator->Evaluate( point );

    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    }        


  std::ofstream output( outputFile, std::ios::out );
  if ( !output.fail() )
    {
    output << "total_mass: " << totalMass << std::endl;
    for( unsigned int j = 0; j < ImageType::ImageDimension; j++ )
      {
      output << "cg_ " << j << ": " << cg[j] << std::endl;
      }
    output << "in_structure: " << inMask << std::endl;
    output.close();
    }

  return EXIT_SUCCESS;
}




int main( int argc, char *argv[] )
{
  if (argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "imageFile outputFile ";
    std::cout << "[format] [dimension] ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string imageFile = argv[1];
  std::string outputFile = argv[2];

  std::string format = "uchar";
  if (argc > 3 )
    {
    format = argv[3];
    }

  unsigned int dimension = 3;
  if ( argc > 4 )
    {
    dimension = atoi( argv[4] );
    }


  typedef itk::Image< unsigned char, 2 >      UCHARImageType;
  typedef itk::Image< unsigned short, 2 >     USHORTImageType;
  typedef itk::Image< unsigned int,  2 >      UINTImageType;
  typedef itk::Image< unsigned char, 3 >      UCHARVolumeType;
  typedef itk::Image< unsigned short, 3 >     USHORTVolumeType;
  typedef itk::Image< unsigned int,  3 >      UINTVolumeType;
  typedef itk::RGBPixel<unsigned char>        RGBPixelType;
  typedef itk::Image< RGBPixelType, 2 >       RGBImageType;
  typedef itk::Image< float, 3>               FLOATVolumeType;
  typedef itk::Transform<double,2,2>          Transform2DType;
  typedef itk::Transform<double,3,3>          Transform3DType;

  int err;

  if( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0 && dimension == 2 )
    {
    err = Pipeline<UCHARImageType>( imageFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 && dimension == 2 )
    {
    err = Pipeline<USHORTImageType>( imageFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 && dimension == 2 )
    {
    err = Pipeline<UINTImageType>( imageFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UCHARVolumeType>( imageFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 && dimension == 3 )
    {
    err = Pipeline<USHORTVolumeType>( imageFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UINTVolumeType>( imageFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "float" ) == 0 && dimension == 3 )
    {
    err = Pipeline<FLOATVolumeType>( imageFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else 
    {
    std::cerr << "format: " << format << " and dimension: " << dimension  << " not supported" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


