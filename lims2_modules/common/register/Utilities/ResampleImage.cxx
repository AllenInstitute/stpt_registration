/*=========================================================================

  ResampleImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"
#include "itkRGBPixel.h"
#include "itksys/String.hxx"
#include <string>
#include "itkExtendedImageIOFactory.h"


template <typename ImageType, typename TransformType > 
int Pipeline( 
 const char * infile, 
 const char * transformfile,
 const char * reffile,
 const char * outfile,
 double       pad,
 const char * interpolation )
{
  typename ImageType::Pointer input;
  typename ImageType::Pointer ref;
  typename ImageType::Pointer output;
  typename TransformType::Pointer trans;

  itk::idp::ReadImage<ImageType>( infile, input );
  itk::idp::ReadImage<ImageType>( reffile, ref );
  itk::idp::ReadTransform<TransformType>( transformfile, trans );
  
  itk::idp::ResampleImage<ImageType>( input, ref, trans, output, pad, interpolation );

  itk::idp::WriteImage<ImageType>( outfile, output ); 

  return EXIT_SUCCESS;
}

template <typename ImageType, typename TransformType > 
int VectorPipeline( 
 const char * infile, 
 const char * transformfile,
 const char * reffile,
 const char * outfile,
 double       pad,
 const char * interpolation )
{
  typename ImageType::Pointer input;
  typename ImageType::Pointer ref;
  typename ImageType::Pointer output;
  typename TransformType::Pointer trans;

  itk::idp::ReadImage<ImageType>( infile, input );
  itk::idp::ReadImage<ImageType>( reffile, ref );
  itk::idp::ReadTransform<TransformType>( transformfile, trans );
  
  itk::idp::VectorResampleImage<ImageType>( input, ref, trans, output, pad, interpolation );

  itk::idp::WriteImage<ImageType>( outfile, output ); 

  return EXIT_SUCCESS;
}


int main( int argc, char *argv[] )
{
  if (argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile transformFile refFile outputFile ";
    std::cout << "[pad] [interpolation] [format] [dimension]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string infile = argv[1];
  std::string transformfile = argv[2];
  std::string reffile = argv[3];
  std::string outfile = argv[4];

  double pad = 0;
  if ( argc > 5 )
    {
    pad = atof( argv[5] );
    }

  std::string interpolation = "linear";
  if ( argc > 6 )
    {
    interpolation = argv[6];
    }

  std::string format = "uchar";
  if (argc > 7 )
    {
    format = argv[7];
    }

  unsigned int dimension = 3;
  if ( argc > 8 )
    {
    dimension = atoi( argv[8] );
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
    err = Pipeline<UCHARImageType,Transform2DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 && dimension == 2 )
    {
    err = Pipeline<USHORTImageType,Transform2DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 && dimension == 2 )
    {
    err = Pipeline<UINTImageType,Transform2DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UCHARVolumeType,Transform3DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 && dimension == 3 )
    {
    err = Pipeline<USHORTVolumeType,Transform3DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UINTVolumeType,Transform3DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "float" ) == 0 && dimension == 3 )
    {
    err = Pipeline<FLOATVolumeType,Transform3DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "rgb" ) == 0 && dimension == 2 )
    {
    err = VectorPipeline<RGBImageType,Transform2DType>( infile.c_str(), transformfile.c_str(), reffile.c_str(), 
                                    outfile.c_str(), pad, interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }    
  else 
    {
    std::cerr << "format: " << format << " and dimension: " << dimension  << " not supported" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


