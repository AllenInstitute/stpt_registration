/*=========================================================================

  OverwriteImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"
#include "itkRGBPixel.h"
#include "itksys/String.hxx"
#include <string>
#include "itkImageRegionIterator.h"


template <typename ImageType > 
int Pipeline( 
 const char * backgroundFile, 
 const char * foregroundFile,
 const char * foregroundMaskFile,
 const char * outputFile )
{
  typename ImageType::Pointer background;
  typename ImageType::Pointer foreground;
  typename ImageType::Pointer mask;

  itk::idp::ReadImage<ImageType>( backgroundFile, background );
  itk::idp::ReadImage<ImageType>( foregroundFile, foreground );
  itk::idp::ReadImage<ImageType>( foregroundMaskFile, mask );

  typename ImageType::RegionType region = background->GetBufferedRegion();
  typedef itk::ImageRegionIterator<ImageType> Iterator;
  Iterator biter( background, region );
  Iterator fiter( foreground, region );
  Iterator miter( mask, region );

  while ( !biter.IsAtEnd() )
    {
    if( miter.Get() != 0 )
      {
      biter.Set( fiter.Get() );
      }
    ++biter;
    ++fiter;
    ++miter;
    }

  itk::idp::WriteImage<ImageType>( outputFile, background ); 

  return EXIT_SUCCESS;
}




int main( int argc, char *argv[] )
{
  if (argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "backgroundFile foregroundFile foregroundMaskFile outputFile ";
    std::cout << "[format] [dimension] ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string backgroundFile = argv[1];
  std::string foregroundFile = argv[2];
  std::string foregroundMaskFile = argv[3];
  std::string outputFile = argv[4];

  std::string format = "uchar";
  if (argc > 5 )
    {
    format = argv[5];
    }

  unsigned int dimension = 3;
  if ( argc > 6 )
    {
    dimension = atoi( argv[6] );
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
    err = Pipeline<UCHARImageType>( backgroundFile.c_str(), foregroundFile.c_str(), 
                                foregroundMaskFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 && dimension == 2 )
    {
    err = Pipeline<USHORTImageType>( backgroundFile.c_str(), foregroundFile.c_str(), 
                                foregroundMaskFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 && dimension == 2 )
    {
    err = Pipeline<UINTImageType>( backgroundFile.c_str(), foregroundFile.c_str(), 
                                foregroundMaskFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UCHARVolumeType>( backgroundFile.c_str(), foregroundFile.c_str(), 
                                foregroundMaskFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 && dimension == 3 )
    {
    err = Pipeline<USHORTVolumeType>( backgroundFile.c_str(), foregroundFile.c_str(), 
                                foregroundMaskFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UINTVolumeType>( backgroundFile.c_str(), foregroundFile.c_str(), 
                                foregroundMaskFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "float" ) == 0 && dimension == 3 )
    {
    err = Pipeline<FLOATVolumeType>( backgroundFile.c_str(), foregroundFile.c_str(), 
                                foregroundMaskFile.c_str(), outputFile.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else 
    {
    std::cerr << "format: " << format << " and dimension: " << dimension  << " not supported" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


