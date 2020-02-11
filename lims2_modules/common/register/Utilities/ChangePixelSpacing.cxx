/*=========================================================================

  ChangePixelSpacing.cxx
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"

#include "itksys/String.hxx"
#include <string>

template <typename ImageType > 
int Pipeline( 
 const char * infile, 
 const itk::FixedArray<double,3> & spacing,
 const itk::FixedArray<double,3> & variance,
 const char * outfile,
 const char * interpolation )
{
  typename ImageType::Pointer input;
  typename ImageType::Pointer output;

  itk::idp::ReadImage<ImageType>( infile, input );
  itk::idp::ChangePixelSpacing<ImageType>( input, spacing, variance, output, interpolation );
  itk::idp::WriteImage<ImageType>( outfile, output );

  return EXIT_SUCCESS;
}


int main( int argc, char *argv[] )
{
  if (argc < 9 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile outputFile ";
    std::cout << "x-spacing y-spacing z-spacing ";
    std::cout << "x-variance y-variance z-variance ";
    std::cout << "[interpolation] [format]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string infile = argv[1];
  std::string outfile = argv[2];
  itk::FixedArray<double,3> spacing;
  itk::FixedArray<double,3> variance;

  for( unsigned int j = 0; j < 3; j++ )
    {
    spacing[j] = atof(argv[3+j]);
    variance[j] = atof(argv[6+j]);
    }

  std::string interpolation = "Linear";
  if ( argc > 9 )
    {
    interpolation = argv[9];
    }

  std::string format = "UCHAR";
  if (argc > 10 )
    {
    format = argv[10];
    }

  unsigned int dimension = 3;

  typedef itk::Image< unsigned char, 3 >      UCHARVolumeType;
  typedef itk::Image< unsigned short, 3 >     USHORTVolumeType;
  typedef itk::Image< unsigned int,  3 >      UINTVolumeType;
  typedef itk::Image< float, 3>               FLOATVolumeType;

  int err;

  if( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UCHARVolumeType>( infile.c_str(), spacing, variance, outfile.c_str(), interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 && dimension == 3 )
    {
    err = Pipeline<USHORTVolumeType>( infile.c_str(), spacing, variance, outfile.c_str(), interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 && dimension == 3 )
    {
    err = Pipeline<UINTVolumeType>( infile.c_str(), spacing, variance, outfile.c_str(), interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "float" ) == 0 && dimension == 3 )
    {
    err = Pipeline<FLOATVolumeType>( infile.c_str(), spacing, variance, outfile.c_str(), interpolation.c_str() );
    if ( err ) return EXIT_FAILURE;
    }
  else 
    {
    std::cerr << "format: " << format << " not supported" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


