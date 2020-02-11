/*=========================================================================

  DownsampleImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpGriddingUtilities.h"
#include "idpRegistrationUtilities.h"

#include "itkImage.h"

int main( int argc, char *argv[] )
{
  if ( argc != 4 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile factor outputFile ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string infile = argv[1];
  double factor = 1.0;
  factor = atof( argv[2] );
  std::string outfile = argv[3];

  typedef itk::Image< float, 3 > ImageType;

  ImageType::Pointer input;
  itk::idp::ReadImage<ImageType>( infile.c_str(), input );

  double sigma_x = input->GetSpacing()[0] * factor / 2.0;
  double sigma_y = input->GetSpacing()[1] * factor / 2.0;
  double variance[] = { sigma_x * sigma_x, sigma_y * sigma_y, 0 };
  std::cout << "blurring image" << std::endl;
  itk::idp::GaussianSmoothImage<ImageType>(input, input, variance );

  std::cout << "downsampling image" << std::endl;
  itk::idp::DownsampleImage<ImageType>(input, input, factor);

  itk::idp::WriteImage<ImageType>( outfile.c_str(), input, true );

  return EXIT_SUCCESS;
}
