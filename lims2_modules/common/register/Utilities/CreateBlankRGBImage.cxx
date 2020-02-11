/*=========================================================================

  CreateBlankRGBImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include <string>

int main( int argc, char *argv[] )
{
  if (argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "outputFile sizex sizey ";
    std::cout << "[red] [blue] [green]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::RGBPixel< unsigned char > PixelType;
  typedef itk::Image< PixelType, 2 >     ImageType;

  std::string outputFile = argv[1];

  ImageType::RegionType region;
  for ( int i = 0; i < 2; i++ )
    {
    region.SetIndex( i, 0 );
    region.SetSize( i, atoi( argv[2+i] ) );
    }  

  PixelType p;

  for ( int k = 0; k < 3; k++ )
    {
    if ( argc > (4+k ) )
      {
      p[k] = atoi( argv[4+k] );
      }
    else
      {
      p[k] = 0;
      }
    }

  
  try
    {
    ImageType::Pointer image = ImageType::New();
    image->SetRegions( region );
    image->Allocate();
    image->FillBuffer( p );

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( image );
    writer->SetFileName( outputFile.c_str() );
    writer->Update();
    
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


