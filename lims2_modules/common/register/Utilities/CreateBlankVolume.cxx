/*=========================================================================

  CreateBlankVolume.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include <string>

int main( int argc, char *argv[] )
{
  if (argc < 5)
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "outputFile ";
    std::cout << "size0 size1 size2 ";
    std::cout << "[spacing0 spacing1 spacing2] ";
    std::cout << "[origin0 origin1 origin2] ";
    std::cout << "[value]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef unsigned char                  PixelType;
  typedef itk::Image< PixelType, 3 >     ImageType;

  std::string outputFile = argv[1];

  ImageType::RegionType region;
  ImageType::SpacingType spacing;
  ImageType::PointType origin;
  PixelType value = 0;

  spacing.Fill( 1.0 );
  origin.Fill( 0.0 );

  for ( int i = 0; i < 3; i++ )
    {
    region.SetIndex( i, 0 );
    region.SetSize( i, atoi( argv[2+i] ) );
    }  

  if ( argc > 7 )
    {
    for ( int i = 0; i < 3; i++ )
      {
      spacing[i] = atof( argv[5+i] );
      }
    }

  if ( argc > 10 )
    {
    for ( int i = 0; i < 3; i++ )
      {
      origin[i] = atof( argv[8+i] );
      }
    }

  if ( argc > 11 )
    {
    value = atoi( argv[11] );
    }
  
  try
    {
    ImageType::Pointer image = ImageType::New();
    image->SetRegions( region );
    image->SetSpacing( spacing );
    image->SetOrigin( origin );
    image->Allocate();
    image->FillBuffer( value );

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


