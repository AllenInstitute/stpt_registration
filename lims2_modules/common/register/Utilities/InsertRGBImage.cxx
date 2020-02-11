/*=========================================================================

  InsertRGBImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include <string>
#include "itkImageRegionIterator.h"

int main( int argc, char *argv[] )
{
  if (argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "backgroundImage insertImage outputImage ";
    std::cout << "startx starty";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::RGBPixel< unsigned char > PixelType;
  typedef itk::Image< PixelType, 2 >     ImageType;

  std::string backgroundFile = argv[1];
  std::string insertFile = argv[2];
  std::string outputFile = argv[3];

  ImageType::IndexType start;
  for ( int i = 0; i < 2; i++ )
    {
    start[i] = atoi( argv[4+i] );
    }  
  
  try
    {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader1 = ReaderType::New();
    ReaderType::Pointer reader2 = ReaderType::New();

    reader1->SetFileName( backgroundFile.c_str() );
    reader2->SetFileName( insertFile.c_str() );

    reader1->Update();
    reader2->Update();

    ImageType::Pointer image = reader1->GetOutput();
    image->DisconnectPipeline();

    ImageType::Pointer insert = reader2->GetOutput();
    insert->DisconnectPipeline();

    ImageType::RegionType bRegion = image->GetBufferedRegion();
    ImageType::RegionType iRegion = insert->GetBufferedRegion();
    ImageType::RegionType cRegion = iRegion;
    cRegion.SetIndex( start );
    
    bool cropped = cRegion.Crop( bRegion );
    if (cropped)
      {
      iRegion.SetSize( cRegion.GetSize() );
      }

    if ( iRegion.GetSize() != cRegion.GetSize() )
      {
      std::cout << "regions of different size" << std::endl;
      std::cout << "insert region: " << iRegion << std::endl;
      std::cout << "target region: " << cRegion << std::endl;
      return EXIT_FAILURE;
      }

    typedef itk::ImageRegionIterator<ImageType> Iterator;
    Iterator bit( image, cRegion );
    Iterator iit( insert, iRegion );

    while( !bit.IsAtEnd() )
      {
      bit.Set( iit.Get() );
      ++bit;
      ++iit;
      }

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

