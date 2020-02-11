/*=========================================================================

  ShiftScaleRGBImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include <string>
#include "itkNthElementPixelAccessor.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionIterator.h"
#include <vector>

int main( int argc, char *argv[] )
{
  if (argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputImage outputImage shift scale";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::RGBPixel< unsigned char > PixelType;
  typedef itk::Image< PixelType, 2 >     ImageType;
  typedef itk::Image< unsigned char, 2 > UCharImageType;
  
  std::string inputFile = argv[1];
  std::string outputFile = argv[2];
  double shift = atof( argv[3] );
  double scale = atof( argv[4] );
  
  try
    {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName( inputFile.c_str() );
    reader->Update();

    ImageType::Pointer image = reader->GetOutput();
    image->DisconnectPipeline();

    ImageType::RegionType region = image->GetBufferedRegion();

    typedef itk::NthElementPixelAccessor< unsigned char, PixelType > AccessorType;
    typedef itk::ImageAdaptor< ImageType, AccessorType > AdaptorType;
    typedef itk::ShiftScaleImageFilter< AdaptorType, UCharImageType > FilterType;

    for ( int k = 0; k < 3; k++ )
      {
      AdaptorType::Pointer adaptor = AdaptorType::New();
      adaptor->SetImage( image );
      adaptor->GetPixelAccessor().SetElementNumber( k );

      FilterType::Pointer filter = FilterType::New();
      filter->SetInput( adaptor );
      filter->SetShift( shift );
      filter->SetScale( scale );
      filter->Update(); 

      typedef itk::ImageRegionIterator< UCharImageType > InputIterator;
      typedef itk::ImageRegionIterator< ImageType > OutputIterator;
      InputIterator it( filter->GetOutput(), region );
      OutputIterator ot( image, region );

      while( !it.IsAtEnd() )
        {
        PixelType p = ot.Get();
        p[k] = it.Get();
        ot.Set( p );
        ++it;
        ++ot;
        }

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

