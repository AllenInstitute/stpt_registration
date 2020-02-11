/*=========================================================================

  OverlayRGBImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include "itkJP2ImageIO.h"
#include "itkExtendedImageIOFactory.h"


int main( int argc, char *argv[] )
{
  if (argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputImage overlayImage outputImage";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  // Need to load in the additional factory which includes jp2 
  itk::ExtendedImageIOFactory::RegisterBuiltInFactories(); 


  typedef itk::RGBPixel< unsigned char > PixelType;
  typedef itk::Image< PixelType, 2 >     ImageType;
  typedef itk::RGBAPixel< unsigned char > OverlayPixelType;
  typedef itk::Image< OverlayPixelType, 2 > OverlayImageType;
  
  std::string inputFile = argv[1];
  std::string overlayFile = argv[2];
  std::string outputFile = argv[3];
  
  try
    {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName( inputFile.c_str() );
    reader->Update();
    
    typedef itk::ImageFileReader< OverlayImageType > OverlayReaderType;
    OverlayReaderType::Pointer overlayReader = OverlayReaderType::New();
    overlayReader->SetFileName( overlayFile.c_str() );
    overlayReader->Update();
    
    ImageType::Pointer output = reader->GetOutput();
    OverlayImageType::Pointer overlay = overlayReader->GetOutput();
            
    ImageType::RegionType region = output->GetBufferedRegion();
    typedef itk::ImageRegionIterator<ImageType> Iterator;
    typedef itk::ImageRegionIterator<OverlayImageType> OverlayIterator;
    Iterator biter( output, region );
    OverlayIterator fiter( overlay, region );

    while ( !biter.IsAtEnd() )
      {
        if( fiter.Get()[3] )
         {
         PixelType pix;
         for ( unsigned i = 0; i < 3; i++ )
            {
            pix[i] = fiter.Get()[i];
            }
         biter.Set( pix );
         }
      ++biter;
      ++fiter;
    }

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();

    writer->SetInput( output );
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

