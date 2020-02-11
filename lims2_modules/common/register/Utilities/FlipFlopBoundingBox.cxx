/*=========================================================================
 
 FlipFlopBoundingBox.cxx

 Copyright (c) Allen Institute for Brain Science. All rights reserved.
 
 =========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkExtendedImageIOFactory.h"
#include "itkJP2ImageIO.h"
#include "itkRGBPixel.h"
#include "itkExtractImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkPasteImageFilter.h"

int main( int argc, char *argv[] )
{
	if (argc < 9 )
    {
		std::cout << "Usage: " << argv[0] << " ";
		std::cout << "inputFile outputFile";
		std::cout << "x y width height vflip hflip";
		std::cout << std::endl;
		return EXIT_FAILURE;
    }
	
	// Include additional IO libraries (e.g. kakadu)
	itk::ExtendedImageIOFactory::RegisterBuiltInFactories(); 
	
	const unsigned int Dimension = 2;
	typedef unsigned char ComponentType;
	typedef itk::RGBPixel< ComponentType > PixelType;
	//typedef ComponentType PixelType;
	typedef itk::Image< PixelType , Dimension >      ImageType;

  int x = atoi( argv[3] );
  int y = atoi( argv[4] );
  int width = atoi( argv[5] );
  int height = atoi( argv[6] );
  bool flip = atoi( argv[7] );
  bool flop = atoi( argv[8] );	
	
	// Read in input image
	typedef itk::ImageFileReader< ImageType >   ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );

  // Extract region of interest
  typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractorType;
  ExtractorType::Pointer extractor = ExtractorType::New();
  extractor->SetInput( reader->GetOutput() );

  ImageType::RegionType region;
  region.SetIndex( 0, x );
  region.SetIndex( 1, y );
  region.SetSize( 0, width );
  region.SetSize( 1, height );

  extractor->SetExtractionRegion( region );

  // Flip/flop region of interest
  typedef itk::FlipImageFilter< ImageType > FlipperType;
  FlipperType::Pointer flipper = FlipperType::New();
  flipper->SetInput( extractor->GetOutput() );

  FlipperType::FlipAxesArrayType faa;
  faa[0] = flop;
  faa[1] = flip;

  flipper->SetFlipAxes( faa );

	try
    {
    flipper->Update();

    ImageType::Pointer image = reader->GetOutput();
    image->DisconnectPipeline();

    ImageType::Pointer roi = flipper->GetOutput();
    roi->DisconnectPipeline();

    typedef itk::PasteImageFilter< ImageType > CopierType;
    CopierType::Pointer copier = CopierType::New();

    copier->SetSourceImage( roi );
    copier->SetSourceRegion( roi->GetBufferedRegion() );

    copier->SetInput( image );
    copier->SetDestinationIndex( region.GetIndex() );

    copier->Update();


    itk::JP2ImageIO::Pointer io = itk::JP2ImageIO::New();
    if ( ! io->CanWriteFile( argv[2] ) )
      {
      std::cout << argv[2];
      std::cout << " is not a valid output jp2 file" << std::endl;
      return EXIT_FAILURE;
      }

    io->SetKeepColorChannelsSeparate( true );

    itk::JP2ImageIO::RatesType rates;
    rates.SetSize( 1 );
    rates[0] = 1.5;
    io->SetRates( rates );


		// Write to file
		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( argv[2] );		
		writer->SetInput( copier->GetOutput() );		
    writer->SetImageIO( io );
		writer->Update();

    }
	catch( itk::ExceptionObject & excp )
    {
		std::cerr << "Exception:" << std::endl;
		std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
	catch( ...)
    {
		std::cerr << "Unknown exception." << std::endl;
    return EXIT_FAILURE;
    }

}