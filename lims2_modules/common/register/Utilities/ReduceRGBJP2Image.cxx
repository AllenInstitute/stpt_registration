/*=========================================================================
 
 ReduceRGBJP2Image.cxx
 
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


int main( int argc, char *argv[] )
{
	if (argc < 4 )
    {
		std::cout << "Usage: " << argv[0] << " ";
		std::cout << "inputFile outputFile factor ";
		std::cout << "[startx] [starty] [sizex] [sizey]";
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
	
	// Read in input image
	typedef itk::ImageFileReader< ImageType >   ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );
	
	typedef itk::JP2ImageIO IOType;
	IOType::Pointer io = IOType::New();
	
	io->SetReduce( atoi( argv[3] ) );
	reader->SetImageIO( io );
	std::cout << "Reduce: " << io->GetReduce() << std::endl;
	
	ImageType::RegionType outputRegion;
	
	if( argc > 7 )
    {
		for( unsigned int j = 0; j < Dimension; j++ )
		{
			outputRegion.SetIndex( j, atoi( argv[4 + j]) );
			outputRegion.SetSize(  j, atoi( argv[6 + j]) );
		}
		io->SetRegionOfInterest( outputRegion );
		
		std::cout << "ROI: " << io->GetRegionOfInterest() << std::endl;
    }
	
	try
    {
		// Write to file
		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( argv[2] );		
		writer->SetInput( reader->GetOutput() );		
		writer->Update();
    }
	catch( itk::ExceptionObject & excp )
    {
		std::cerr << "Exception:" << std::endl;
		std::cerr << excp << std::endl;
    }
	catch( ...)
    {
		std::cerr << "Unknown exception." << std::endl;
    }
	
	return EXIT_SUCCESS;
}
