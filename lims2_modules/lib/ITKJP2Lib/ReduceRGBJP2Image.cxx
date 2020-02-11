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

template <typename ComponentType>
int Pipeline(
  const char * infile,
  const char * outfile,
  signed int   factor,
  bool         useRegionOfInterest,
  signed int   startx,
  signed int   starty,
  signed int   sizex,
  signed int   sizey ) 
{
	const unsigned int Dimension = 2;
	typedef itk::RGBPixel< ComponentType > PixelType;
	typedef itk::Image< PixelType , Dimension >      ImageType;
	
	// Read in input image
	typedef itk::ImageFileReader< ImageType >   ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( infile );
	
	typedef itk::JP2ImageIO IOType;
	typename IOType::Pointer io = IOType::New();
	
	io->SetReduce( factor );
	reader->SetImageIO( io );
	std::cout << "Reduce: " << io->GetReduce() << std::endl;
	
	typename ImageType::RegionType outputRegion;
	
	if( useRegionOfInterest )
    {
        outputRegion.SetIndex( 0, startx );
        outputRegion.SetIndex( 1, starty );
        outputRegion.SetSize( 0, sizex );
        outputRegion.SetSize( 1, sizey );
		io->SetRegionOfInterest( outputRegion );		
		std::cout << "ROI: " << io->GetRegionOfInterest() << std::endl;
    }
	
	try
    {
		// Write to file
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( outfile );		
		writer->SetInput( reader->GetOutput() );		
		writer->Update();
    }
	catch( itk::ExceptionObject & excp )
    {
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
    }
	catch( ...)
    {
        std::cout << "caught unknown exception" << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
    
}

int main( int argc, char *argv[] )
{
	if (argc < 4 )
    {
		std::cout << "Usage: " << argv[0] << " ";
		std::cout << "inputFile outputFile factor ";
		std::cout << "[startx] [starty] [sizex] [sizey] [format]";
		std::cout << std::endl;
		return EXIT_FAILURE;
    }
    
    std::string infile = argv[1];
    std::string outfile = argv[2];
    
    unsigned int factor = atoi( argv[3] );
    
    bool useRegionOfInterest = false;
    signed int startx = 0;
    signed int starty = 0;
    signed int sizex = 0;
    signed int sizey = 0;
    
    if ( argc > 7 )
        {
        useRegionOfInterest = true;
        startx = atoi( argv[4] );
        starty = atoi( argv[5] );
        sizex  = atoi( argv[6] );
        sizey  = atoi( argv[7] );
        }
        
    std::string format = "uchar";
    if ( argc > 8 )
        {
        format = argv[8];
        }
    
	
	// Include additional IO libraries (e.g. kakadu)
	itk::ExtendedImageIOFactory::RegisterBuiltInFactories(); 
    
	typedef unsigned char UCHARComponentType;
    typedef unsigned short USHORTComponentType;
    
    int err;
    
    if ( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0  )
        {
        err = Pipeline<UCHARComponentType>( infile.c_str(), outfile.c_str(), factor, 
                                            useRegionOfInterest, startx, starty, sizex, sizey );
        if (err) return EXIT_FAILURE;                                    
        }
    else if ( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 )
        {
        err = Pipeline<USHORTComponentType>( infile.c_str(), outfile.c_str(), factor, 
                                            useRegionOfInterest, startx, starty, sizex, sizey );
        if (err) return EXIT_FAILURE;         
        }
    else 
        {
        std::cerr << "format: " << format << " not supported" << std::endl;
        return EXIT_FAILURE;
    }
    
	return EXIT_SUCCESS;
}
