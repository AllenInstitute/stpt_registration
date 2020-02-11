/*=========================================================================

  CropImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itksys/String.hxx"

#include <string>


template <typename ImageType > 
int Pipeline( 
 const char * infile, 
 const char * outfile,
 const typename ImageType::IndexType & index,
 const typename ImageType::SizeType & size )
{

  typedef itk::ImageFileReader< ImageType >   ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( infile );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image " << infile << std::endl;
    std::cerr << excp << std::endl;
    return -1;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    return -1;
    }

  // Setup meta information for output image
  typename ImageType::Pointer input  = reader->GetOutput();
  typename ImageType::Pointer output = ImageType::New();

  typename ImageType::RegionType outputRegion;
  outputRegion.SetIndex( index );
  outputRegion.SetSize( size );
 
  bool ok = outputRegion.Crop( input->GetBufferedRegion() );

  if( !ok )
    {
    std::cout << "Error: Crop box is completely outside input image domain." << std::endl;
    return -1;
    }

  output->CopyInformation( input );
  output->SetRegions( outputRegion );
  output->Allocate();

  // Populate output image
  typedef itk::ImageRegionConstIterator< ImageType > InputIterator;
  typedef itk::ImageRegionIterator< ImageType > OutputIterator;

  InputIterator iiter( input, outputRegion );
  OutputIterator oiter( output, outputRegion );

  while( !oiter.IsAtEnd() )
    {
    oiter.Set( iiter.Get() );
    ++iiter;
    ++oiter;
    }


  // Write to file
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();

  writer->SetInput( output );
  writer->SetFileName( outfile );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the image " << outfile << std::endl;
    std::cerr << excp << std::endl;
    return -1;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    return -1;
    }
   return 0;
}


int main( int argc, char *argv[] )
{
  if (argc < 7 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile outputFile ";
    std::cout << "startx starty ";
    std::cout << "sizex sizey ";
    std::cout << "[format] ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string infile = argv[1];
  std::string outfile = argv[2];

  const unsigned int Dimension = 2;
  typedef itk::Image< unsigned char, Dimension >      UCHARImageType;
  typedef itk::Image< unsigned short, Dimension >     USHORTImageType;
  typedef itk::Image< unsigned int,   Dimension >     UINTImageType;
  typedef itk::RGBPixel< unsigned char >              RGBPixelType;
  typedef itk::Image< RGBPixelType, Dimension >       RGBImageType;

  UCHARImageType::IndexType index;
  UCHARImageType::SizeType  size;
  for( unsigned int j = 0; j < Dimension; j++ )
    {
    index[j] = atoi( argv[3 + j] );
    size[j]  = atoi( argv[3 + Dimension + j] );
    }

  std::string format = "uchar";
  if (argc > 7 )
    {
    format = argv[7];
    }

  int err;

  if( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0 )
    {
    err = Pipeline<UCHARImageType>( infile.c_str(), outfile.c_str(), index, size );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 )
    {
    err = Pipeline<USHORTImageType>( infile.c_str(), outfile.c_str(), index, size );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 )
    {
    err = Pipeline<UINTImageType>( infile.c_str(), outfile.c_str(), index, size );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "rgb" ) == 0 )
    {
    err = Pipeline<RGBImageType>( infile.c_str(), outfile.c_str(), index, size );
    if ( err ) return EXIT_FAILURE;
    }
  else 
    {
    std::cerr << "format type " << format << "not supported" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


