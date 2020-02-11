/*=========================================================================

  CropVolume.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

int main( int argc, char *argv[] )
{
  if (argc < 9 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile outputFile ";
    std::cout << "startx starty startz ";
    std::cout << "sizex sizey sizez ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned char, Dimension >      ImageType;

  // Read in input image
  typedef itk::ImageFileReader< ImageType >   ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    }

  // Setup meta information for output image
  ImageType::Pointer input  = reader->GetOutput();
  ImageType::Pointer output = ImageType::New();

  ImageType::RegionType outputRegion;
  for( unsigned int j = 0; j < Dimension; j++ )
    {
    outputRegion.SetIndex( j, atoi( argv[3 + j]) );
    outputRegion.SetSize(  j, atoi( argv[6 + j]) );
    }

  bool ok = outputRegion.Crop( input->GetBufferedRegion() );

  if( !ok )
    {
    std::cout << "Error: Crop box is completely outside input image domain." << std::endl;
    return EXIT_FAILURE;
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
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( output );
  writer->SetFileName( argv[2] );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    }

  return EXIT_SUCCESS;
}


