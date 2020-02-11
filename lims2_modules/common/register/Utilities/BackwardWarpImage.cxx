/*=========================================================================

  BackwardWarpImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <string>

int main( int argc, char *argv[] )
{
  if (argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile ";
    std::cout << "fieldX fieldY fieldZ ";
    std::cout << "outputFile [defaultValue] ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned short, Dimension >      ImageType;
  typedef itk::Image< float, Dimension > FieldType;

  std::string inputFile = argv[1];
  std::string fieldFile[3];
  fieldFile[0] = argv[2];
  fieldFile[1] = argv[3];
  fieldFile[2] = argv[4];
  std::string outputFile = argv[5];

  unsigned char defaultValue = 0;
  if (argc > 6)
    {
    defaultValue = atoi( argv[6] );
    }
    

  ImageType::Pointer input;
  ImageType::Pointer output;
  FieldType::Pointer field[3];

  try
    {

    // Read in input image
    typedef itk::ImageFileReader< ImageType >   ReaderType;

    ReaderType::Pointer ir = ReaderType::New();
    ir->SetFileName( inputFile.c_str() );

    typedef itk::ImageFileReader<FieldType> FieldReaderType;
    FieldReaderType::Pointer fr[3];
    for ( unsigned int k = 0; k < 3; k++ )
      { 
      fr[k] = FieldReaderType::New();
      fr[k]->SetFileName( fieldFile[k].c_str() );
      }

    ir->Update();
    input = ir->GetOutput();
    input->DisconnectPipeline();

    for ( unsigned int k = 0; k < 3; k++ )
      {
      fr[k]->Update();
      field[k] = fr[k]->GetOutput();
      field[k]->DisconnectPipeline();
      }
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    }

  // Setup meta information for output image
  ImageType::RegionType outputRegion = field[0]->GetBufferedRegion();
  output = ImageType::New();
  output->CopyInformation( field[0] );
  output->SetRegions( outputRegion );
  output->Allocate();
  output->FillBuffer( defaultValue );

  ImageType::RegionType inputRegion = input->GetBufferedRegion();

  // Populate output image
  typedef itk::ImageRegionIterator< ImageType > Iterator;
  typedef itk::ImageRegionConstIterator< FieldType > FieldIterator;
  Iterator oiter( output, outputRegion );
  oiter.GoToBegin();

  FieldIterator fiter[3];
  for( unsigned int k = 0; k < 3; k++ )
    {
    fiter[k] = FieldIterator( field[k], outputRegion );
    fiter[k].GoToBegin();
    }

  while( !oiter.IsAtEnd() )
    {


    ImageType::PointType point;
    ImageType::IndexType index;
    for ( unsigned int k = 0; k < 3; k++ )
      {
      point[k] = fiter[k].Get();
      }

    input->TransformPhysicalPointToIndex( point, index );
    if ( inputRegion.IsInside( index ) )
      {
      oiter.Set( input->GetPixel( index ) );
      }

    ++oiter;
    for ( unsigned int k = 0; k < 3; k++ )
      {
      ++fiter[k];
      }
    }


  // Write to file
  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( output );
  writer->SetFileName( outputFile );

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


