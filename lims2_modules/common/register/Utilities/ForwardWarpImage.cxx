/*=========================================================================

  FowardWarpImage.cxx

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
  if (argc < 7 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile refFile ";
    std::cout << "fieldX fieldY fieldZ ";
    std::cout << "outputFile [defaultValue] ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned short, Dimension >      ImageType;
  typedef itk::Image< float, Dimension > FieldType;

  std::string inputFile = argv[1];
  std::string refFile = argv[2];
  std::string fieldFile[3];
  fieldFile[0] = argv[3];
  fieldFile[1] = argv[4];
  fieldFile[2] = argv[5];
  std::string outputFile = argv[6];

  unsigned char defaultValue = 0;
  if (argc > 7)
    {
    defaultValue = atoi( argv[7] );
    }
    

  ImageType::Pointer input;
  ImageType::Pointer ref;
  ImageType::Pointer output;
  FieldType::Pointer field[3];

  try
    {

    // Read in input image
    typedef itk::ImageFileReader< ImageType >   ReaderType;

    ReaderType::Pointer ir = ReaderType::New();
    ir->SetFileName( inputFile.c_str() );

    ReaderType::Pointer rr = ReaderType::New();
    rr->SetFileName( refFile.c_str() );

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

    rr->Update();
    ref = rr->GetOutput();
    ref->DisconnectPipeline();

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
  ImageType::RegionType outputRegion = ref->GetBufferedRegion();
  output = ImageType::New();
  output->CopyInformation( ref );
  output->SetRegions( outputRegion );
  output->Allocate();
  output->FillBuffer( defaultValue );

  // Populate output image
  typedef itk::ImageRegionConstIterator< ImageType > Iterator;
  typedef itk::ImageRegionConstIterator< FieldType > FieldIterator;
  Iterator iiter( input, input->GetBufferedRegion() );
  iiter.GoToBegin();

  FieldIterator fiter[3];
  for( unsigned int k = 0; k < 3; k++ )
    {
    fiter[k] = FieldIterator( field[k], input->GetBufferedRegion() );
    fiter[k].GoToBegin();
    }

  while( !iiter.IsAtEnd() )
    {

    ImageType::PointType point;
    ImageType::IndexType index;
    for ( unsigned int k = 0; k < 3; k++ )
      {
      point[k] = fiter[k].Get();
      }

    output->TransformPhysicalPointToIndex( point, index );
    if ( outputRegion.IsInside( index ) )
      {
      output->SetPixel( index, iiter.Get() );
      }

    ++iiter;
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


