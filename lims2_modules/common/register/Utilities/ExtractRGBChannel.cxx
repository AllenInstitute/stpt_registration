/*=========================================================================

  ExtractRGBChannel.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include "itkRGBToGrayscalePixelAccessor.h"
#include "itkImageAdaptor.h"
#include "itkCastImageFilter.h"

int main( int argc, char *argv[] )
{
  if (argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputImage outputImage channelIndex";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::RGBPixel< unsigned char > PixelType;
  typedef itk::Image< PixelType, 2 >     ImageType;
  typedef itk::Image< unsigned char, 2 > OutputImageType;
  
  std::string inputFile = argv[1];
  std::string outputFile = argv[2];
  unsigned int channelIndex = atoi( argv[3] );
  
  try
    {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName( inputFile.c_str() );
    reader->Update();

    typedef itk::ImageAdaptor< ImageType, itk::RGBToGrayscalePixelAccessor<unsigned char> > ImageAdaptorType;
    itk::RGBToGrayscalePixelAccessor<unsigned char> accessor;

    if ( channelIndex == 0 )
      {
      accessor.RedOnly();
      }
    else if ( channelIndex == 1 )
      {
      accessor.GreenOnly();
      }
    else if ( channelIndex == 2 )
      {
      accessor.BlueOnly();
      }

    ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SetImage( reader->GetOutput() );
    adaptor->SetPixelAccessor( accessor );

    typedef itk::CastImageFilter< ImageAdaptorType, OutputImageType > CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput( adaptor );

    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();

    writer->SetInput( caster->GetOutput() );
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

