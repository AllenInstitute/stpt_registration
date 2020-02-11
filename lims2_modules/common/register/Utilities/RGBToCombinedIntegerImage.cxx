/*=========================================================================

  RGBToCombinedIntegerImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include "itkRGBToCombinedIntegerPixelAccessor.h"
#include "itkImageAdaptor.h"
#include "itkCastImageFilter.h"

int main( int argc, char *argv[] )
{
  if (argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputImage outputImage";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::RGBPixel< unsigned char > PixelType;
  typedef itk::Image< PixelType, 2 >     ImageType;
  typedef itk::Image< unsigned short, 2 > OutputImageType;
  
  std::string inputFile = argv[1];
  std::string outputFile = argv[2];
  
  try
    {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName( inputFile.c_str() );
    reader->Update();

    typedef itk::ImageAdaptor< ImageType, itk::RGBToCombinedIntegerPixelAccessor<unsigned short> > ImageAdaptorType;
    itk::RGBToCombinedIntegerPixelAccessor<unsigned short> accessor;

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

