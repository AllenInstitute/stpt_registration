/*=========================================================================

  ImageSeriesReadWrite.cxx
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.
  
=========================================================================*/

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"

template <class TImage>
int Pipeline(
const char * inputPrefix,
const char * inputExtension,
unsigned int first,
unsigned int last,
const char * outputFilename,
unsigned int numPlaces
)
{

  const unsigned int Dimension = TImage::ImageDimension;
  typedef TImage  ImageType;
  typedef itk::ImageSeriesReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter<   ImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  char buffer[10];
  sprintf(buffer, "%%0%dd", numPlaces );

  std::string pattern = inputPrefix;
  pattern += buffer;
  pattern += inputExtension;

  typedef itk::NumericSeriesFileNames    NameGeneratorType;

  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
   
  nameGenerator->SetStartIndex( first );
  nameGenerator->SetEndIndex( last );
  nameGenerator->SetIncrementIndex( 1 );

  nameGenerator->SetSeriesFormat( pattern.c_str() );
 
  reader->SetFileNames( nameGenerator->GetFileNames()  );

  writer->SetFileName( outputFilename );

  try
    { 
    reader->Update();
    writer->SetInput( reader->GetOutput() );
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error reading the series " << std::endl;
    std::cerr << excp << std::endl;
    return -1;
    }


  return 0;
}



int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "inputPrefix inputExtension first last outputImageFile [numPlaces] [type]" << std::endl;
    return -1;
    }

  std::string inputPrefix = argv[1];
  std::string inputExtension = argv[2];
  unsigned int first = atoi( argv[3] );
  unsigned int last = atoi( argv[4] );
  std::string outputImageFile = argv[5];
  unsigned int numPlaces = 3;
  std::string format = "UCHAR";

  if ( argc > 6 )
    {
    numPlaces = atoi( argv[6] );
    }

  if ( argc > 7 )
    {
    format = argv[7];
    }

  typedef itk::Image< unsigned char, 3 >      UCHARVolumeType;
  typedef itk::Image< unsigned short, 3 >     USHORTVolumeType;
  typedef itk::Image< unsigned int,  3 >      UINTVolumeType;
  typedef itk::Image< float, 3>               FLOATVolumeType;

  int err;

  if( itksys::SystemTools::Strucmp( format.c_str(), "uchar" ) == 0 )
    {
    err = Pipeline<UCHARVolumeType>( inputPrefix.c_str(), inputExtension.c_str(),
                                     first, last, outputImageFile.c_str(), numPlaces );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "ushort" ) == 0 )
    {
    err = Pipeline<USHORTVolumeType>( inputPrefix.c_str(), inputExtension.c_str(),
                                     first, last, outputImageFile.c_str(), numPlaces );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "uint" ) == 0 )
    {
    err = Pipeline<UINTVolumeType>( inputPrefix.c_str(), inputExtension.c_str(),
                                     first, last, outputImageFile.c_str(), numPlaces );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( format.c_str(), "float" ) == 0 )
    {
    err = Pipeline<FLOATVolumeType>( inputPrefix.c_str(), inputExtension.c_str(),
                                     first, last, outputImageFile.c_str(), numPlaces );
    if ( err ) return EXIT_FAILURE;
    }

  else 
    {
    std::cerr << "format: " << format << " not supported" << std::endl;
    return EXIT_FAILURE;
    }

  return 0;
}



