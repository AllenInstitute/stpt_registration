/*=========================================================================

  ImageReadImageSeriesWrite.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"

#include "itksys/SystemTools.hxx"
#include <string>

template <class PixelType>
class Pipeline 
{
public:
    typedef itk::Image< PixelType, 2 > ImageType;
    typedef itk::Image< PixelType, 3 > VolumeType;

    int static Execute( const char *   inputFile,
                 const char *   outputPrefix,
                 const char *   outputExtension,
                 unsigned int   numberOfPlaces )
      {
      typedef itk::ImageFileReader< VolumeType > ReaderType;
      typedef itk::ImageSeriesWriter< VolumeType, ImageType > WriterType;

      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( inputFile );

      try
        {
        reader->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cerr << "Exception thrown while reading the image " << inputFile << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
        }
      catch( ... )
        {
        std::cerr << "Caught unknown exception" << std::endl;
        return EXIT_FAILURE;
        }

      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( reader->GetOutput() );
      
      typedef itk::NumericSeriesFileNames NameGeneratorType;
      NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

      std::string format = outputPrefix;
      format += "%0";
      char buffer[10];
      sprintf( buffer, "%d", numberOfPlaces );
      format += buffer;
      format += "d.";
      format += outputExtension;

      nameGenerator->SetSeriesFormat( format.c_str() );

      typename VolumeType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
      nameGenerator->SetStartIndex( region.GetIndex()[2] );
      nameGenerator->SetEndIndex( region.GetIndex()[2] + region.GetSize()[2] - 1 );
      nameGenerator->SetIncrementIndex( 1 );

      writer->SetFileNames( nameGenerator->GetFileNames() );
      
      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cerr << "Exception thrown while reading the image " << inputFile << std::endl;
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

};


int main( int argc, char *argv[] )
{
  if (argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile outputPrefix ";
    std::cout << "outputExtension ";
    std::cout << "[pixelType] [numberOfPlaces]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string infile = argv[1];
  std::string outputPrefix = argv[2];
  std::string outputExtension = argv[3];
  std::string pixelType = "uchar";
  unsigned int numberOfPlaces = 3;

  if ( argc > 4 )
    {
    pixelType = argv[4];
    }
    
  if ( argc > 5 )
    {
    numberOfPlaces = atoi( argv[5] );
    }

  int err;

  if( itksys::SystemTools::Strucmp( pixelType.c_str(), "uchar" ) == 0 )
    {
    err = Pipeline<unsigned char>::Execute( infile.c_str(), outputPrefix.c_str(), 
                                            outputExtension.c_str(), numberOfPlaces );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( pixelType.c_str(), "ushort" ) == 0 )
    {
    err = Pipeline<unsigned short>::Execute( infile.c_str(), outputPrefix.c_str(), 
                                             outputExtension.c_str(), numberOfPlaces );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( pixelType.c_str(), "uint" ) == 0 )
    {
    err = Pipeline<unsigned int>::Execute( infile.c_str(), outputPrefix.c_str(), 
                                           outputExtension.c_str(), numberOfPlaces );
    if ( err ) return EXIT_FAILURE;
    }
  else if( itksys::SystemTools::Strucmp( pixelType.c_str(), "float" ) == 0 )
    {
    err = Pipeline<float>::Execute( infile.c_str(), outputPrefix.c_str(), 
                                    outputExtension.c_str(), numberOfPlaces );
    if ( err ) return EXIT_FAILURE;
    }
  else 
    {
    std::cerr << "format type " << pixelType << "not supported" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


