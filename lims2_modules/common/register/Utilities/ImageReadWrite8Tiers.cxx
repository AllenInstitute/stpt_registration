/*=========================================================================

  ImageReadWrite8Tiers.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImage.h"
#include "itkExtendedImageIOFactory.h"
#include "itkTextOutput.h"

#include "itkJP2ImageIO.h"

//---------------- code for handling kdu messaging --------------------------//
#include "kdu_messaging.h"

class kdu_stream_message : public kdu_message {
  public: // Member classes
    kdu_stream_message(std::ostream *stream)
      { this->stream = stream; }
    void put_text(const char *string)
      { (*stream) << string; }
    void flush(bool end_of_message=false)
      { stream->flush(); }
  private: // Data
    std::ostream *stream;
  };

static kdu_stream_message cout_message(&std::cout);
static kdu_stream_message cerr_message(&std::cerr);
static kdu_message_formatter pretty_cout(&cout_message);
static kdu_message_formatter pretty_cerr(&cerr_message);

//--------------------------------------------------------------------//


int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile  outputImageFile " << std::endl;
    return -1;
    }

   itk::OutputWindow::SetInstance(itk::TextOutput::New().GetPointer());


  // Custom messaging services
  kdu_customize_warnings(&pretty_cout);
  kdu_customize_errors(&pretty_cerr);

  
  // Need to load in the additional factory which includes jp2 
  itk::ExtendedImageIOFactory::RegisterBuiltInFactories(); 


  typedef itk::RGBPixel<unsigned char>   PixelType;
  const   unsigned int        Dimension = 2;
  typedef itk::Image< PixelType, Dimension >    ImageType;

  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
 
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
 
  const char * inputFilename  = argv[1];
  const char * outputFilename = argv[2];

  reader->SetFileName( inputFilename  );
  writer->SetFileName( outputFilename );

  itk::JP2ImageIO::Pointer io = itk::JP2ImageIO::New();
  if ( ! io->CanWriteFile( outputFilename ) )
    {
    std::cout << outputFilename;
    std::cout << " is not a valid output jp2 file" << std::endl;
    return EXIT_FAILURE;
    }

  const char * param = "Clayers=20 Corder=RPCL Clevels=8 Creversible=yes Ckernels=W5X3 Cuse_precincts=yes Cblk={32,32} Cprecincts={256,256},{256,256},{128,128}";
  io->SetCompressionParameterString( param );

  itk::JP2ImageIO::RatesType rates;
  rates.SetSize( 1 );
  rates[0] = 1.5;
  io->SetRates( rates );

  writer->SetImageIO( io );

  writer->SetInput( reader->GetOutput() );

  try 
    { 
    writer->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return -1;
    } 

  return 0;
}



