/*=========================================================================

  ImageReadWrite.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImage.h"
#include "itkExtendedImageIOFactory.h"
#include "itkTextOutput.h"

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
		std::cerr << argv[0] << " inputImageFile outputImageFile " << std::endl;
		return -1;
    }
	
	itk::OutputWindow::SetInstance(itk::TextOutput::New().GetPointer());
	
	
	// Custom messaging services
	kdu_customize_warnings(&pretty_cout);
	kdu_customize_errors(&pretty_cerr);
	
	
	// Need to load in the additional factory which includes jp2 
	itk::ExtendedImageIOFactory::RegisterBuiltInFactories(); 
	
	typedef unsigned short ComponentType;
	typedef itk::RGBPixel<ComponentType>   PixelType;
	// typedef ComponentType       PixelType;
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
