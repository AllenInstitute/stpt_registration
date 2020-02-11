/*=========================================================================

  ConvertVaa3dRaw.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtendedImageIOFactory.h"
#include "itkJP2ImageIO.h"

int main(int argc, char** argv) 
{
    if (argc != 3) {
	std::cout << "Usage: ConvertVaa3dRaw <input file> <output file>" << std::endl;
	return 1;
    }

    typedef itk::Image<float, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<ImageType> WriterType;

    itk::ExtendedImageIOFactory::RegisterBuiltInFactories(); 

    itk::JP2ImageIO io;
    io.SetReduce(2);
    std::cout << "reading " << argv[1] << std::endl;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetImageIO(&io);
    reader->SetFileName(argv[1]);
    reader->Update();

    std::cout << "writing " << argv[2] << std::endl;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(argv[2]);
    writer->SetInput(reader->GetOutput());
    writer->Update();

    return 0;
}
