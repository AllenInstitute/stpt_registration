/*=========================================================================

  itkVaa3dRawImageIO.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#include "itkVaa3dRawImageIO.h"
#include <fstream>

#include "itkByteSwapper.h"

std::string VAA3D_HEADER_FORMAT_KEY = "raw_image_stack_by_hpeng";
unsigned int VAA3D_HEADER_SIZE = VAA3D_HEADER_FORMAT_KEY.size() + 2 + 4*4 + 1;

namespace itk
{
    Vaa3dRawImageIO
    ::Vaa3dRawImageIO()
    {
    }

    Vaa3dRawImageIO
    ::~Vaa3dRawImageIO()
    {
    }

    bool Vaa3dRawImageIO
    ::CanReadFile( const char* filename )
    {
	std::string fname = filename;
	if ( fname == "" )
	{
	    itkDebugMacro(<< "No filename specified.");
	}
		
	std::ifstream f(filename, std::ios::binary);
	if (!f) 
	{
	    itkExceptionMacro(<< "Could not open file: " << filename);
	}
		
	unsigned int headerLength = VAA3D_HEADER_FORMAT_KEY.size();
	char * headerKey = new char [headerLength+1];
	if (!headerKey)
	{
	    itkExceptionMacro(<< "Failed to allocate memory.");
	}
		
	f.read(headerKey, headerLength);
	if (!f)
	{
	    itkExceptionMacro(<< "File unrecognized or corrupted.");
	}
		
	headerKey[headerLength] = '\0';

	if (std::string(headerKey) != VAA3D_HEADER_FORMAT_KEY)
	{
	    return false;
	}

	f.close();

	return true;
    }

    void Vaa3dRawImageIO
    ::ReadImageInformation()
    {
	const char* filename = this->GetFileName();

	std::ifstream f(filename, std::ios::binary);
	if (!f) 
	{
	    itkExceptionMacro(<< "Could not open file: " << filename);
	}

	// jump past the header key
	f.seekg(VAA3D_HEADER_FORMAT_KEY.size(), f.beg);

	char fileByteOrderCode;
	f.read(&fileByteOrderCode, 1);

	if (!f) 
	{
	    itkExceptionMacro(<<"Could not read endian-ness code.");
	}
	else if (fileByteOrderCode == 'B')
	{
	    m_FileByteOrder = itk::ImageIOBase::BigEndian;
	}
	else if (fileByteOrderCode == 'L')
	{
	    m_FileByteOrder = itk::ImageIOBase::LittleEndian;
	}
	else
	{
	    itkExceptionMacro(<<"Endian-ness code not supported: " << fileByteOrderCode);
	}

	typedef itk::ByteSwapper<short> ShortSwapper;
	typedef itk::ByteSwapper<unsigned int> IntSwapper;
	typedef itk::ByteSwapper<long> LongSwapper;

	short dcode = 0;
	f.read(reinterpret_cast<char*>(&dcode), 2); 

	if (!f) 
	{
	    itkExceptionMacro(<<"Could not read raw data type.");
	}
	else if (m_FileByteOrder == itk::ImageIOBase::LittleEndian)
	{
	    ShortSwapper::SwapFromSystemToLittleEndian(&dcode);
	}
	else
	{
	    ShortSwapper::SwapFromSystemToBigEndian(&dcode);
	}

	switch (dcode)
	{
	case 1:
	    m_ComponentType = UCHAR;
	    break;
			
	case 2:
	    m_ComponentType = USHORT;
	    break;
			
	case 4:
	    m_ComponentType = FLOAT;
	    break;
			
	default:
	{
	    itkExceptionMacro(<<"Unrecognized data type code: " << dcode << ". The file type is incorrect or this code is not supported in this version.\n");
	    break;
	}
	}

		
	unsigned int dims[4] = { 0, 0, 0, 0 };
	f.read(reinterpret_cast<char*>(&dims[0]), sizeof(unsigned int)*4);

	if (!f)
	{
	    itkExceptionMacro(<< "Incorrect number of bytes read while reading image dimensions>");
	}

	if (m_FileByteOrder == itk::ImageIOBase::LittleEndian) 
	{
	    IntSwapper::SwapRangeFromSystemToLittleEndian(dims, 4);
	}
	else
	{
	    IntSwapper::SwapRangeFromSystemToBigEndian(dims, 4);
	}

	unsigned int numDims = 0;
	unsigned int numComps = 0;
	unsigned int c = 0;
	for (c = 0; c < 4; c++) 
	{
	    if (dims[c] == 0) 
		break;
	}

	numDims = c-1;
	numComps = dims[c-1];

	this->SetNumberOfComponents(numComps);
	this->SetNumberOfDimensions(numDims);
	for (unsigned int i = 0; i < numDims; i++) 
	{
	    m_Dimensions[i] = dims[i];
	}
    }

    void Vaa3dRawImageIO
    ::Read(void* buffer)
    {
	const char* filename = this->GetFileName();

	std::ifstream f(filename, std::ios::binary);
	if (!f) 
	{
	    itkExceptionMacro(<< "Could not open file: " << filename);
	}

	// get file length
	f.seekg(0, f.end);
	unsigned long fileSize = f.tellg();

	// seek to the beginning of the raw data
	f.seekg(VAA3D_HEADER_SIZE, f.beg);

	unsigned long bytesRemaining = fileSize - VAA3D_HEADER_SIZE;

	unsigned int componentSize = 0;
	switch (m_ComponentType)
	{
	case UCHAR:
	    componentSize = 1;
	    break;
			
	case USHORT:
	    componentSize = 2;
	    break;
			
	case FLOAT:
	    componentSize = 4;
	    break;
	default:
	{
	    itkExceptionMacro(<<"Unrecognized data type: " << m_ComponentType)
		break;
	}
	}

	unsigned long expectedBytesRemaining = componentSize * this->GetNumberOfComponents();
	for (unsigned int i = 0; i < this->GetNumberOfDimensions(); i++)
	    expectedBytesRemaining *= this->GetDimensions(i);

	if (expectedBytesRemaining != bytesRemaining) 
	{
	    itkExceptionMacro(<< "Incorrect file size based on dimensions in header.");
	}

	f.read(reinterpret_cast<char*>(buffer), bytesRemaining);

	if (!f) 
	{
	    itkExceptionMacro(<< "Incorrect number of bytes read while reading image dimensions.");
	}


	unsigned long length = expectedBytesRemaining / componentSize;

	switch (m_ComponentType) 
	{
	case USHORT:
	    if (m_FileByteOrder == itk::ImageIOBase::LittleEndian) 
		itk::ByteSwapper<unsigned short>::SwapRangeFromSystemToLittleEndian(reinterpret_cast<unsigned short*>(buffer), length);
	    else
		itk::ByteSwapper<unsigned short>::SwapRangeFromSystemToBigEndian(reinterpret_cast<unsigned short*>(buffer), length);
	    break;
			
	case FLOAT:
	    if (m_FileByteOrder == itk::ImageIOBase::LittleEndian) 
		itk::ByteSwapper<float>::SwapRangeFromSystemToLittleEndian(reinterpret_cast<float*>(buffer), length);
	    else
		itk::ByteSwapper<float>::SwapRangeFromSystemToBigEndian(reinterpret_cast<float*>(buffer), length);
	    break;
	default:
	    break;
	}
    }

    bool Vaa3dRawImageIO
    ::CanWriteFile( const char* filename )
    {
	std::string fname = filename;
	if ( fname == "" )
	{
	    itkDebugMacro(<< "No filename specified.");
	}

	bool extensionFound = false;
	std::string::size_type pos = fname.rfind(".v3draw");
	if ((pos != std::string::npos)
	    && (pos == fname.length() - 7))
	{
	    extensionFound = true;
	}

	return extensionFound;
    }

    void Vaa3dRawImageIO
    ::WriteImageInformation()
    {
    }
	

    void Vaa3dRawImageIO
    ::Write(const void* buffer)
    {
	const char* filename = this->GetFileName();
	
	std::ofstream f(filename, std::ios::binary | std::ios::out);

	// write the header
	long lenkey = VAA3D_HEADER_FORMAT_KEY.size();
	f.write(VAA3D_HEADER_FORMAT_KEY.c_str(), lenkey); 
/*
	if (nwrite != lenkey)
	{
	    itkExceptionMacro(<< "failed to write format key: ");
	}
*/

	// assuming that this will be set to system endianness
	char endianCode = 0;
	if (m_ByteOrder == itk::ImageIOBase::BigEndian)
	    endianCode = 'B';
	else 
	    endianCode = 'L';
	
	if (endianCode != 'B' && endianCode != 'L')
	{
	    itkExceptionMacro(<< "invalid endianness code: " << endianCode);
	}

	f.write(&endianCode, 1);
/*
	if (nwrite != 1)
	{
	    itkExceptionMacro(<< "failed to write endianness code: " << endianCode);
	}
*/
	
	short int dcode	= 0;
	switch (m_ComponentType)
	{
	case UCHAR:
	    dcode=1;
	    break;
	case USHORT:
	    dcode=2;
	    break;
	case FLOAT:
	    dcode=4;
	    break;
	default:
	    itkExceptionMacro(<< "unsupported component type: " << m_ComponentType );
	    break;
	}
	
	f.write(reinterpret_cast<const char*>(&dcode), 2);
/*
	if (nwrite != 2)
	{
	    itkExceptionMacro(<< "failed to write data type code");
	}
*/
	long unitSize = dcode;

	unsigned int size[4] = { 1, 1, 1, 1 };
	for (unsigned int i = 0 ; i < this->GetNumberOfDimensions(); i++) 
	{
	    size[i] = this->GetDimensions(i);
	    std::cout << "size " << i << ": " << size[i] << std::endl;
	}

	f.write(reinterpret_cast<const char*>(size), 4 * sizeof(unsigned int));

	long totalUnit = 1;
	for (unsigned int i = 0; i < 4; i++) 
	{
	    totalUnit *= size[i];
	}
	
	std::cout << "unit size: " << unitSize << std::endl;
	std::cout << "totalUnit: " << totalUnit << std::endl;

	f.write(reinterpret_cast<const char*>(buffer), unitSize * totalUnit);
/*
	if (nwrite != unitSize * totalUnit)
	{
	    itkExceptionMacro(<< "failed to write buffer");
	}
*/	
	f.close();
    }

    void Vaa3dRawImageIO
    ::PrintSelf(std::ostream& os, Indent indent) const
    {
		
    }
}


