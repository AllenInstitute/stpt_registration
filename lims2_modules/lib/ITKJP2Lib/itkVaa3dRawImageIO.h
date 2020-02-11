/*=========================================================================

  itkVaa3DRawImageIO.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkVaa3dRawImageIO_h
#define __itkVaa3dRawImageIO_h

#include <stddef.h>
#include "itkImageIOBase.h"

namespace itk
{
    class ITK_EXPORT Vaa3dRawImageIO : public ImageIOBase
    {
    public:
	typedef Vaa3dRawImageIO    Self;
	typedef ImageIOBase        Superclass;
	typedef SmartPointer<Self> Pointer;

	itkNewMacro(Self);
		
	itkTypeMacro(Vaa3dRawImageIO, Superclass);

	virtual bool CanReadFile(const char*);
	virtual void ReadImageInformation();
	virtual void Read(void* buffer);
		
	virtual bool CanWriteFile(const char*);
	virtual void WriteImageInformation();
	virtual void Write(const void* buffer);
		
	Vaa3dRawImageIO();
	~Vaa3dRawImageIO();
	void PrintSelf(std::ostream& os, Indent indent) const;

	ImageIOBase::ByteOrder m_FileByteOrder;

    private:
	Vaa3dRawImageIO(const Self&);
	void operator=(const Self&);
    };
}
#endif
