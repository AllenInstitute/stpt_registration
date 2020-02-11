/*=========================================================================

  itkVaa3DRawImageIOFactory.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#include "itkVaa3dRawImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkVaa3dRawImageIO.h"
#include "itkVersion.h"


namespace itk
{
	
    Vaa3dRawImageIOFactory::Vaa3dRawImageIOFactory()
    {
	this->RegisterOverride("itkImageIOBase",
			       "itkVaa3dRawImageIO",
			       "Vaa3dRaw Image IO",
			       1,
			       CreateObjectFunction<Vaa3dRawImageIO>::New());
    }
	
    Vaa3dRawImageIOFactory::~Vaa3dRawImageIOFactory()
    {
    }
	
    const char* 
    Vaa3dRawImageIOFactory::GetITKSourceVersion(void) const
    {
	return ITK_SOURCE_VERSION;
    }
	
    const char* 
    Vaa3dRawImageIOFactory::GetDescription() const
    {
	return "Vaa3dRaw ImageIO Factory, allows the loading of Vaa3dRaw images into Insight";
    }
	
} // end namespace itk

