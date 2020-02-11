/*=========================================================================

  itkJP2ImageIOFactory.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkJP2ImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkJP2ImageIO.h"
#include "itkVersion.h"


namespace itk
{
	
	JP2ImageIOFactory::JP2ImageIOFactory()
	{
		this->RegisterOverride("itkImageIOBase",
							   "itkJP2ImageIO",
							   "Kakadu JPEG2000 Image IO",
							   1,
							   CreateObjectFunction<JP2ImageIO>::New());
	}
	
	JP2ImageIOFactory::~JP2ImageIOFactory()
	{
	}
	
	const char* 
	JP2ImageIOFactory::GetITKSourceVersion(void) const
	{
		return ITK_SOURCE_VERSION;
	}
	
	const char* 
	JP2ImageIOFactory::GetDescription() const
	{
		return "Kakadu JPEG2000 ImageIO Factory, allows the loading of JP2 images into Insight";
	}
	
} // end namespace itk

