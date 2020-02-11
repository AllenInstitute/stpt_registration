/*=========================================================================

  itkExtendedImageIOFactory.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#ifndef __itkExtendedImageIOFactory_h
#define __itkExtendedImageIOFactory_h

#include <stddef.h>
#include "itkImageIOFactory.h"

namespace itk
{
	/** \class ExtendedImageIOFactory
	 * \brief Create instances of ImageIO objects using an object factory.
	 */
	class ITK_EXPORT ExtendedImageIOFactory : public ImageIOFactory
	{
	public:  
		/** Standard class typedefs. */
		typedef ExtendedImageIOFactory   Self;
		typedef ImageIOFactory  Superclass;
		typedef SmartPointer<Self>  Pointer;
		typedef SmartPointer<const Self>  ConstPointer;
		
		/** Run-time type information (and related methods). */
		itkTypeMacro(ExtendedImageIOFactory, ImageIOFactory);
		
		/** Register Built-in factories */
		static void RegisterBuiltInFactories();
		
	protected:
		ExtendedImageIOFactory();
		~ExtendedImageIOFactory();
		
	private:
		ExtendedImageIOFactory(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
	};
	
	
} // end namespace itk

#endif
