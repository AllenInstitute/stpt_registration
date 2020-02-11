/*=========================================================================

  itkJP2ImageIOFactory.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#ifndef __itkJP2ImageIOFactory_h
#define __itkJP2ImageIOFactory_h

#include <stddef.h>
#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
	/** \class JP2ImageIOFactory
	 * \brief Create instances of BMPImageIO objects using an object factory.
	 */
	class ITK_EXPORT JP2ImageIOFactory : public ObjectFactoryBase
	{
	public:  
		/** Standard class typedefs. */
		typedef JP2ImageIOFactory   Self;
		typedef ObjectFactoryBase  Superclass;
		typedef SmartPointer<Self>  Pointer;
		typedef SmartPointer<const Self>  ConstPointer;
		
		/** Class methods used to interface with the registered factories. */
		virtual const char* GetITKSourceVersion(void) const;
		virtual const char* GetDescription(void) const;
		
		/** Method for class instantiation. */
		itkFactorylessNewMacro(Self);
		
		/** Run-time type information (and related methods). */
		itkTypeMacro(JP2ImageIOFactory, ObjectFactoryBase);
		
		/** Register one factory of this type  */
		static void RegisterOneFactory(void)
		{
			JP2ImageIOFactory::Pointer BMPFactory = JP2ImageIOFactory::New();
			ObjectFactoryBase::RegisterFactory(BMPFactory);
		}
		
	protected:
		JP2ImageIOFactory();
		~JP2ImageIOFactory();
		
	private:
		JP2ImageIOFactory(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
	};
	
	
} // end namespace itk

#endif
