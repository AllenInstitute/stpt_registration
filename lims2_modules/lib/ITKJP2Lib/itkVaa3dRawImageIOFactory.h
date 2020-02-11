/*=========================================================================

  itkVaa3DRawImageIOFactory.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkVaa3dRawImageIOFactory_h
#define __itkVaa3dRawImageIOFactory_h

#include <stddef.h>
#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
	/** \class Vaa3dRawImageIOFactory
	 * \brief Create instances of BMPImageIO objects using an object factory.
	 */
	class ITK_EXPORT Vaa3dRawImageIOFactory : public ObjectFactoryBase
	{
	public:  
		/** Standard class typedefs. */
		typedef Vaa3dRawImageIOFactory   Self;
		typedef ObjectFactoryBase  Superclass;
		typedef SmartPointer<Self>  Pointer;
		typedef SmartPointer<const Self>  ConstPointer;
		
		/** Class methods used to interface with the registered factories. */
		virtual const char* GetITKSourceVersion(void) const;
		virtual const char* GetDescription(void) const;
		
		/** Method for class instantiation. */
		itkFactorylessNewMacro(Self);
		
		/** Run-time type information (and related methods). */
		itkTypeMacro(Vaa3dRawImageIOFactory, ObjectFactoryBase);
		
		/** Register one factory of this type  */
		static void RegisterOneFactory(void)
		{
			Vaa3dRawImageIOFactory::Pointer BMPFactory = Vaa3dRawImageIOFactory::New();
			ObjectFactoryBase::RegisterFactory(BMPFactory);
		}
		
	protected:
		Vaa3dRawImageIOFactory();
		~Vaa3dRawImageIOFactory();
		
	private:
		Vaa3dRawImageIOFactory(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
	};
	
	
} // end namespace itk

#endif
