/*=========================================================================

  itkExtendedImageIOFactory.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkExtendedImageIOFactory.h"
#include "itkJP2ImageIOFactory.h"
#include "itkVaa3dRawImageIOFactory.h"
#include "itkMutexLock.h"
#include "itkMutexLockHolder.h"

namespace itk
{
	
	void
	ExtendedImageIOFactory::RegisterBuiltInFactories()
	{
		static bool firstTime = true;
		
		static SimpleMutexLock mutex;
		{
			// This helper class makes sure the Mutex is unlocked 
			// in the event an exception is thrown.
			MutexLockHolder<SimpleMutexLock> mutexHolder( mutex );
			if( firstTime )
			{
				ObjectFactoryBase::RegisterFactory( JP2ImageIOFactory::New() ); 
				ObjectFactoryBase::RegisterFactory( Vaa3dRawImageIOFactory::New() ); 
				firstTime = false;
			}
		}
		
	}	
		
} // end namespace itk
