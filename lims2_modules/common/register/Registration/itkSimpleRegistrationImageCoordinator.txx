/*=========================================================================

  PitkSimpleRegistrationImageCoordinator.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkSimpleRegistrationImageCoordinator_txx
#define _itkSimpleRegistrationImageCoordinator_txx

#include "itkSimpleRegistrationImageCoordinator.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"

namespace itk
{

template< typename TFixedImage, typename TMovingImage >
SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::SimpleRegistrationImageCoordinator()
{
  m_FixedImage = NULL;
  m_FixedImages.resize( 0 );
  m_FixedImageMask = NULL;
  m_FixedMaskImage = NULL;
  m_FixedImageStartingFactors.Fill( 1 );

  m_MovingImage = NULL;
  m_MovingImages.resize( 0 );
  m_MovingImageMask = NULL;
  m_MovingMaskImage = NULL;
  m_MovingImageStartingFactors.Fill( 1 );

  m_Verbose = false;
  
}

template< typename TFixedImage, typename TMovingImage >
void
SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );

}



template< typename TFixedImage, typename TMovingImage >
void
SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::Initialize() throw (ExceptionObject)
{
  // Create pyramid from fixed image
  if ( m_FixedImage.IsNull() )
    {
    itkExceptionMacro( << "FixedImage not set" );
    }

   if ( this->GetNumberOfLevels() )
    {
    typedef Image<float,FixedImageDimension>              FloatImageType;
    typedef RecursiveMultiResolutionPyramidImageFilter<
                                            FixedImageType,
                                            FloatImageType > PyramidType;
    typename PyramidType::Pointer pyramid = PyramidType::New();
  
    pyramid->SetInput( m_FixedImage );
    pyramid->SetNumberOfLevels( this->GetNumberOfLevels() );

    pyramid->SetStartingShrinkFactors( m_FixedImageStartingFactors.GetDataPointer() );
    pyramid->UpdateLargestPossibleRegion();

    m_FixedImages.resize( this->GetNumberOfLevels() );

    for ( unsigned int k = 0; k < this->GetNumberOfLevels(); k++ )
      {
      typedef CastImageFilter<FloatImageType,FixedImageType> CasterType;
      typename CasterType::Pointer caster = CasterType::New();
      caster->SetInput( pyramid->GetOutput( k ) );
      caster->Update();
      m_FixedImages[k] = caster->GetOutput();
      m_FixedImages[k]->DisconnectPipeline();
      }

    }

  // Create pyramid from moving image  
  if ( m_MovingImage.IsNull() )
    {
    itkExceptionMacro( << "MovingImage not set" );
    }

   if ( this->GetNumberOfLevels() )
    {
    typedef Image<float,FixedImageDimension>              FloatImageType;
    typedef RecursiveMultiResolutionPyramidImageFilter<
                                            MovingImageType,
                                            FloatImageType > PyramidType;
    typename PyramidType::Pointer pyramid = PyramidType::New();
  
    pyramid->SetInput( m_MovingImage );
    pyramid->SetNumberOfLevels( this->GetNumberOfLevels() );

    pyramid->SetStartingShrinkFactors( m_MovingImageStartingFactors.GetDataPointer() );
    pyramid->UpdateLargestPossibleRegion();

    m_MovingImages.resize( this->GetNumberOfLevels() );

    for ( unsigned int k = 0; k < this->GetNumberOfLevels(); k++ )
      {
      typedef CastImageFilter<FloatImageType,FixedImageType> CasterType;
      typename CasterType::Pointer caster = CasterType::New();
      caster->SetInput( pyramid->GetOutput( k ) );
      caster->Update();
      m_MovingImages[k] = caster->GetOutput();
      m_MovingImages[k]->DisconnectPipeline();
      }

    }                       
  // Create fixed image mask
  if ( m_FixedMaskImage.IsNotNull() )
    {
    m_FixedImageMask = FixedImageMaskType::New();
    m_FixedImageMask->SetImage( m_FixedMaskImage );
    }

  // Create moving image mask
  if ( m_MovingMaskImage.IsNotNull() )
    {
    m_MovingImageMask = MovingImageMaskType::New();
    m_MovingImageMask->SetImage( m_MovingMaskImage );
    }

  
}

template< typename TFixedImage, typename TMovingImage >
typename SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::FixedImageType *
SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::GetFixedImage( unsigned long idx ) const
{
  if ( idx > this->GetNumberOfLevels() - 1 ) 
    {
    idx = this->GetNumberOfLevels() - 1;
    }
  if ( idx > m_FixedImages.size() - 1 )
    {
    return NULL;
    }
  
  return m_FixedImages[idx];
 
}

template< typename TFixedImage, typename TMovingImage >
typename SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::FixedImageMaskType *
SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::GetFixedImageMask() const
{
  return m_FixedImageMask;
}

template< typename TFixedImage, typename TMovingImage >
typename SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::MovingImageMaskType *
SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::GetMovingImageMask() const
{
  return m_MovingImageMask;
}


template< typename TFixedImage, typename TMovingImage >
typename SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::MovingImageType *
SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
::GetMovingImage( unsigned long idx ) const
{

 if ( idx > this->GetNumberOfLevels() - 1 ) 
    {
    idx = this->GetNumberOfLevels() - 1;
    }
  if ( idx > m_MovingImages.size() - 1 )
    {
    return NULL;
    }
  
  return m_MovingImages[idx];
 
}

} // namespace itk

#endif
