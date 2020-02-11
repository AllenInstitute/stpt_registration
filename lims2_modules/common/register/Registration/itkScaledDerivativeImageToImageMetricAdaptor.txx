/*=========================================================================

  itkScaledDerivativeImageToImageMetricAdaptor.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkScaledDerivativeImageToImageMetricAdaptor_txx
#define __itkScaledDerivativeImageToImageMetricAdaptor_txx

#include "itkScaledDerivativeImageToImageMetricAdaptor.h"

namespace itk
{

template <typename TFixedImage, typename TMovingImage, typename TMetric>
ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::ScaledDerivativeImageToImageMetricAdaptor()
{
  m_Metric = NULL;
  m_ScalesInitialized = false;
}


template <typename TFixedImage, typename TMovingImage, typename TMetric>
void
ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::GetDerivative(
const ParametersType & parameters,
DerivativeType & derivative ) const
{

  if ( !m_Metric )
    {
    itkExceptionMacro(<<"Metric is not present");
    }

  // Delegate computation to the underlying Metric
  m_Metric->GetDerivative( parameters, derivative );

  const unsigned int spaceDimension = 
    vnl_math_min( m_Metric->GetNumberOfParameters(), m_Scales.GetSize() );

  if ( m_ScalesInitialized )
    {
    for( unsigned int j = 0; j < spaceDimension; j++)
      {
      derivative[j] /= m_Scales[j];
      }
    }

}


template <typename TFixedImage, typename TMovingImage, typename TMetric>
void
ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::GetValueAndDerivative(
const ParametersType & parameters,
MeasureType & value,
DerivativeType & derivative ) const
{

  if ( !m_Metric )
    {
    itkExceptionMacro(<<"Metric is not present");
    }

  // Delegate computation to the underlying Metric
  m_Metric->GetValueAndDerivative( parameters, value, derivative );

  const unsigned int spaceDimension = 
    vnl_math_min( m_Metric->GetNumberOfParameters(), m_Scales.GetSize() );

  if ( m_ScalesInitialized )
    {
    for( unsigned int j = 0; j < spaceDimension; j++)
      {
      derivative[j] /= m_Scales[j];
      }
    }

}


template <typename TFixedImage, typename TMovingImage, typename TMetric>
typename ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::MeasureType
ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::GetValue( const ParametersType & parameters ) const
{
  if ( !m_Metric )
    {
    itkExceptionMacro(<<"Metric is not present");
    }

  // Delegate computation to the underlying Metric
  return m_Metric->GetValue( parameters );

}

template <typename TFixedImage, typename TMovingImage, typename TMetric>
void
ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::SetScales( const ScalesType & scales )
{
  if ( !m_ScalesInitialized || scales != m_Scales )
    {
    m_Scales = scales;
    m_ScalesInitialized = true;
    this->Modified();
    }
}

template <typename TFixedImage, typename TMovingImage, typename TMetric>
void
ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::Initialize() throw ( ExceptionObject )
{
  if ( !m_Metric )
    {
    itkExceptionMacro(<<"Metric is not present");
    }

  // Setup the inputs of underlying metric
  m_Metric->SetFixedImage( this->GetFixedImage() );
  m_Metric->SetMovingImage( this->GetMovingImage() );
  m_Metric->SetFixedImageRegion( this->GetFixedImageRegion() );

  TransformType * tptr = const_cast<TransformType *>( this->m_Transform.GetPointer() );
  m_Metric->SetTransform( tptr );
  
  InterpolatorType * iptr = const_cast<InterpolatorType *>( this->m_Interpolator.GetPointer() );
  m_Metric->SetInterpolator( iptr );
  
  MovingImageMaskType * mptr = const_cast<MovingImageMaskType *>( this->m_MovingImageMask.GetPointer() );
  m_Metric->SetMovingImageMask( mptr );
  
  FixedImageMaskType * fptr = const_cast<FixedImageMaskType *>( this->m_FixedImageMask.GetPointer() );
  m_Metric->SetFixedImageMask( fptr );
  
  
  // Delegate to underlying metric
  m_Metric->Initialize();

}


template <typename TFixedImage, typename TMovingImage, typename TMetric>
const unsigned long &
ScaledDerivativeImageToImageMetricAdaptor<TFixedImage,TMovingImage,TMetric>
::GetNumberOfPixelsCounted() const
{
  if ( !m_Metric )
    {
    itkExceptionMacro(<<"Metric is not present");
    }

  return m_Metric->GetNumberOfPixelsCounted();
}


} // end namespace itk

#endif
