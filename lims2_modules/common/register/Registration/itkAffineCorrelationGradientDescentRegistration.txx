/*=========================================================================

  itkAffineCorrelationGradientDescentRegistration.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkAffineCorrelationGradientDescentRegistration_txx
#define _itkAffineCorrelationGradientDescentRegistration_txx

#include "itkAffineCorrelationGradientDescentRegistration.h"

namespace itk
{

/*
 * Constructor
 */
template < typename TFixedImage, typename TMovingImage >
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::AffineCorrelationGradientDescentRegistration()
{

  // setup default values
  m_Coordinator = NULL;

  m_MatrixScale              = 250.0 * 250.0 * 100;
  m_MatrixScaleRate          = 0.1;
  m_CenterScale             = 5000 * 5000;
  m_StepSize                = 25;
  m_StepSizeRate            = 0.25;
  m_NumberOfIterations      = 400;
  m_IterationRate           = 0.5;

  m_Verbose                 = true;  

  // setup the default components
  typename AffineTransformType::Pointer          transform    = AffineTransformType::New();
  typename CorrelationMetricType::Pointer metric       = CorrelationMetricType::New();
  typename MetricAdaptorType::Pointer           adaptor      = MetricAdaptorType::New();
  typename LinearInterpolatorType::Pointer      interpolator = LinearInterpolatorType::New();
  GradientDescentOptimizer::Pointer             optimizer    = GradientDescentOptimizer::New();

  this->SetTransform( transform );
  adaptor->SetMetric( metric );
  this->SetMetric( adaptor );
  this->SetInterpolator( interpolator );
  this->SetOptimizer( optimizer );

  // set up default values for metric
  metric->SubtractMeanOn();

  // set up default values for adaptor
  typename MetricAdaptorType::ScalesType scales( transform->GetNumberOfParameters() );
  scales.Fill( 1.0 );
  for( unsigned int j = 0; j < FixedImageDimension * FixedImageDimension; j++ )
    {
    scales[j] = m_MatrixScale;
    }
  for( unsigned int j = FixedImageDimension * FixedImageDimension; j < (FixedImageDimension * FixedImageDimension + FixedImageDimension); j++ )
    {
    scales[j] = m_CenterScale;
    }
  adaptor->SetScales( scales );

  // set up default values for optimizer
  GradientDescentOptimizer::ScalesType optimizerScales( transform->GetNumberOfParameters() );
  optimizerScales.Fill( 1.0 );
  optimizer->SetScales( optimizerScales );

  optimizer->SetMaximize( false );

  // Setup an optimizer observer
  typedef SimpleMemberCommand<Self> CommandType;
  typename CommandType::Pointer command = CommandType::New();
  command->SetCallbackFunction( this, &Self::OptimizerIterationHandler );

  m_ObserverTag = this->GetGradientDescentOptimizer()->AddObserver( IterationEvent(), command );

}

template < typename TFixedImage, typename TMovingImage >
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::~AffineCorrelationGradientDescentRegistration()
{
  this->GetGradientDescentOptimizer()->RemoveObserver( m_ObserverTag );
}


template < typename TFixedImage, typename TMovingImage >
void
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Coordinator: " 
               << m_Coordinator.GetPointer() << std::endl;
  os << indent << "MatrixScale: "
               << m_MatrixScale << std::endl;
  os << indent << "MatrixScaleRate: "
               << m_MatrixScaleRate << std::endl;
}

/*
 * Do setup before registration
 */
template < typename TFixedImage, typename TMovingImage >
void
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::BeforeRegistration( )
{
  if ( m_Coordinator.IsNull() )
    {
    itkExceptionMacro( << "Coordinator not set" );
    }

  m_Coordinator->Initialize();

  // Set the fixed amd moving image masks
  this->GetMetricAdaptor()->SetFixedImageMask( m_Coordinator->GetFixedImageMask() );
  this->GetMetricAdaptor()->SetMovingImageMask( m_Coordinator->GetMovingImageMask() );

   // Set parameter scales
   typedef typename MetricAdaptorType::ScalesType ScalesType;
   ScalesType scales = this->GetMetricAdaptor()->GetScales();
   scales.Fill( 1.0 );
  for( unsigned int j = 0; j < FixedImageDimension * FixedImageDimension; j++ )
    {
    scales[j] = m_MatrixScale;
    }
  for( unsigned int j = FixedImageDimension * FixedImageDimension; j < (FixedImageDimension * FixedImageDimension + FixedImageDimension); j++ )
    {
    scales[j] = m_CenterScale;
    }
   this->GetMetricAdaptor()->SetScales( scales );

   // Setup optimizer
  typedef GradientDescentOptimizer  OptimizerType;
  typename OptimizerType::Pointer optimizer = this->GetGradientDescentOptimizer();
  optimizer->SetLearningRate( m_StepSize );
  optimizer->SetNumberOfIterations( m_NumberOfIterations );


}


/*
 * Do setup before registration
 */
template < typename TFixedImage, typename TMovingImage >
void
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::BeforeIteration( )
{
   if ( m_Verbose )
    {
    std::cout << "Level: " << this->GetCurrentLevel() << std::endl;
    }
    
   // Set fixed and moving images for this level
   this->SetFixedImage( m_Coordinator->GetFixedImage( this->GetCurrentLevel() ) );

   typedef typename TFixedImage::RegionType RegionType;
   RegionType region = this->GetFixedImage()->GetBufferedRegion();

   this->SetFixedImageRegion( region );

   this->SetMovingImage( m_Coordinator->GetMovingImage( this->GetCurrentLevel() ) );

   // Setup optimizer
  if ( this->GetCurrentLevel() )
    {
    typedef GradientDescentOptimizer  OptimizerType;
    typename OptimizerType::Pointer optimizer = this->GetGradientDescentOptimizer();
    double rate = optimizer->GetLearningRate();
    rate *= m_StepSizeRate;
    optimizer->SetLearningRate( rate );
    double n = optimizer->GetNumberOfIterations();
    n *= m_IterationRate;
    optimizer->SetNumberOfIterations( static_cast<unsigned int>( n ) );
    }

   // Set parameter scales
  if ( this->GetCurrentLevel() )
    {
    typedef typename MetricAdaptorType::ScalesType ScalesType;
    ScalesType scales = this->GetMetricAdaptor()->GetScales();
    for( unsigned int j = 0; j < FixedImageDimension * FixedImageDimension; j++ )
      {
      scales[j] *= m_MatrixScaleRate;
      }
    this->GetMetricAdaptor()->SetScales( scales );
    }
 
}


/*
 * Do setup after registration
 */
template < typename TFixedImage, typename TMovingImage >
void
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::AfterIteration( )
{
 
}

/*
 * Do setup before registration
 */
template < typename TFixedImage, typename TMovingImage >
void
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::AfterRegistration( )
{


}

/*
 * Print out info each optimization iteration
 */
template < typename TFixedImage, typename TMovingImage >
void
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::OptimizerIterationHandler()
{
   if ( m_Verbose )
    {
    GradientDescentOptimizer::Pointer optimizer = this->GetGradientDescentOptimizer();

    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
    }

}

/*
 * Return a pointer to the internal metric adaptor
 */
template < typename TFixedImage, typename TMovingImage >
typename AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::MetricAdaptorType *
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::GetMetricAdaptor()
{
  MetricAdaptorType * ptr =
    dynamic_cast<MetricAdaptorType *>( this->GetMetric() );

  if ( !ptr ) 
    { 
    itkExceptionMacro( << "Metric is not of type ScaledDerivativeImageToImageMetricAdaptor" ); 
    }

  return ptr;
}


/*
 * Return a pointer to the internal metric
 */
template < typename TFixedImage, typename TMovingImage >
typename AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::CorrelationMetricType *
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::GetCorrelationMetric()
{
  CorrelationMetricType * ptr =
    dynamic_cast<CorrelationMetricType *>( this->GetMetricAdaptor()->GetMetric() );

  if ( !ptr ) 
    { 
    itkExceptionMacro( << "Metric is not of type NormalizedCorrelationImageToImageMetric" ); 
    }

  return ptr;
}


/*
 * Return a point to the internal transform
 */
template < typename TFixedImage, typename TMovingImage >
typename AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::AffineTransformType *
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::GetAffineTransform()
{
  AffineTransformType * ptr =
    dynamic_cast<AffineTransformType *>( this->GetTransform() );

  if ( !ptr ) 
    { 
    itkExceptionMacro( << "Metric is not of type Euler3DTransform" ); 
    }

  return ptr;
}

/*
 * Return a point to the internal interpolator
 */
template < typename TFixedImage, typename TMovingImage >
typename AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::LinearInterpolatorType *
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::GetLinearInterpolator()
{
  LinearInterpolatorType * ptr =
    dynamic_cast<LinearInterpolatorType *>( this->GetInterpolator() );

  if ( !ptr ) 
    { 
    itkExceptionMacro( << "Metric is not of type LinearInterpolateImageFunction" ); 
    }

  return ptr;
}

/*
 * Return a point to the internal optimizer
 */
template < typename TFixedImage, typename TMovingImage >
GradientDescentOptimizer *
AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
::GetGradientDescentOptimizer()
{
  GradientDescentOptimizer * ptr =
    dynamic_cast<GradientDescentOptimizer *>( this->GetOptimizer() );

  if ( !ptr ) 
    { 
    itkExceptionMacro( << "Metric is not of type GradientDescentOptimizer" ); 
    }

  return ptr;
}

} // end namespace itk


#endif
