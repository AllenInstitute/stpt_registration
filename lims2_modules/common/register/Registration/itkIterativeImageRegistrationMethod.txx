/*=========================================================================

  itkIterativeImageRegistrationMethod.txx
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkIterativeImageRegistrationMethod_txx
#define _itkIterativeImageRegistrationMethod_txx

#include "itkIterativeImageRegistrationMethod.h"

namespace itk
{

/*
 * Constructor
 */
template < typename TFixedImage, typename TMovingImage >
IterativeImageRegistrationMethod<TFixedImage,TMovingImage>
::IterativeImageRegistrationMethod()
{

  m_FixedImage   = 0; // has to be provided by the user.
  m_MovingImage  = 0; // has to be provided by the user.
  m_Transform    = 0; // has to be provided by the user.
  m_Interpolator = 0; // has to be provided by the user.
  m_Metric       = 0; // has to be provided by the user.
  m_Optimizer    = 0; // has to be provided by the user.

  m_NumberOfLevels = 1;
  m_CurrentLevel = 0;

  m_Stop = false;

  m_InitialTransformParameters = ParametersType(1);
  m_InitialTransformParametersOfNextLevel = ParametersType(1);
  m_LastTransformParameters = ParametersType(1);

  m_InitialTransformParameters.Fill( 0.0f );
  m_InitialTransformParametersOfNextLevel.Fill( 0.0f );
  m_LastTransformParameters.Fill( 0.0f );

}


/*
 * Initialize by setting the interconnects between components. 
 */
template < typename TFixedImage, typename TMovingImage >
void
IterativeImageRegistrationMethod<TFixedImage,TMovingImage>
::Initialize() throw (ExceptionObject)
{

  // Sanity checks
  if ( !m_Metric )
    {
    itkExceptionMacro(<<"Metric is not present" );
    }

  if ( !m_Optimizer )
    {
    itkExceptionMacro(<<"Optimizer is not present" );
    }

  if( !m_Transform )
    {
    itkExceptionMacro(<<"Transform is not present");
    }

  if( !m_Interpolator )
    {
    itkExceptionMacro(<<"Interpolator is not present");
    }

  // Setup the metric
  m_Metric->SetMovingImage( m_MovingImage );
  m_Metric->SetFixedImage( m_FixedImage );
  m_Metric->SetTransform( m_Transform );
  m_Metric->SetInterpolator( m_Interpolator );
  m_Metric->SetFixedImageRegion( m_FixedImageRegion );
  m_Metric->Initialize();

  // Setup the optimizer
  m_Optimizer->SetCostFunction( m_Metric );
  m_Optimizer->SetInitialPosition( m_InitialTransformParametersOfNextLevel );


}


/*
 * Stop the Registration Process
 */
template < typename TFixedImage, typename TMovingImage >
void
IterativeImageRegistrationMethod<TFixedImage,TMovingImage>
::StopRegistration( void )
{
  m_Stop = true;
}


/*
 * Starts the Registration Process
 */
template < typename TFixedImage, typename TMovingImage >
void
IterativeImageRegistrationMethod<TFixedImage,TMovingImage>
::StartRegistration( void )
{ 

  this->BeforeRegistration();

  m_InitialTransformParametersOfNextLevel = m_InitialTransformParameters;

  if ( m_InitialTransformParametersOfNextLevel.Size() != 
       m_Transform->GetNumberOfParameters() )
    {
    itkExceptionMacro(<<"Size mismatch between initial parameter and transform"); 
    }

  m_Stop = false;

  for ( m_CurrentLevel = 0; m_CurrentLevel < m_NumberOfLevels;
        m_CurrentLevel++ )
    {

    this->BeforeIteration();

    // Invoke an iteration event.
    // This allows a UI to reset any of the components between
    // resolution level.
    this->InvokeEvent( IterationEvent() );

    // Check if there has been a stop request
    if ( m_Stop ) 
      {
      break;
      }

    try
      {
      // initialize the interconnects between components
      this->Initialize();
      }
    catch( ExceptionObject& err )
      {
      m_LastTransformParameters = ParametersType(1);
      m_LastTransformParameters.Fill( 0.0f );

      // pass exception to caller
      throw err;
      }

    try
      {
      // do the optimization
      m_Optimizer->StartOptimization();
      }
    catch( ExceptionObject& err )
      {
      // An error has occurred in the optimization.
      // Update the parameters
      m_LastTransformParameters = m_Optimizer->GetCurrentPosition();

      // Pass exception to caller
      throw err;
      }

    // get the results
    m_LastTransformParameters = m_Optimizer->GetCurrentPosition();
    m_Transform->SetParameters( m_LastTransformParameters );

    this->AfterIteration();

    // setup the initial parameters for next level
    if ( m_CurrentLevel < m_NumberOfLevels - 1 )
      {
      m_InitialTransformParametersOfNextLevel =
        m_LastTransformParameters;
      }
    }

    this->AfterRegistration();

}


/*
 * PrintSelf
 */
template < typename TFixedImage, typename TMovingImage >
void
IterativeImageRegistrationMethod<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Metric: " << m_Metric.GetPointer() << std::endl;
  os << indent << "Optimizer: " << m_Optimizer.GetPointer() << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "FixedImage: " << m_FixedImage.GetPointer() << std::endl;
  os << indent << "MovingImage: " << m_MovingImage.GetPointer() << std::endl;

  os << indent << "NumberOfLevels: ";
  os << m_NumberOfLevels << std::endl;

  os << indent << "CurrentLevel: ";
  os << m_CurrentLevel << std::endl;  

  os << indent << "InitialTransformParameters: ";
  os << m_InitialTransformParameters << std::endl;
  os << indent << "InitialTransformParametersOfNextLevel: ";
  os << m_InitialTransformParametersOfNextLevel << std::endl;
  os << indent << "LastTransformParameters: ";
  os << m_LastTransformParameters << std::endl;
  os << indent << "FixedImageRegion: ";
  os << m_FixedImageRegion << std::endl;

}




} // end namespace itk


#endif
