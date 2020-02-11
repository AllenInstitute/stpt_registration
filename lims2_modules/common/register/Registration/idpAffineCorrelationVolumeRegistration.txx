/*=========================================================================

  idpAffineCorrelationVolumeRegistration.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpAffineCorrelationVolumeRegistration_txx
#define __idpAffineCorrelationVolumeRegistration_txx

#include "idpAffineCorrelationVolumeRegistration.h"
#include "tinyxml.h"

namespace itk
{
namespace idp
{

/**
 * Constructor
 */
template< typename TVolumeType >
AffineCorrelationVolumeRegistration<TVolumeType>::
AffineCorrelationVolumeRegistration ()
{
  m_FixedVolume = 0;
  m_FixedMask = 0;
  m_MovingVolume = 0;
  m_MovingMask = 0;

  m_Helper = HelperType::New();
  m_Worker = WorkerType::New();
  m_InputTransform = TransformType::New();
  m_OutputTransform = TransformType::New();

  m_Helper->SetVerbose( false );
  m_Worker->SetVerbose( false );
  
}

/**
 * Destructor
 */
template< typename TVolumeType >
AffineCorrelationVolumeRegistration<TVolumeType>::
~AffineCorrelationVolumeRegistration ()
{

}

/**
 * PrintSelf
 */
template< typename TVolumeType >
void 
AffineCorrelationVolumeRegistration<TVolumeType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << m_FixedVolume << std::endl;
  os << indent << m_FixedMask << std::endl;
  os << indent << m_MovingVolume << std::endl;
  os << indent << m_MovingMask << std::endl;
  os << indent << m_Helper << std::endl;
  os << indent << m_Worker << std::endl;
  os << indent << m_InputTransform << std::endl;
  os << indent << m_OutputTransform << std::endl;

}


/**
 * Load registration parameters from XML
 */
template< typename TVolumeType >
void
AffineCorrelationVolumeRegistration<TVolumeType>::
LoadParametersFromXML( const char * xmlFile )
{

  // Open xml file
  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load file: " << xmlFile );
    }

  TiXmlNode *moduleNode;
  TiXmlNode *helperNode;
  TiXmlNode *registrationNode;

  // Locate the root node
  moduleNode = doc.FirstChild("Module");
  if ( !moduleNode )
    {
    itkExceptionMacro( << "Cannot find Module node" );
    }

  // Locate the helper node
  helperNode = moduleNode->FirstChild("Helper");
  if ( !helperNode )
    {
    itkExceptionMacro( << "Cannot find Helper node" );
    }
  
  // Locate the registration node
  registrationNode = moduleNode->FirstChild("Registration");
  if ( !registrationNode )
    {
    itkExceptionMacro( << "Cannot find Registration node" );
    }

  try
    {
    m_Helper->PopulateFromXML( helperNode );
    m_Worker->PopulateFromXML( registrationNode );
    }
  catch( itk::ExceptionObject & err )
    {
    throw err;
    }
  catch( ... )
    {
    itkExceptionMacro( << "Caught unknown exception" );
    }


}

/**
 * Initiate registration
 */
template< typename TVolumeType >
void
AffineCorrelationVolumeRegistration<TVolumeType>::
Compute()
{

  try
    {
    m_Helper->SetFixedImage( m_FixedVolume );
    m_Helper->SetMovingImage( m_MovingVolume );
    m_Helper->SetFixedMaskImage( m_FixedMask );
    m_Helper->SetMovingMaskImage( m_MovingMask );

    m_Worker->SetCoordinator( m_Helper );
    m_Worker->SetTransform( m_InputTransform );
    m_Worker->SetInitialTransformParameters( m_InputTransform->GetParameters() );

    m_Worker->StartRegistration();
    m_OutputTransform = m_Worker->GetAffineTransform();
    m_OutputTransform->SetParameters( m_Worker->GetLastTransformParameters() );

    }
  catch( itk::ExceptionObject & err )
    {
    throw err;
    }
  catch( ... )
    {
    itkExceptionMacro( << "Caught unknown exception" );
    }

}



} // end namespace idp
} //end namespace itk

#endif
