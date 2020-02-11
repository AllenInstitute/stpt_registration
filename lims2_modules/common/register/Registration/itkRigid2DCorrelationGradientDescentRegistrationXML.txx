/*=========================================================================

  itkRigid2DCorrelationGradientDescentRegistrationXML.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkRigid2DCorrelationGradientDescentRegistrationXML_txx
#define _itkRigid2DCorrelationGradientDescentRegistrationXML_txx

#include "itkRigid2DCorrelationGradientDescentRegistrationXML.h"
#include <string>
#include <vector>

namespace itk
{

template< typename TFixedImage, typename TMovingImage >
Rigid2DCorrelationGradientDescentRegistrationXML<TFixedImage,TMovingImage>
::Rigid2DCorrelationGradientDescentRegistrationXML()
{
  
}

template< typename TFixedImage, typename TMovingImage >
void
Rigid2DCorrelationGradientDescentRegistrationXML<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );

}


template< typename TFixedImage, typename TMovingImage >
void
Rigid2DCorrelationGradientDescentRegistrationXML<TFixedImage,TMovingImage>
::PopulateFromXML(TiXmlNode * node) throw (ExceptionObject)
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXMLNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while( ( n = node->IterateChildren( n ) ) )
    {
    e = n->ToElement();


#define _XMLSetMacro(name,type) \
  ( param.compare( #name ) == 0 ) \
    { \
    this->Set##name( static_cast<type>( atof(e->GetText() ) ) ); \
    if( this->GetVerbose() ) { std::cout << this->Get##name(); } \
    }

    if ( e )
      {
      std::string param = e->Value();
      if ( this->GetVerbose() ) { std::cout << param << ": "; }
      if _XMLSetMacro( NumberOfLevels, unsigned long )
      else if _XMLSetMacro( AngleScale, double )
      else if _XMLSetMacro( AngleScaleRate, double )
      else if _XMLSetMacro( CenterScale, double )
      else if _XMLSetMacro( StepSize, double )
      else if _XMLSetMacro( StepSizeRate, double )
      else if _XMLSetMacro( NumberOfIterations, unsigned long )
      else if _XMLSetMacro( IterationRate, double )
      else if _XMLSetMacro( Verbose, bool )
      else
        {
        if ( this->GetVerbose() ) { std::cout << e->GetText() << " [Not used]"; }
        }
      if ( this->GetVerbose() ) { std::cout << std::endl; }
      }

#undef _XMLSetMacro

    }

}


} // namespace itk

#endif
