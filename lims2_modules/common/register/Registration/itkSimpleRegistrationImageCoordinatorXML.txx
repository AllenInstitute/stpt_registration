/*=========================================================================

  itkSimpleRegistrationImageCoordinatorXML.txx
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkSimpleRegistrationImageCoordinatorXML_txx
#define _itkSimpleRegistrationImageCoordinatorXML_txx

#include "itkSimpleRegistrationImageCoordinatorXML.h"
#include <string>
#include <vector>

namespace itk
{

template< typename TFixedImage, typename TMovingImage >
SimpleRegistrationImageCoordinatorXML<TFixedImage,TMovingImage>
::SimpleRegistrationImageCoordinatorXML()
{
  
}

template< typename TFixedImage, typename TMovingImage >
void
SimpleRegistrationImageCoordinatorXML<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );

}

template<typename A>
void
PopulateFixedArray( 
TiXmlNode * node,
A & arr )
{
  if ( !node ) { return; }

  TiXmlNode *n;
  TiXmlElement *e;
  int        i;

  n = NULL;
  while( ( n = node->IterateChildren( "Element", n ) ) )
    {
    e = n->ToElement();
    const char * att = e->Attribute("index", &i);
    if ( !att || i < 0 ) { continue; }
    arr[i] = static_cast<typename A::ValueType>( atof( e->GetText() ) );
    }
}

template< typename TFixedImage, typename TMovingImage >
void
SimpleRegistrationImageCoordinatorXML<TFixedImage,TMovingImage>
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

    if ( e )
      {
      std::string param = e->Value();
      if ( this->GetVerbose() ) { std::cout << param << ": "; };
      if ( param.compare( "NumberOfLevels" ) == 0 )
        {
        this->SetNumberOfLevels( atoi( e->GetText() ) );
        if ( this->GetVerbose() ) { std::cout << this->GetNumberOfLevels(); };
        }
      else if ( param.compare( "FixedImageStartingFactors" ) == 0 )
        {
        FactorsType arr;
        PopulateFixedArray( n, arr );
        this->SetFixedImageStartingFactors( arr );
        if ( this->GetVerbose() ) { std::cout << this->GetFixedImageStartingFactors(); };
        }
      else if ( param.compare( "MovingImageStartingFactors" ) == 0 )
        {
        FactorsType arr;
        PopulateFixedArray( n, arr );
        this->SetMovingImageStartingFactors( arr );
        if ( this->GetVerbose() ) { std::cout << this->GetMovingImageStartingFactors(); };
        }
      else
        {
        if ( this->GetVerbose() ) { std::cout << e->GetText() << " [Not used]"; };
        }
      if ( this->GetVerbose() ) { std::cout << std::endl; };
      }

    }

}


} // namespace itk

#endif
