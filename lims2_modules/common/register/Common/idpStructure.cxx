/*=========================================================================

  idpStructure.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpStructure_cxx
#define __idpStructure_cxx

#include "idpStructure.h"
#include "itkNumericTraits.h"
#include "tinyxml.h"
#include <algorithm>
#include <list>


namespace itk
{
namespace idp
{


/**
 * Initialize static variable
 */
unsigned long  Structure::m_NodeCount = 0;

/**
 * Constructor
 */
Structure::Structure()
{
  m_Id = 0;
  m_Acronym = "";
  m_Name = "";
  m_Red = 0;
  m_Green = 0;
  m_Blue = 0;
  m_Order = -1;
  m_Level = -1;
  m_Parent = 0;
  m_Children.resize( 0 );
  m_IdMap.clear();
  m_AcronymMap.clear();  
}

/**
 * Destructor
 */
Structure::~Structure()
{

}

/**
 * PrintSelf
 */
void 
Structure::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Id: " << this->GetId() << std::endl;
  os << indent << "Acronym: " << this->GetAcronym() << std::endl;
  os << indent << "Name: " << this->GetName() << std::endl;
  os << indent << "Red: " << static_cast<unsigned long>( this->GetRed() )<< std::endl;
  os << indent << "Green: " << static_cast<unsigned long>( this->GetGreen() )<< std::endl;
  os << indent << "Blue: " << static_cast<unsigned long>( this->GetBlue() )<< std::endl;
  os << indent << "Order: " << this->GetOrder() << std::endl;
  os << indent << "Level: " << this->GetLevel() << std::endl;
  os << indent << "NumberOfChildren: " << m_Children.size() << std::endl;
  os << indent << "Parent: " << m_Parent.GetPointer() << std::endl;
}


/** Operator to comparing the order of two structures
 *
 */
bool StructureLessThan(
  Structure::Pointer s1,
  Structure::Pointer s2 )
{
  return s1->GetOrder() < s2->GetOrder();
}


/** 
 * Load data from XML file
 */
void
Structure::LoadFromXML( const char * xmlFile )
{
  // Open xml file
  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile );
    }

  TiXmlNode * graphNode;

  // Locate the root node
  graphNode = doc.FirstChild( "structure_graph" );
  if ( !graphNode )
    {
    itkExceptionMacro( << "Can not find structure_graph node in file " << xmlFile );
    }
    
  // Locate the root structure
  TiXmlNode * n;
  n = graphNode->FirstChild( "structure" );
  if ( !n )
    {
    itkExceptionMacro( << "Can not find structure node in graph " );     
    }
   
  m_NodeCount = 0;
  this->SetOrder( m_NodeCount++ ); 
  this->LoadFromXML( n );

}

/** 
 * Load data from XML node
 */
void
Structure::LoadFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();
      if ( parameter.compare( "id" ) == 0 )
        {
        this->SetId( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "acronym" ) == 0 )
        {
        this->SetAcronym( e->GetText() );
        //std::cout << this->GetAcronym() << std::endl;
        }
       else if ( parameter.compare( "name" ) == 0 )
        {
        this->SetName( e->GetText() );
        }
       else if ( parameter.compare( "red" ) == 0 )
        {
        this->SetRed( atoi( e->GetText() ) );
        }
       else if ( parameter.compare( "green" ) == 0 )
        {
        this->SetGreen( atoi( e->GetText() ) );
        }
       else if ( parameter.compare( "blue" ) == 0 )
        {
        this->SetBlue( atoi( e->GetText() ) );
        }
      // else if ( parameter.compare( "st_order" ) == 0 )
      //  {
      //  this->SetOrder( atoi( e->GetText() ) );
      //  }
       else if ( parameter.compare( "st_level" ) == 0 )
        {
        if (e->GetText())
          {
          this->SetLevel( atoi( e->GetText() ) );
          }
        }
       else if ( parameter.compare( "children" ) == 0 )
        {
        this->LoadChildrenFromXML( n );
        }

      } // end if (e)   

    } // end while

}

/** 
 * Load data from XML node
 */
void
Structure::LoadChildrenFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();
      if ( parameter.compare( "structure" ) == 0 )
        {
        this->AddChild();
        m_Children.back()->LoadFromXML( n );
        }

      } // end if (e)   

    } // end while

}

/** 
 * Add new child to the children vector
 */
void
Structure::AddChild()
{
  Pointer c = Self::New();
  c->SetOrder( m_NodeCount++ );
  c->SetParent( this );
  m_Children.push_back( c );
}


/** 
 * Populate Id map
 */
void
Structure::PopulateIdMap()
{
  m_IdMap.clear();

  typedef std::list<Pointer> ListType;
  ListType list;
  list.push_back( this );

  while( !list.empty() )
    {
    Pointer s = list.front();
    list.pop_front();
    m_IdMap[ s->GetId() ] = s;

    for( unsigned int j = 0; j < s->GetChildren().size(); j++ )
      {
      list.push_back( s->GetChildren()[j] );
      }    
    }

}

/** 
 * Populate Acronym map
 */
void
Structure::PopulateAcronymMap()
{
  m_AcronymMap.clear();

  typedef std::list<Pointer> ListType;
  ListType list;
  list.push_back( this );

  while( !list.empty() )
    {
    Pointer s = list.front();
    list.pop_front();
    m_AcronymMap[ s->GetAcronym() ] = s;

    for( unsigned int j = 0; j < s->GetChildren().size(); j++ )
      {
      list.push_back( s->GetChildren()[j] );
      }    
    }

}



} // end namespace idp
} //end namespace itk

#endif
