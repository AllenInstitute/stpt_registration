/*=========================================================================

  idpReferenceSpaceUtilities.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpReferenceSpaceUtilities_cxx
#define __idpReferenceSpaceUtilities_cxx

#include "idpReferenceSpaceUtilities.h"
#include <utility>
#include "idpXMLUtilities.h"

namespace itk
{
namespace idp
{

/**
 * Constructor
 */
ReferenceSpaceUtilities::
ReferenceSpaceUtilities ()
{
  m_PointMap.clear();
  m_RegionOfInterestMap.clear();  
}

/**
 * Destructor
 */
ReferenceSpaceUtilities::
~ReferenceSpaceUtilities ()
{

}

/**
 * PrintSelf
 */
void 
ReferenceSpaceUtilities::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

}


/**
 * Insert point
 */
void 
ReferenceSpaceUtilities::
InsertPoint(
const char * label,
const PointType & point )
{
  m_PointMap.insert( std::make_pair( label, point ) );
}

/**
 * Insert region of interest
 */
void 
ReferenceSpaceUtilities::
InsertRegionOfInterest(
const char * label,
const RegionOfInterestType & roi )
{
  m_RegionOfInterestMap.insert( std::make_pair( label, roi ) );
}

/**
 * Find point from label
 */
ReferenceSpaceUtilities::
PointType
ReferenceSpaceUtilities::
GetPoint(
const char * label )
{
  return m_PointMap[label];
}


/**
 * Find region of interest from label
 */
ReferenceSpaceUtilities::
RegionOfInterestType
ReferenceSpaceUtilities::
GetRegionOfInterest(
const char * label )
{
  return m_RegionOfInterestMap[label];
}

/** 
 * Check if a label exists in the point map
 */
bool
ReferenceSpaceUtilities::
PointLabelExists( 
const char * label )
{
  return ( m_PointMap.count(label) > 0 );
}


/** 
 * Check if a label exists in the point of interest map
 */
bool
ReferenceSpaceUtilities::
RegionOfInterestLabelExists( 
const char * label )
{
  return ( m_RegionOfInterestMap.count(label) > 0 );
}

/** 
 * Insert points and region of interests from xml file
 */
void
ReferenceSpaceUtilities::
InsertFromXML( 
const char * xmlFile )
{
  // open xml file
  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile );
    }

  TiXmlNode * node;

  // Locate the root node
  node = doc.FirstChild( "fiducials" );
  if ( !node )
    {
    itkExceptionMacro( << "Can not find model node in file " << xmlFile );
    }

  TiXmlNode *n;
  TiXmlElement *e;
  n = NULL;

  while( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {

      std::string parameter = e->Value();
      if ( parameter.compare( "point" ) == 0 )
        {
        this->InsertPoint( n );
        }
      else if ( parameter.compare( "region_of_interest" ) == 0 )
        {
        this->InsertRegionOfInterest( n );
        }

      } // end if (e)

    } // end while
}

/** 
 * Insert point from xml node
 */
void
ReferenceSpaceUtilities::
InsertPoint( 
TiXmlNode * node) throw (ExceptionObject)
{

  if ( !node )
    {
    itkExceptionMacro( << "TiXMLNode is NULL" );
    }

  std::string label;
  PointType   point;
  point.Fill( 0.0 );

  TiXmlNode *n;
  TiXmlElement *e;
   
  n = NULL;

  while( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {

      std::string parameter = e->Value();
      if ( parameter.compare( "id" ) == 0 )
        {
        label = e->GetText();
        }
      else if ( parameter.compare( "x" ) == 0 )
        {
        point[0] = atof( e->GetText() );
        }
      else if ( parameter.compare( "y" ) == 0 )
        {
        point[1] = atof( e->GetText() );
        }
      else if ( parameter.compare( "z" ) == 0 )
        {
        point[2] = atof( e->GetText() );
        }

      } // end if (e)

    } // end while

   this->InsertPoint( label.c_str(), point );

}


/** 
 * Insert region of interests from xml node
 */
void
ReferenceSpaceUtilities::
InsertRegionOfInterest( 
TiXmlNode * node) throw (ExceptionObject)
{

  if ( !node )
    {
    itkExceptionMacro( << "TiXMLNode is NULL" );
    }

  std::string label;
  RegionOfInterestType  roi;
  roi.StartPoint.Fill( 0.0 );
  roi.EndPoint.Fill( 0.0 );

  TiXmlNode *n;
  TiXmlElement *e;
   
  n = NULL;

  while( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {

      std::string parameter = e->Value();
      if ( parameter.compare( "id" ) == 0 )
        {
        label = e->GetText();
        }
      else if ( parameter.compare( "start_x" ) == 0 )
        {
        roi.StartPoint[0] = atof( e->GetText() );
        }
      else if ( parameter.compare( "start_y" ) == 0 )
        {
        roi.StartPoint[1] = atof( e->GetText() );
        }
      else if ( parameter.compare( "start_z" ) == 0 )
        {
        roi.StartPoint[2] = atof( e->GetText() );
        }
      else if ( parameter.compare( "end_x" ) == 0 )
        {
        roi.EndPoint[0] = atof( e->GetText() );
        }
      else if ( parameter.compare( "end_y" ) == 0 )
        {
        roi.EndPoint[1] = atof( e->GetText() );
        }
      else if ( parameter.compare( "end_z" ) == 0 )
        {
        roi.EndPoint[2] = atof( e->GetText() );
        }

      } // end if (e)

    } // end while

   this->InsertRegionOfInterest( label.c_str(), roi );


}


/** 
 * Write to xml
 */
void
ReferenceSpaceUtilities::
WriteToXML( 
const char * xmlFile )
{
  // open the doc
  TiXmlDocument doc;
  TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "UTF-8", "" );
  doc.LinkEndChild( decl );

  // set up root (fiducials) node
  TiXmlElement * fiducials = new TiXmlElement( "fiducials" );
  doc.LinkEndChild( fiducials );

  PointMap::const_iterator pit = m_PointMap.begin();
  PointMap::const_iterator pend = m_PointMap.end();

  for ( ; pit != pend; pit++ )
    {
    TiXmlElement * e = new TiXmlElement( "point" );
    XMLUtilities::InsertTextNode<const char *>( e, "id", (*pit).first.c_str() );
    XMLUtilities::InsertTextNode<double>( e, "x", (*pit).second[0] );
    XMLUtilities::InsertTextNode<double>( e, "y", (*pit).second[1] );
    XMLUtilities::InsertTextNode<double>( e, "z", (*pit).second[2] );
    fiducials->LinkEndChild( e );
    }

  RegionOfInterestMap::const_iterator rit = m_RegionOfInterestMap.begin();
  RegionOfInterestMap::const_iterator rend = m_RegionOfInterestMap.end();

  for ( ; rit != rend; rit++ )
    {
    TiXmlElement * e = new TiXmlElement( "region_of_interest" );
    XMLUtilities::InsertTextNode<const char *>( e, "id", (*rit).first.c_str() );
    XMLUtilities::InsertTextNode<double>( e, "start_x", (*rit).second.StartPoint[0] );
    XMLUtilities::InsertTextNode<double>( e, "start_y", (*rit).second.StartPoint[1] );
    XMLUtilities::InsertTextNode<double>( e, "start_z", (*rit).second.StartPoint[2] );
    XMLUtilities::InsertTextNode<double>( e, "end_x", (*rit).second.EndPoint[0] );
    XMLUtilities::InsertTextNode<double>( e, "end_y", (*rit).second.EndPoint[1] );
    XMLUtilities::InsertTextNode<double>( e, "end_z", (*rit).second.EndPoint[2] );
    fiducials->LinkEndChild( e );
    }
  
  // close file
  doc.SaveFile( xmlFile );

}



} // end namespace idp
} //end namespace itk

#endif

  
