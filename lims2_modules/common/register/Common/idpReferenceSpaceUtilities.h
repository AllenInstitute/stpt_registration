/*=========================================================================

  idpReferenceSpaceUtilities.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpReferenceSpaceUtilities_h
#define __idpReferenceSpaceUtilities_h


#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkPoint.h"
#include "itkVector.h"
#include <string>
#include <map>
#include "tinyxml.h"

namespace itk
{
namespace idp
{

template< class TPoint>
class RegionOfInterest
{
  public:
    TPoint  StartPoint;
    TPoint  EndPoint;
};

/** \class ReferenceSpaceUtilities
 *
 */
class ReferenceSpaceUtilities : public Object
{
public:
  /** Standard typedefs. */
  typedef ReferenceSpaceUtilities      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(ReferenceSpaceUtilities, Object);
  
  /** Other types. */
  static const unsigned int Dimension = 3;
  typedef double CoordRepType;
  typedef Point<CoordRepType,Dimension> PointType;
  typedef RegionOfInterest<PointType> RegionOfInterestType;

  typedef std::map< std::string, PointType> PointMap;
  typedef std::map< std::string, RegionOfInterestType> RegionOfInterestMap;

  /** Inserting points  and region of interests */
  void InsertPoint( const char * label, const PointType & point );
  void InsertRegionOfInterest( const char * label, const RegionOfInterestType & roi );
  void InsertFromXML( const char * filename );
  void InsertPoint( TiXmlNode * node ) throw (ExceptionObject);
  void InsertRegionOfInterest( TiXmlNode * node ) throw (ExceptionObject);

  /** Output to xml */
  void WriteToXML( const char * filename );

  /** Fetching points and region of interests */
  PointType GetPoint( const char * label );
  RegionOfInterestType GetRegionOfInterest( const char * label );

  /** Label */
  bool PointLabelExists( const char * label );
  bool RegionOfInterestLabelExists( const char * label );

  /** Get PointMap */
  const PointMap & GetPointMap() const
   { return m_PointMap; }

  /** Get RegionOfInterestMap */
  const RegionOfInterestMap & GetRegionOfInterestMap() const
   { return m_RegionOfInterestMap; }

protected:
  ReferenceSpaceUtilities();
  ~ReferenceSpaceUtilities();
  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  ReferenceSpaceUtilities(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  PointMap                m_PointMap;
  RegionOfInterestMap     m_RegionOfInterestMap;

};

} // end namespace idp
} //end namespace itk
           
#endif

  
