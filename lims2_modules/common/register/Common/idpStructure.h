/*=========================================================================

  idpStructure.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpStructure_h
#define __idpStructure_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include <string>
#include <algorithm>
#include <vector>
#include <map>

#include "tinyxml.h"

namespace itk
{
namespace idp
{

/** \class Structure
 *  \brief Encapsulate structure object
 *
 */
class Structure : public Object
{
public:
  /** Standard typedefs. */
  typedef Structure      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(Structure, Object);

  /** Other types. */
  typedef std::vector<Pointer>  StructureVectorType;
  typedef std::map< std::string, Pointer > AcronymMapType;
  typedef std::map< unsigned long, Pointer > IdMapType;

  /** Set/Get the structure id. */
  itkSetMacro( Id, unsigned long );
  itkGetConstMacro( Id, unsigned long );
  
  /** Set/Get the acronym. */
  itkSetStringMacro( Acronym );
  itkGetStringMacro( Acronym );
  
  /** Set/Get the name. */
  itkSetStringMacro( Name );
  itkGetStringMacro( Name );

  /** Set/Get the red value. */
  itkSetMacro( Red, unsigned char );
  itkGetConstMacro( Red, unsigned char );

  /** Set/Get the green value. */
  itkSetMacro( Green, unsigned char );
  itkGetConstMacro( Green, unsigned char );

  /** Set/Get the blue value. */
  itkSetMacro( Blue, unsigned char );
  itkGetConstMacro( Blue, unsigned char );
  
  /** Set/Get the structure order. */
  itkSetMacro( Order, signed long );
  itkGetConstMacro( Order, signed long );

  /** Set/Get the structure level. */
  itkSetMacro( Level, signed long );
  itkGetConstMacro( Level, signed long );

  /** Set/Get the parent. */
  itkSetObjectMacro( Parent, Self );
  itkGetConstObjectMacro( Parent, Self );
  itkGetObjectMacro( Parent, Self );

  /** Get the children. */
  StructureVectorType & GetChildren()
    { return m_Children; }
  const StructureVectorType & GetChildren() const
    { return m_Children; }

  /** Add a new child to the back of the child vector. */
  virtual void AddChild();

      
  /** Load information from XML. */
  virtual void LoadFromXML( const char * xmlFile );
  virtual void LoadFromXML( TiXmlNode * node );
  virtual void LoadChildrenFromXML( TiXmlNode * node );

  /** Populate maps. */
  virtual void PopulateAcronymMap();
  virtual void PopulateIdMap();
  
  /** Get maps. */ 
  AcronymMapType & GetAcronymMap()
    { return m_AcronymMap; }
  const AcronymMapType & GetAcronymMap() const
    { return m_AcronymMap; }
  IdMapType & GetIdMap()
    { return m_IdMap; }
  const IdMapType & GetIdMap() const
    { return m_IdMap; }
  

protected:
  Structure();
  ~Structure();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  Structure(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  unsigned long         m_Id;
  std::string           m_Acronym; 
  std::string           m_Name;
  unsigned char         m_Red;
  unsigned char         m_Green;
  unsigned char         m_Blue;
  signed long           m_Order;
  signed long           m_Level;
  StructureVectorType   m_Children;
  Pointer               m_Parent;
  AcronymMapType        m_AcronymMap;
  IdMapType             m_IdMap;

  static unsigned long  m_NodeCount;


};

bool StructureLessThan (
  Structure::Pointer s1,
  Structure::Pointer s2 );


} // end namespace idp
} //end namespace itk

#endif
