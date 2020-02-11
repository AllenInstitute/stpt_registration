/*=========================================================================

  idpXMLUtilities.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpXMLUtilities_h
#define __idpXMLUtilities_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include <tinyxml.h>

namespace itk
{
namespace idp
{

/** \class XMLUtilities
 *
 */
class XMLUtilities : public Object
{
public:

  /** Standard typedefs. */
  typedef XMLUtilities      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(XMLUtilities, Object);
  
  /** Helper to attach new text node  . */
  template<class T>
  static TiXmlElement * InsertTextNode( TiXmlNode * node, 
                                        const char * value,
                                        const T& input )
   {
   std::ostringstream ostr;
   TiXmlElement * e = new TiXmlElement( value );
   ostr << input;
   TiXmlText * t = new TiXmlText( ostr.str().c_str() );
   e->LinkEndChild( t );
   node->LinkEndChild( e );
   return e;
   }


protected:
  XMLUtilities(){};
  ~XMLUtilities(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  XMLUtilities(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace idp
} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "idpImageSeriesUtilities.txx"
#endif

#endif
