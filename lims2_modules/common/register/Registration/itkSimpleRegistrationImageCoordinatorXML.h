/*=========================================================================

  itkSimpleRegistrationImageCoordinatorXML.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkSimpleRegistrationImageCoordinatorXML_h
#define _itkSimpleRegistrationImageCoordinatorXML_h

#include "itkSimpleRegistrationImageCoordinator.h"
#include "tinyxml.h"

namespace itk
{

/** \class SimpleRegistrationImageCoordinatorXML
 * \brief Helper class to coordinate images for IterativeRegistrationMethod
 */
template < typename TFixedImage, typename TMovingImage >
class ITK_EXPORT SimpleRegistrationImageCoordinatorXML : 
  public SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef SimpleRegistrationImageCoordinatorXML Self;
  typedef SimpleRegistrationImageCoordinator<TFixedImage,TMovingImage> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(SimpleRegistrationImageCoordinatorXML, SimpleRegistrationImageCoordinator);

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Type of defining the shrink factors. */
  typedef typename Superclass::FactorsType FactorsType;

  /** Populate parameters from XML node. */
  virtual void PopulateFromXML( TiXmlNode * node ) throw (ExceptionObject);

protected:
  SimpleRegistrationImageCoordinatorXML ();
  virtual ~SimpleRegistrationImageCoordinatorXML () {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  SimpleRegistrationImageCoordinatorXML(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
    
};


} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSimpleRegistrationImageCoordinatorXML.txx"
#endif

#endif
