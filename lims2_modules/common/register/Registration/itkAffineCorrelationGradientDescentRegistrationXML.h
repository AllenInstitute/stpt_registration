/*=========================================================================

  itkAffineCorrelationGradientDescentRegistrationXML.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkAffineCorrelationGradientDescentRegistrationXML_h
#define _itkAffineCorrelationGradientDescentRegistrationXML_h

#include "itkAffineCorrelationGradientDescentRegistration.h"
#include "tinyxml.h"

namespace itk
{

/** \class AffineCorrelationGradientDescentRegistrationXML
 * \brief Affine registration using correlation and gradient descent optimization.
 */
template < typename TFixedImage, typename TMovingImage >
class ITK_EXPORT AffineCorrelationGradientDescentRegistrationXML : 
  public AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef AffineCorrelationGradientDescentRegistrationXML Self;
  typedef AffineCorrelationGradientDescentRegistration<TFixedImage,TMovingImage> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(AffineCorrelationGradientDescentRegistrationXML, AffineCorrelationGradientDescentRegistration);

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Populate parameters from XML node. */
  virtual void PopulateFromXML( TiXmlNode * node ) throw (ExceptionObject);

protected:
  AffineCorrelationGradientDescentRegistrationXML ();
  virtual ~AffineCorrelationGradientDescentRegistrationXML () {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  AffineCorrelationGradientDescentRegistrationXML(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
    
};


} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAffineCorrelationGradientDescentRegistrationXML.txx"
#endif

#endif
