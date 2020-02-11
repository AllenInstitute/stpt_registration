/*=========================================================================

  itkRigid2DCorrelationGradientDescentRegistrationXML.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkRigid2DCorrelationGradientDescentRegistrationXML_h
#define _itkRigid2DCorrelationGradientDescentRegistrationXML_h

#include "itkRigid2DCorrelationGradientDescentRegistration.h"
#include "tinyxml.h"

namespace itk
{

/** \class Rigid2DCorrelationGradientDescentRegistrationXML
 * \brief Rigid2D registration using correlation and gradient descent optimization.
 */
template < typename TFixedImage, typename TMovingImage >
class ITK_EXPORT Rigid2DCorrelationGradientDescentRegistrationXML : 
  public Rigid2DCorrelationGradientDescentRegistration<TFixedImage,TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef Rigid2DCorrelationGradientDescentRegistrationXML Self;
  typedef Rigid2DCorrelationGradientDescentRegistration<TFixedImage,TMovingImage> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(Rigid2DCorrelationGradientDescentRegistrationXML, Rigid2DCorrelationGradientDescentRegistration);

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Populate parameters from XML node. */
  virtual void PopulateFromXML( TiXmlNode * node ) throw (ExceptionObject);

protected:
  Rigid2DCorrelationGradientDescentRegistrationXML ();
  virtual ~Rigid2DCorrelationGradientDescentRegistrationXML () {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  Rigid2DCorrelationGradientDescentRegistrationXML(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
    
};


} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRigid2DCorrelationGradientDescentRegistrationXML.txx"
#endif

#endif
