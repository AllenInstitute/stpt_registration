/*=========================================================================

  idpAffineCorrelationVolumeRegistration.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpAffineCorrelationVolumeRegistration_h
#define __idpAffineCorrelationVolumeRegistration_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkSimpleRegistrationImageCoordinatorXML.h"
#include "itkAffineCorrelationGradientDescentRegistrationXML.h"


namespace itk
{
namespace idp
{

/** \class AffineCorrelationVolumeRegistration
 *
 */
template< typename TVolumeType >
class AffineCorrelationVolumeRegistration: public Object
{
public:

  /** Standard typedefs. */
  typedef AffineCorrelationVolumeRegistration Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(AffineCorrelationVolumeRegistration, Object);
  
  /** Other types. */
  typedef TVolumeType VolumeType;
  typedef TVolumeType MaskType;

  typedef SimpleRegistrationImageCoordinatorXML< VolumeType,
                                                 VolumeType > HelperType;
  typedef AffineCorrelationGradientDescentRegistrationXML< VolumeType,
                                                           VolumeType > WorkerType;
  typedef typename WorkerType::AffineTransformType TransformType;

  /** Get the helper. **/
  itkGetObjectMacro( Helper, HelperType );

  /** Get the worker. **/
  itkGetObjectMacro( Worker, WorkerType );

  /** Set/Get fixed and moving volume/mask. */
  itkSetObjectMacro( FixedVolume, VolumeType );
  itkGetConstObjectMacro( FixedVolume, VolumeType );
  itkSetObjectMacro( FixedMask, MaskType );
  itkGetConstObjectMacro( FixedMask, MaskType );
  itkSetObjectMacro( MovingVolume, VolumeType );
  itkGetConstObjectMacro( MovingVolume, VolumeType );
  itkSetObjectMacro( MovingMask, MaskType );
  itkGetConstObjectMacro( MovingMask, MaskType );

  /** Set/Get the input transform. */
  itkSetObjectMacro( InputTransform, TransformType );
  itkGetConstObjectMacro( InputTransform, TransformType );

  /** Get the output transform. */
  itkGetConstObjectMacro( OutputTransform, TransformType );

  /** Load registration parameters from XML file. */
  void LoadParametersFromXML( const char * fileName );

  /** Initiate registration. */
  void Compute();

protected:
  AffineCorrelationVolumeRegistration();
  ~AffineCorrelationVolumeRegistration();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  AffineCorrelationVolumeRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename VolumeType::Pointer     m_FixedVolume;
  typename MaskType::Pointer       m_FixedMask;
  typename VolumeType::Pointer     m_MovingVolume;
  typename MaskType::Pointer       m_MovingMask;

  typename HelperType::Pointer     m_Helper;
  typename WorkerType::Pointer     m_Worker;
  typename TransformType::Pointer  m_InputTransform;
  typename TransformType::Pointer  m_OutputTransform;

};

} // end namespace idp
} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "idpAffineCorrelationVolumeRegistration.txx"
#endif

#endif
