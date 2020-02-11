/*=========================================================================

  idpTransformUtilities.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpTransformUtilities_h
#define __idpTransformUtilities_h

#include "itkTranslationTransform.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkCenteredRigid2DTransform.h"


namespace itk
{
namespace idp
{

/** \class TransformUtilities
 *
 */
class TransformUtilities
{
public:

  /* transform types */
  typedef TranslationTransform<double,3>      TranslationTransformType;
  typedef TranslationTransformType::Pointer   TranslationTransformPointer;

  typedef CenteredEuler3DTransform<double>    RigidTransformType;
  typedef RigidTransformType::Pointer         RigidTransformPointer;

  typedef CenteredRigid2DTransform<double>    Rigid2DTransformType;
  typedef Rigid2DTransformType::Pointer       Rigid2DTransformPointer;

  typedef TransformBase::ParametersType       ParametersType;

  static void Compose( const TranslationTransformType * first,
                       const TranslationTransformType * second,
                       TranslationTransformPointer & output );

  static void Compose( const TranslationTransformType * first,
                       const RigidTransformType * second,
                       RigidTransformPointer & output );

  static void Compose( const RigidTransformType * first,
                       const TranslationTransformType * second,
                       RigidTransformPointer & output );

  static void Compose( const RigidTransformType * first,
                       const RigidTransformType * second,
                       RigidTransformPointer & output );

  static void Compose( const RigidTransformType * first,
                       const Rigid2DTransformType * second,
                       RigidTransformPointer & output );

  static void Compose( const Rigid2DTransformType * first,
                       const RigidTransformType * second,
                       RigidTransformPointer & output );


};

} // end namespace idp
} //end namespace itk
           
#endif
  
