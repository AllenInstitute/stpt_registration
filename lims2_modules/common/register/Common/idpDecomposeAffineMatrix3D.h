/*=========================================================================

  idpDecomposeAffineMatrix3D.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpDecomposeAffineMatrix3D_h
#define __idpDecomposeAffineMatrix3D_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkMatrix.h"
#include "itkVector.h"

namespace itk
{
namespace idp
{

/** \class DecomposeAffineMatrix3D
 *
 */
class DecomposeAffineMatrix3D : public Object
{
public:

  /** Standard typedefs. */
  typedef DecomposeAffineMatrix3D      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(DecomposeAffineMatrix3D, Object);
  
  /** Other types. */
  typedef Matrix<double,3,3> MatrixType;
  typedef Vector<double,3>   VectorType;

  /** Set/Get input matrix. */
  itkSetMacro( InputMatrix, MatrixType );
  itkGetConstMacro( InputMatrix, MatrixType );

  /** Set/Get output matrix. */
  itkSetMacro( OutputMatrix, MatrixType );
  itkGetConstMacro( OutputMatrix, MatrixType );

  /** Set/Get the scale. */
  itkSetMacro( Scale, VectorType );
  itkGetConstMacro( Scale, VectorType );

  /** Set/Get the shear. */
  itkSetMacro( Shear, VectorType );
  itkGetConstMacro( Shear, VectorType );

  /** Set/Get the rotation. */
  itkSetMacro( Rotation, VectorType );
  itkGetConstMacro( Rotation, VectorType );

  /** Set/Get use scale flag. */
  itkSetMacro( UseScale, bool );
  itkGetConstMacro( UseScale, bool );

  /** Set/Get use shear flag. */
  itkSetMacro( UseShear, bool );
  itkGetConstMacro( UseShear, bool );

  /** Set/Get use rotation flag. */
  itkSetMacro( UseRotation, bool );
  itkGetConstMacro( UseRotation, bool );

  /** Perform decomposition and recomposition */
  void Compute();

protected:
  DecomposeAffineMatrix3D();
  ~DecomposeAffineMatrix3D();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  DecomposeAffineMatrix3D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MatrixType              m_InputMatrix;
  MatrixType              m_OutputMatrix;
  VectorType              m_Scale;
  VectorType              m_Shear;
  VectorType              m_Rotation;
  bool                    m_UseScale;
  bool                    m_UseShear;
  bool                    m_UseRotation;
 
};

} // end namespace idp
} //end namespace itk

#endif
