/*=========================================================================
  
  idpDecomposeAffineMatrix2D.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpDecomposeAffineMatrix2D_h
#define __idpDecomposeAffineMatrix2D_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkMatrix.h"
#include "itkVector.h"

namespace itk
{
namespace idp
{

/** \class DecomposeAffineMatrix2D
 *
 */
class DecomposeAffineMatrix2D : public Object
{
public:

  /** Standard typedefs. */
  typedef DecomposeAffineMatrix2D      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(DecomposeAffineMatrix2D, Object);
  
  /** Other types. */
  typedef Matrix<double,2,2> MatrixType;
  typedef Vector<double,2>   VectorType;

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
  itkSetMacro( Shear, double );
  itkGetConstMacro( Shear, double );

  /** Set/Get the rotation. */
  itkSetMacro( Rotation, double );
  itkGetConstMacro( Rotation, double );

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
  DecomposeAffineMatrix2D();
  ~DecomposeAffineMatrix2D();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  DecomposeAffineMatrix2D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MatrixType              m_InputMatrix;
  MatrixType              m_OutputMatrix;
  Vector<double,2>        m_Scale;
  double                  m_Shear;
  double                  m_Rotation;
  bool                    m_UseScale;
  bool                    m_UseShear;
  bool                    m_UseRotation;
 
};

} // end namespace idp
} //end namespace itk

#endif
