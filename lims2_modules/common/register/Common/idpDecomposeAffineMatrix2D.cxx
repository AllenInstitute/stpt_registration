/*=========================================================================

  idpDecomposeAffineMatrix2D.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpDecomposeAffineMatrix2D_cxx
#define __idpDecomposeAffineMatrix2D_cxx

#include "idpDecomposeAffineMatrix2D.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageMaskSpatialObject.h"
#include "itkTranslationTransform.h"

#include <vnl/vnl_math.h>

namespace itk
{
namespace idp
{

/**
 * Constructor
 */
DecomposeAffineMatrix2D::
DecomposeAffineMatrix2D ()
{
  m_InputMatrix.SetIdentity();
  m_OutputMatrix.SetIdentity();

  m_Scale.Fill( 1.0 );
  m_Shear = 0.0;
  m_Rotation = 0.0;

  m_UseScale = true;
  m_UseShear = true;
  m_UseRotation = true;
  
}

/**
 * Destructor
 */
DecomposeAffineMatrix2D::
~DecomposeAffineMatrix2D ()
{

}

/**
 * PrintSelf
 */
void 
DecomposeAffineMatrix2D::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << m_InputMatrix << std::endl;
  os << indent << m_OutputMatrix << std::endl;
  os << indent << m_Scale << std::endl;
  os << indent << m_Shear << std::endl;
  os << indent << m_Rotation << std::endl;
  os << indent << m_UseScale << std::endl;
  os << indent << m_UseShear << std::endl;
  os << indent << m_UseRotation << std::endl;

}


/**
 * Initiate registration
 */
void
DecomposeAffineMatrix2D::
Compute()
{
    MatrixType A;
    MatrixType S;
    MatrixType G;
    MatrixType R;

    A = m_InputMatrix.GetTranspose();
    A *= m_InputMatrix;

    m_Scale[0] = vcl_sqrt( A[0][0] );
    m_Scale[1] = vcl_sqrt( A[1][1] - ( vnl_math_sqr( A[0][1] ) / A[0][0] ));
    m_Shear = A[0][1] / ( m_Scale[0] * m_Scale[1] );

    S.SetIdentity();
    S[0][0] = m_Scale[0];
    S[1][1] = m_Scale[1];

    G.SetIdentity();
    G[0][1] = m_Shear;
    
    R = G * S;
    R = m_InputMatrix * R.GetInverse();
    
    m_Rotation = asin( R[1][0] );

    R[0][0] = cos( m_Rotation );
    R[0][1] = -1 * sin( m_Rotation );
    R[1][0] = sin( m_Rotation );
    R[1][1] = R[0][0];

    m_OutputMatrix.SetIdentity();
    
    if ( m_UseRotation )
      {
      m_OutputMatrix *= R;
      }
    if ( m_UseShear )
      {
      m_OutputMatrix *= G;
      }
    if ( m_UseScale )
      {
      m_OutputMatrix *= S;
      }
                                                
}



} // end namespace idp
} //end namespace itk

#endif
