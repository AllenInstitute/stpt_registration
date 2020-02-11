/*=========================================================================

  idpDecomposeAffineMatrix3D.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpDecomposeAffineMatrix3D_cxx
#define __idpDecomposeAffineMatrix3D_cxx

#include "idpDecomposeAffineMatrix3D.h"
#include <vnl/vnl_math.h>

namespace itk
{
namespace idp
{

/**
 * Constructor
 */
DecomposeAffineMatrix3D::
DecomposeAffineMatrix3D ()
{
  m_InputMatrix.SetIdentity();
  m_OutputMatrix.SetIdentity();

  m_Scale.Fill( 1.0 );
  m_Shear.Fill( 0.0 );;
  m_Rotation.Fill( 0.0 );

  m_UseScale = true;
  m_UseShear = true;
  m_UseRotation = true;
  
}

/**
 * Destructor
 */
DecomposeAffineMatrix3D::
~DecomposeAffineMatrix3D ()
{

}

/**
 * PrintSelf
 */
void 
DecomposeAffineMatrix3D::
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
DecomposeAffineMatrix3D::
Compute()
{
    MatrixType A;
    MatrixType S;
    MatrixType G, G0, G1, G2;
    MatrixType R;

    A = m_InputMatrix.GetTranspose();
    A *= m_InputMatrix;

    m_Scale[1] = A[1][1];
    m_Scale[1] -= vnl_math_sqr(A[1][2]) / A[2][2];
    m_Scale[1] = vcl_sqrt( m_Scale[1] );

    m_Scale[0] = A[0][0];
    m_Scale[0] -= vnl_math_sqr( A[0][1] * A[2][2] - A[0][1] * A[1][2] ) / vnl_math_sqr( A[2][2] * m_Scale[1] );
    m_Scale[0] = vcl_sqrt( m_Scale[0] );

    m_Scale[2] = A[2][2];
    m_Scale[2] -= vnl_math_sqr( A[0][2] ) / vnl_math_sqr( m_Scale[0] );
    m_Scale[2] = vcl_sqrt( m_Scale[0] );
    
    double D = m_Scale[0] * m_Scale[1] * m_Scale[2];
  
    m_Shear[0] = A[0][2] * m_Scale[1] / D;
    m_Shear[1] = ( A[0][1]*A[2][2] - A[0][2]*A[1][2] ) * m_Scale[2] / A[2][2] / D;
    m_Shear[2] = A[1][2] * vnl_math_sqr( m_Scale[2] ) * m_Scale[0] / A[2][2] / D;

    S.SetIdentity();
    S[0][0] = m_Scale[0];
    S[1][1] = m_Scale[1];
    S[2][2] = m_Scale[2];

    G0.SetIdentity();
    G0[0][2] = m_Shear[0];

    G1.SetIdentity();
    G1[1][0] = m_Shear[1];

    G2.SetIdentity();
    G2[2][1] = m_Shear[2];

    G = G0 * G1 * G2;
    
    R = G * S;
    R = m_InputMatrix * R.GetInverse();
    
    m_Rotation[1] = - asin( R[0][2] );
    m_Rotation[0] = asin( R[1][2] / cos(m_Rotation[1]) );
    m_Rotation[2] = asin( R[0][1] / cos(m_Rotation[1]) );

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
