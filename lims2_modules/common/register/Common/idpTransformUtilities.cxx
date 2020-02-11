/*=========================================================================

  idpTransformUtilities.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpTransformUtilities_cxx
#define __idpTransformUtilities_cxx

#include "idpTransformUtilities.h"


namespace itk
{
namespace idp
{

void
TransformUtilities::
Compose( 
const TranslationTransformType * first,
const TranslationTransformType * second,
TranslationTransformPointer & output )
{
  output = TranslationTransformType::New();
  ParametersType op( output->GetNumberOfParameters() );
  op = first->GetParameters() + second->GetParameters();
  output->SetParameters( op );
}

void
TransformUtilities::
Compose( 
const TranslationTransformType * first,
const RigidTransformType * second,
RigidTransformPointer & output )
{
  typedef RigidTransformType::InputVectorType VectorType;
  VectorType v1;
  
  ParametersType fp = first->GetParameters();

  for ( int k = 0; k < 3; k++ )
    {
    v1[k] = fp[k];
    }

  VectorType v2 = second->TransformVector( v1 );

  ParametersType sp = second->GetParameters();

  for ( int k = 0; k < 3; k++ )
    {
    sp[6+k] += v2[k];
    }

  output = RigidTransformType::New();
  output->SetParameters( sp );

}

void
TransformUtilities::
Compose( 
const RigidTransformType * first,
const RigidTransformType * second,
RigidTransformPointer & output )
{
  ParametersType fp = first->GetParameters();
  ParametersType sp = second->GetParameters();
  
  output = RigidTransformType::New();
  ParametersType op( output->GetNumberOfParameters() );
  op.Fill( 0.0 );

  for ( int k = 0; k < 3; k++ )
    {
    op[k] = fp[k] + sp[k];
    }

  for( int k = 3; k < 6; k++ )
    {
    op[k] = sp[k];
    }

  typedef RigidTransformType::InputVectorType VectorType;
  VectorType v1;
  VectorType v2;

  for( int k = 0; k < 3; k++ )
    {
    v1[k] = second->GetCenter()[k] - first->GetCenter()[k];
    v2[k] = fp[k+6] - second->GetCenter()[k] + first->GetCenter()[k];
    }

  VectorType v3 = first->TransformVector( v1 );
  VectorType v4 = second->TransformVector( v2 );
  VectorType v5 = second->TransformVector( v3 );

  for( int k = 0; k < 3; k++ )
    {
    op[k+6] = sp[k+6] + v5[k] + v4[k];
    }  

  output->SetParameters( op );

 // std::cout << first->GetParameters() << std::endl;
 // std::cout << second->GetParameters() << std::endl;
 // std::cout << output->GetParameters() << std::endl;

}


void
TransformUtilities::
Compose( 
const RigidTransformType * first,
const Rigid2DTransformType * second,
RigidTransformPointer & output )
{
  RigidTransformType::Pointer dummy = RigidTransformType::New();
  ParametersType dp( dummy->GetNumberOfParameters() );
  dp.Fill( 0.0 );

  ParametersType sp = second->GetParameters();
  ParametersType fp = first->GetParameters();
  dp[2] = sp[0];

  for( int k = 0; k < 2; k++ )
    {
    dp[k+3] = sp[k+1];
    }
  dp[5] = fp[5];
  
  for ( int k = 0; k < 2; k++ )
    {
    dp[k+6] = sp[k+3];
    }
  dummy->SetParameters( dp );

  Compose( first, dummy, output );

}

void
TransformUtilities::
Compose( 
const Rigid2DTransformType * first,
const RigidTransformType * second,
RigidTransformPointer & output )
{
  RigidTransformType::Pointer dummy = RigidTransformType::New();
  ParametersType dp( dummy->GetNumberOfParameters() );
  dp.Fill( 0.0 );

  ParametersType sp = second->GetParameters();
  ParametersType fp = first->GetParameters();
  dp[2] = fp[0];

  for( int k = 0; k < 2; k++ )
    {
    dp[k+3] = fp[k+1];
    }
  dp[5] = sp[5];
  
  for ( int k = 0; k < 2; k++ )
    {
    dp[k+6] = fp[k+3];
    }
  dummy->SetParameters( dp );

  Compose( dummy, second, output );

}


} // end namespace idp
} //end namespace itk

#endif

  
