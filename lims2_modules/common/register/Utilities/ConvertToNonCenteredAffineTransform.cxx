/*=========================================================================

  ConvertToNonCenteredAffineTransform.cxx
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkMatrixOffsetTransformBase.h"
#include "idpRegistrationUtilities.h"
#include <string>

int main( int argc, char *argv[] )
{
  if (argc < 3)
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "inputFile outputFile [dim]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFile = argv[1];
  std::string outputFile = argv[2];

  unsigned int dimension = 3;

  if ( argc > 3 )
    {
    dimension = atoi( argv[3] );
    }

  typedef itk::MatrixOffsetTransformBase<double,2,2> Base2DTransformType;
  typedef itk::AffineTransform<double,2> Affine2DTransformType;
  
  typedef itk::MatrixOffsetTransformBase<double,3,3> Base3DTransformType;
  typedef itk::AffineTransform<double,3> Affine3DTransformType;
  
  try
    {
    if ( dimension == 2 )
      {
      Base2DTransformType::Pointer it;
      itk::idp::ReadTransform<Base2DTransformType>( inputFile.c_str(), it );

      Affine2DTransformType::Pointer ot = Affine2DTransformType::New();
      Affine2DTransformType::InputPointType center;
      center.Fill( 0.0 );
      ot->SetCenter( center );
      ot->SetMatrix( it->GetMatrix() );
      ot->SetOffset( it->GetOffset() );

      itk::idp::WriteTransform<Affine2DTransformType>( outputFile.c_str(), ot );
      }
    else if ( dimension == 3 )
      {
      Base3DTransformType::Pointer it;
      itk::idp::ReadTransform<Base3DTransformType>( inputFile.c_str(), it );

      Affine3DTransformType::Pointer ot = Affine3DTransformType::New();
      Affine3DTransformType::InputPointType center;
      center.Fill( 0.0 );
      ot->SetCenter( center );
      ot->SetMatrix( it->GetMatrix() );
      ot->SetOffset( it->GetOffset() );

      itk::idp::WriteTransform<Affine3DTransformType>( outputFile.c_str(), ot );
      }
    else
      {
      std::cerr << "dimension = " << dimension << " not supported" << std::endl;
      }
    
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


