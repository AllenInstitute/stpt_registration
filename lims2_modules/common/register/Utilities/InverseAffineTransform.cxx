/*=========================================================================

  InverseAffineTransform.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkAffineTransform.h"
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

  typedef itk::AffineTransform<double,2> Affine2DTransformType; 
  typedef itk::AffineTransform<double,3> Affine3DTransformType;
  
  try
    {
    if ( dimension == 2 )
      {
      Affine2DTransformType::Pointer it, ot;
      itk::idp::ReadTransform<Affine2DTransformType>( inputFile.c_str(), it );

      ot = Affine2DTransformType::New();
      ot->SetCenter( it->GetCenter() );
      it->GetInverse( ot );
      itk::idp::WriteTransform<Affine2DTransformType>( outputFile.c_str(), ot );
      }
    else if ( dimension == 3 )
      {
      Affine3DTransformType::Pointer it, ot;
      itk::idp::ReadTransform<Affine3DTransformType>( inputFile.c_str(), it );

      ot = Affine3DTransformType::New();
      ot->SetCenter( it->GetCenter() );
      it->GetInverse( ot );
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


