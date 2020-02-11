/*=========================================================================

  ComposeAffineTransform.cxx 

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkMatrixOffsetTransformBase.h"
#include "idpRegistrationUtilities.h"
#include <string>

int main( int argc, char *argv[] )
{
  if (argc < 4)
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "input1File input2File outputFile [pre] [dim]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string input1File = argv[1];
  std::string input2File = argv[2];
  std::string outputFile = argv[3];
  bool pre = false;
  unsigned int dimension = 3;

  if ( argc > 4 )
    {
    pre = atoi( argv[4] );
    }

  if ( argc > 5 )
    {
    dimension = atoi( argv[5] );
    }

  typedef itk::MatrixOffsetTransformBase<double,2,2> Base2DTransformType;
  typedef itk::AffineTransform<double,2> Affine2DTransformType;
  
  typedef itk::MatrixOffsetTransformBase<double,3,3> Base3DTransformType;
  typedef itk::AffineTransform<double,3> Affine3DTransformType;
  
  try
    {
    if ( dimension == 2 )
      {
      Base2DTransformType::Pointer it1, it2;
      itk::idp::ReadTransform<Base2DTransformType>( input1File.c_str(), it1 );
      itk::idp::ReadTransform<Base2DTransformType>( input2File.c_str(), it2 );
      it1->Compose( it2, pre );
      itk::idp::WriteTransform<Base2DTransformType>( outputFile.c_str(), it1 );
      }
    else if ( dimension == 3 )
      {
      Base3DTransformType::Pointer it1, it2;
      itk::idp::ReadTransform<Base3DTransformType>( input1File.c_str(), it1 );
      itk::idp::ReadTransform<Base3DTransformType>( input2File.c_str(), it2 );
      it1->Compose( it2, pre );
      itk::idp::WriteTransform<Base3DTransformType>( outputFile.c_str(), it1 );
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


