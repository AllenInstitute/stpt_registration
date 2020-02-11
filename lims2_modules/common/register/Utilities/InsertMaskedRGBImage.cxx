/*=========================================================================

  InsertMaskedRGBImage.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include <string>
#include "itkImageRegionIterator.h"
#include "idpRegistrationUtilities.h"
#include "itkAffineTransform.h"

int main( int argc, char *argv[] )
{
  if (argc < 8 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "backgroundImage insertImage outputImage ";
    std::cout << "maskImage foregroundRed foregroundGreen foregroundBlue";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::RGBPixel< unsigned char > PixelType;
  typedef itk::Image< PixelType, 2 >     ImageType;
  typedef itk::Image<unsigned char, 2>   MaskImageType;
  typedef itk::ImageRegionIterator<ImageType> Iterator;
  typedef itk::ImageRegionIterator<MaskImageType> MaskIterator;


  std::string backgroundFile = argv[1];
  std::string insertFile = argv[2];
  std::string outputFile = argv[3];
  std::string maskFile = argv[4];

  PixelType foreground;
  for ( int i = 0; i < 3; i++ )
    {
    foreground[i] = atoi( argv[5+i] );
    }  

  // --- read in backgroundImage and downsample mask
  ImageType::Pointer backgroundImage;
  ImageType::Pointer insertImage;
  ImageType::Pointer maskImage;
  double scale[2];

  try
    {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader1 = ReaderType::New();
    ReaderType::Pointer reader2 = ReaderType::New();
    ReaderType::Pointer reader3 = ReaderType::New();

    reader1->SetFileName( backgroundFile.c_str() );
    reader2->SetFileName( maskFile.c_str() );
    reader3->SetFileName( insertFile.c_str() );

    reader1->Update();
    reader2->Update();
    reader3->Update();

    backgroundImage = reader1->GetOutput();
    backgroundImage->DisconnectPipeline();

    maskImage = reader2->GetOutput();
    maskImage->DisconnectPipeline();

    insertImage = reader3->GetOutput();
    insertImage->DisconnectPipeline();

    for ( unsigned int j = 0; j < 2; j++ )
      {
      scale[j] = static_cast<double>( maskImage->GetBufferedRegion().GetSize(j) ) /
                 static_cast<double>( backgroundImage->GetBufferedRegion().GetSize(j) );
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

  // ---- Compute mask region interest and upsample mask to match background image
  typedef ImageType::RegionType RegionType;
  RegionType region;
  MaskImageType::Pointer binaryMask;

  try
    {
    ImageType::IndexType minIndex;
    ImageType::IndexType maxIndex;
    RegionType downsampleRegion;

    minIndex.Fill( 1000000000 );
    maxIndex.Fill( 0 );

    binaryMask = MaskImageType::New();
    binaryMask->CopyInformation( maskImage );
    binaryMask->SetRegions( maskImage->GetBufferedRegion() );
    binaryMask->Allocate();
    binaryMask->FillBuffer( 0 );

    Iterator iter( maskImage, maskImage->GetBufferedRegion() );

    MaskIterator miter( binaryMask, binaryMask->GetBufferedRegion() );

    while( !iter.IsAtEnd() )
      {
      if ( iter.Get() == foreground )
        {
        miter.Set( 255 );

        ImageType::IndexType index = iter.GetIndex();
        for ( unsigned int j = 0; j < 2; j++ )
          {
          if ( index[j] < minIndex[j] ) { minIndex[j] = index[j]; }
          if ( index[j] > maxIndex[j] ) { maxIndex[j] = index[j]; }
          }
        }
      ++iter;
      ++miter;
      }

    for ( unsigned int j = 0; j < 2; j++ )
      {
      downsampleRegion.SetIndex(j, minIndex[j]);
      downsampleRegion.SetSize(j , maxIndex[j] - minIndex[j] + 1);
      region.SetIndex(j, vnl_math_rnd( static_cast<double>( downsampleRegion.GetIndex(j) ) / scale[0] ));
      region.SetSize(j, vnl_math_rnd( static_cast<double>( downsampleRegion.GetSize(j) ) / scale[0] ));
      }

    std::cout << region << std::endl;

    typedef itk::AffineTransform<double,2> TransformType;
    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    transform->Scale( scale[0] );

    MaskImageType::Pointer dummy = MaskImageType::New();
    dummy->CopyInformation( backgroundImage );
    dummy->SetRegions( backgroundImage->GetBufferedRegion() );
    dummy->Allocate();

    itk::idp::ResampleImage<MaskImageType>( 
      binaryMask, dummy, transform, binaryMask, 0, "NearestNeighbor" );

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
  
  try
    {
    Iterator biter( backgroundImage, region );
    Iterator iiter( insertImage, region );
    MaskIterator miter( binaryMask, region );

    while( !miter.IsAtEnd() )
      {
      if ( miter.Get() )
        {
        biter.Set( iiter.Get() );
        }
      ++biter;
      ++iiter;
      ++miter;
      }

    itk::idp::WriteImage<ImageType>( outputFile.c_str(), backgroundImage );

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

