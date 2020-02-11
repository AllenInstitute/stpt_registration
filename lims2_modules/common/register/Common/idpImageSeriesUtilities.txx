/*=========================================================================

  idpImageSeriesUtilities.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpImageSeriesUtilities_txx
#define __idpImageSeriesUtilities_txx

#include "idpImageSeriesUtilities.h"
#include "idpRegistrationUtilities.h"
#include "idpTransformUtilities.h"

#include "itkRGBPixel.h"
#include "itkExtractImageFilter.h"
#include "itkRGBToGrayscalePixelAccessor.h"
#include "itkImageAdaptor.h"
#include "itkImageMomentsCalculator.h"
#include "itkTranslationTransform.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImportImageFilter.h"
#include "itkSliceImageConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkJP2ImageIO.h"

#include "tinyxml.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <itksys/SystemTools.hxx>

#include "idpXMLUtilities.h"
#include "idpDecomposeAffineMatrix3D.h"

namespace itk
{
namespace idp
{

/**
 * Constructor
 */
template<class TPixel>
ImageSeriesUtilities<TPixel>
::ImageSeriesUtilities ()
{
  m_ImageSeries = 0;
  m_Spacing.Fill( 1.0 );
  m_Origin.Fill( 0.0 );
  m_Size.Fill( 0 );
  m_Center.Fill( 0.0 );
  m_OutputWidth = 1.0;
  m_OutputHeight = 1.0;
  m_ModelDirectory = "";
  m_Verbose = false;
  m_DownsampleFactor = 16;
  m_GenerateMask = true;
  m_UseStandardSize = true;
  m_InvertDarkFieldImages = false;
}

/**
 * Destructor
 */
template<class TPixel>
ImageSeriesUtilities<TPixel>
::~ImageSeriesUtilities ()
{

}

/**
 * PrintSelf
 */
template<class TPixel>
void 
ImageSeriesUtilities<TPixel>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Size: " << this->GetSize() << std::endl;
  os << indent << "Spacing: " << this->GetSpacing() << std::endl;
  os << indent << "Origin: " << this->GetOrigin() << std::endl;
  os << indent << "Center: " << this->GetCenter() << std::endl;
  os << indent << "OutputWidth: " << this->GetOutputWidth() << std::endl;
  os << indent << "OutputHeight: " << this->GetOutputHeight() << std::endl;
  os << indent << "ModelDirectory: " << this->GetModelDirectory() << std::endl;

}


/**
 * Standard size per (age,plane) combination
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>
::GetStandardSizeFromXML()
{
  std::string xmlFile = this->GetModelDirectory();
  xmlFile += "/model.xml";

  // open xml file
  TiXmlDocument doc( xmlFile.c_str() );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile );
    }

  TiXmlNode * node;

  // Locate the root node
  node = doc.FirstChild( "model" );
  if ( !node )
    {
    itkExceptionMacro( << "Can not find model node in file " << xmlFile );
    }

  TiXmlNode *n;
  TiXmlElement *e;
  n = NULL;

  while( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {

      std::string parameter = e->Value();
      if ( parameter.compare( "output-width" ) == 0 )
        {
        this->SetOutputWidth( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "output-height" ) == 0 )
        {
        this->SetOutputHeight( atof( e->GetText() ) );
        }

      } // end if (e)

    } // end while
}


/**
 * Determine the volume spacing, origin and size
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>
::ComputeVolumeMetaInformation()
{
  
  const SubImageArray & subimages = m_ImageSeries->GetSubImages(); 

  if ( m_UseStandardSize )
    {
    this->GetStandardSizeFromXML();
    }
  else
    {
    // set size of image-series meta data.
    Vector<double,2> extend;
    extend.Fill( 0.0 );
    for ( unsigned int k = 0; k < subimages.size(); k++ )
      {
      if ( subimages[k]->GetFailed() )
        {
        continue;
        }
      for ( unsigned int i = 0; i < 2; i++ )
          {
          double f = subimages[k]->GetSize()[i] * subimages[k]->GetImageResolution();
          if ( f > extend[i] )
            {
            extend[i] = f;
            }
          }
      }
    this->SetOutputWidth( extend[0] );
    this->SetOutputHeight( extend[1] );
    }

  // Assume that subimages are already sorted by SpecimenTissueIndex
  signed int interval = m_ImageSeries->GetInterSectionInterval();

  // Get the first and last index with data
  signed int firstIndex = subimages[0]->GetSpecimenTissueIndex();
  signed int lastIndex = subimages[subimages.size()-1]->GetSpecimenTissueIndex();

  // Compute spacing, origin and region
  m_Spacing.Fill( m_ImageSeries->GetImageResolution() * static_cast<double>( m_DownsampleFactor ) );
  m_Spacing[2] = m_ImageSeries->GetSectionThickness() * static_cast<double>( interval );

  m_Size[0] = static_cast<long>( this->GetOutputWidth() / m_Spacing[0] );
  m_Size[1] = static_cast<long>( this->GetOutputHeight() / m_Spacing[1] );
  m_Size[2] = (( lastIndex - firstIndex ) / interval + 1 );

  m_Origin.Fill( 0.0 );
  m_Origin[2] = firstIndex * m_ImageSeries->GetSectionThickness();

  for ( unsigned int j = 0; j < Dimension; j++ )
    {
    m_Center[j] = m_Origin[j] + 0.5*(m_Size[j] -1) * m_Spacing[j];
    }

}

/**
 * Initialize the volume to the right size and defaults
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
InitializeVolume(
VolumePointer & volume,
MaskVolumePointer & mask,
MaskVolumePointer & sliceMask,
PixelType imgBkgnd )
{
  volume = VolumeType::New();
  volume->SetSpacing( this->GetSpacing() );
  volume->SetOrigin( this->GetOrigin() );
  volume->SetRegions( this->GetSize() );
  volume->Allocate();
  volume->FillBuffer( imgBkgnd );

  mask = MaskVolumeType::New();
  mask->SetSpacing( this->GetSpacing() );
  mask->SetOrigin( this->GetOrigin() );
  mask->SetRegions( this->GetSize() );
  mask->Allocate();
  mask->FillBuffer( 0 );

  sliceMask = MaskVolumeType::New();
  sliceMask->SetSpacing( this->GetSpacing() );
  sliceMask->SetOrigin( this->GetOrigin() );
  sliceMask->SetRegions( this->GetSize() );
  sliceMask->Allocate();
  sliceMask->FillBuffer( 0 );

}

/**
 * Populate volume and mask by resampling each section
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
PopulateVolume(
VolumePointer & volume,
MaskVolumePointer & mask,
MaskVolumePointer & sliceMask )
{

  double redScale = 0.0;
  double greenScale = 0.0;
  double blueScale = 1.0;

  if( this->GetImageSeries()->HasTreatment( "NISSL" ) )
    {
    redScale = 0.0;
    greenScale = 1.0;
    blueScale = 0.0;
    }

 this->PopulateVolume( volume, mask, sliceMask, redScale, greenScale, blueScale );

}

/**
 * Populate volume and mask by resampling each section
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
PopulateVolume(
VolumePointer & volume,
MaskVolumePointer & mask,
MaskVolumePointer & sliceMask,
double redScale,
double greenScale,
double blueScale,
PixelType imgBkgnd )
{

  this->InitializeVolume( volume, mask, sliceMask, imgBkgnd );  
  
  double reduce = log( static_cast<double>(m_DownsampleFactor) ) / log( 2.0 );

  std::cout << "reduce: " << reduce << std::endl;

  // resample each section and insert into volume
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }
  
    typename VolumeType::Pointer sImage;
    typename MaskVolumeType::Pointer sMask;
    typename VolumeType::Pointer sRef;
    typename MaskVolumeType::Pointer sMaskRef;
    typename VolumeType::Pointer rImage;
    typename MaskVolumeType::Pointer rMask;

    // Setup reference image
    try
      {
      typename VolumeType::SizeType sSize;
      PointType sOrigin;

      sSize = this->GetSize();
      sSize[2] = 1;
      sOrigin = this->GetOrigin();
      sOrigin[2] = static_cast<double>( subimages[k]->GetSpecimenTissueIndex() * 
                                        m_ImageSeries->GetSectionThickness() );

      sRef = VolumeType::New();
      sRef->SetSpacing( this->GetSpacing() );
      sRef->SetOrigin( sOrigin );
      sRef->SetRegions( sSize );

      sRef->Allocate();
      sRef->FillBuffer( 0 );
      MakeImage<MaskVolumeType,VolumeType>( sRef, sMaskRef, 0 );
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      continue;
      }
    catch( ... )
      {
      std::cerr << "Caught unknown exception" << std::endl;
      continue;
      }    
    
    // Compute downsampled 16 bounding box
    typename VolumeType::RegionType sRegion;
    for ( unsigned int j = 0; j < Dimension - 1; j++ )
      {
      sRegion.SetIndex( j, subimages[k]->GetIndex()[j] / m_DownsampleFactor );
      sRegion.SetSize( j , subimages[k]->GetSize()[j] / m_DownsampleFactor );
      }
    sRegion.SetIndex( 2, 0 );
    sRegion.SetSize( 2, 1 );

    // extract the downsampled image and mask
    try
      {
      typedef itk::Image< RGBPixelType, Dimension > RGBVolumeType;
      typename RGBVolumeType::Pointer tImage;

      std::string fileName = subimages[k]->GetDownsampleFilename();
      if ( fileName.empty() )
        {
        fileName = subimages[k]->GetTifFilename();
        } 

      if ( !fileName.empty() )
        {
        if ( !itksys::SystemTools::FileExists( fileName.c_str(), true ) )
          {
          itkExceptionMacro( << "Missing file " << fileName << std::endl );
          }

        typedef itk::ImageFileReader< RGBVolumeType > ReaderType;
        typename ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( fileName.c_str() );
        reader->Update();
        sRegion.Crop( reader->GetOutput()->GetBufferedRegion() );
        tImage = reader->GetOutput();
        }
      else
        {
        typename RGBImageType::Pointer t2Image;
        this->ReadImage<RGBImageType>( subimages[k], t2Image, static_cast<int>(reduce) );
        CastImage<RGBImageType,RGBVolumeType>( t2Image, tImage );
        sRegion.SetIndex(0,0); sRegion.SetIndex(1,0); sRegion.SetIndex(2,0); // offset + crop has already be taken care off by ReadImage
        sRegion.Crop( tImage->GetBufferedRegion() );
        }

      bool isDarkField = false;
      if ( this->GetInvertDarkFieldImages() )
        {
        // extract blue channel
        typename VolumeType::Pointer blue;
        VectorIndexSelection<RGBVolumeType,VolumeType>( tImage, 2, blue );
        // compute statistics 
        double min, max, sum, mean, variance, sigma;
        ImageStatistics<VolumeType>( blue, min, max, sum, mean, variance, sigma );
        //std::cout << subimages[k]->GetSpecimenTissueIndex() << ": " << mean << std::endl;
        if ( mean < 100.0 )
          {
          isDarkField = true;
          }
        }
      

      typedef itk::ImageAdaptor< RGBVolumeType, itk::RGBToGrayscalePixelAccessor<PixelType> > ImageAdaptorType;
      itk::RGBToGrayscalePixelAccessor<PixelType> accessor;
      accessor.SetScales( redScale, greenScale, blueScale );

      if ( this->GetInvertDarkFieldImages() && isDarkField )
        {
        accessor.SetScales( 0.0, 0.0, 1.0 );
        }

      typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
      adaptor->SetImage( tImage );
      adaptor->SetPixelAccessor( accessor );

      typedef itk::ExtractImageFilter<ImageAdaptorType,VolumeType> FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetInput( adaptor );
      filter->SetExtractionRegion( sRegion );

      filter->Update();
      sImage = filter->GetOutput();
      sImage->DisconnectPipeline();
      typename VolumeType::SpacingType spacing = this->GetSpacing();
      spacing[0] = this->GetDownsampleFactor() * subimages[k]->GetImageResolution();
      spacing[1] = spacing[0];
      sImage->SetSpacing( spacing );
      sImage->SetOrigin( sRef->GetOrigin() );
      sImage->SetRegions( sImage->GetBufferedRegion().GetSize() );

      if ( this->GetInvertDarkFieldImages() && isDarkField )
        {
        InvertIntensity<VolumeType>( sImage, sImage );
        }

      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      throw err;
      }
    catch( ... )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      itkExceptionMacro( << "Caught unknown exception" << std::endl );
      }

    for ( unsigned int j = 0; j < Dimension - 1; j++ )
      {
      sRegion.SetIndex( j, subimages[k]->GetIndex()[j] / m_DownsampleFactor );
      sRegion.SetSize( j , subimages[k]->GetSize()[j] / m_DownsampleFactor );
      }
    sRegion.SetIndex( 2, 0 );
    sRegion.SetSize( 2, 1 );

    try
      {

      if ( m_GenerateMask )
        {
        std::string fileName = subimages[k]->GetMaskFilename();
        ChangePaths( fileName );

        if ( fileName.empty() )
          {
          itkExceptionMacro( << "No mask file found for sub-image " << subimages[k]->GetId() << std::endl );
          }

        if ( !itksys::SystemTools::FileExists( fileName.c_str(), true ) )
          {
          itkExceptionMacro( << "Missing file " << fileName << std::endl );
          }

        typedef itk::ImageFileReader< VolumeType > ReaderType;
        typename ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( fileName.c_str() );
        reader->Update();
        sRegion.Crop( reader->GetOutput()->GetBufferedRegion() );

        typedef itk::ExtractImageFilter<VolumeType,MaskVolumeType> FilterType;
        typename FilterType::Pointer filter = FilterType::New();
        filter->SetInput( reader->GetOutput() );
        filter->SetExtractionRegion( sRegion );

        filter->Update();
        sMask = filter->GetOutput();
        sMask->DisconnectPipeline();
        typename VolumeType::SpacingType spacing = this->GetSpacing();
        spacing[0] = this->GetDownsampleFactor() * subimages[k]->GetImageResolution();
        spacing[1] = spacing[0];
        sMask->SetSpacing( spacing );
        sMask->SetOrigin( sRef->GetOrigin() );
        sMask->SetRegions( sMask->GetBufferedRegion().GetSize() );
        }
 

      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      throw err;
      }
    catch( ... )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      itkExceptionMacro( << "Caught unknown exception" << std::endl );
      }

    try
      {
      // resample section
      ResampleImage<VolumeType>( sImage, sRef, subimages[k]->GetTransform(),
                     rImage, imgBkgnd, "Linear" );

      if ( m_GenerateMask )
        {
        // resample mask
        ResampleImage<MaskVolumeType>( sMask, sMaskRef, subimages[k]->GetTransform(),
                     rMask, 0, "NearestNeighbor" );
        }
      else
        {
        MakeImage<MaskVolumeType,MaskVolumeType>( sMaskRef, rMask, 255 );
        }

      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      throw err;
      }
    catch( ... )
      {
      itkExceptionMacro( << "Caught unknown exception" << std::endl );
      }   

     // insert into volume 
    try
      {
      typename VolumeType::PointType dummyPt;
      dummyPt.Fill( 0.0 );
      dummyPt[2] = subimages[k]->GetSpecimenTissueIndex() * m_ImageSeries->GetSectionThickness();

      typename VolumeType::IndexType index;
      volume->TransformPhysicalPointToIndex( dummyPt, index );

      typename VolumeType::RegionType cRegion;
      cRegion.SetSize( sRef->GetBufferedRegion().GetSize() );
      cRegion.SetIndex( index );

      CopyRegion<VolumeType>( rImage, rImage->GetBufferedRegion(),
                              volume, cRegion );

      CopyRegion<MaskVolumeType>( rMask, rMask->GetBufferedRegion(),
                                  mask, cRegion );


      typedef itk::ImageRegionIterator<MaskVolumeType> Iterator;
      Iterator iter( sliceMask, cRegion );
      while( !iter.IsAtEnd() )
        {
        iter.Set( 255 );
        ++iter;
        }

      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      throw err;
      }
    catch( ... )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      itkExceptionMacro( << "Caught unknown exception" << std::endl );
      }   

    }

}

/**
 * Compute transforms to align centroid of each section
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
CentroidAlignment()
{

  // compute centroid at each section
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }

    // Compute downsampled 16 bounding box
    typename VolumeType::RegionType sRegion;
    for ( unsigned int j = 0; j < Dimension - 1; j++ )
      {
      sRegion.SetIndex( j, subimages[k]->GetIndex()[j] / m_DownsampleFactor );
      sRegion.SetSize( j , subimages[k]->GetSize()[j] / m_DownsampleFactor );
      }
    sRegion.SetIndex( 2, 0 );
    sRegion.SetSize( 2, 1 );

    typename VolumeType::Pointer sMask;

    try
      {
      std::string fileName = subimages[k]->GetMaskFilename();
      ChangePaths( fileName );

      if ( !itksys::SystemTools::FileExists( fileName.c_str(), true ) )
        {
        itkExceptionMacro( << "Missing file " << fileName << std::endl );
        }

      typedef itk::ImageFileReader< VolumeType > ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( fileName.c_str() );
      reader->Update();
      sRegion.Crop( reader->GetOutput()->GetBufferedRegion() );

      typedef itk::ExtractImageFilter<VolumeType,VolumeType> FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetInput( reader->GetOutput() );
      filter->SetExtractionRegion( sRegion );

      typedef itk::BinaryThresholdImageFilter<VolumeType,VolumeType> ThresholdType;
      typename ThresholdType::Pointer threshold = ThresholdType::New();
      threshold->SetInput( filter->GetOutput() );
      threshold->SetLowerThreshold( 1 );
      threshold->SetUpperThreshold( 255 );
      threshold->SetInsideValue( 1 );
      threshold->SetOutsideValue( 0 );

      threshold->Update();
      sMask = threshold->GetOutput();
      sMask->DisconnectPipeline();
      typename VolumeType::SpacingType spacing = this->GetSpacing();
      spacing[0] = this->GetDownsampleFactor() * subimages[k]->GetImageResolution();
      spacing[1] = spacing[0];
      sMask->SetSpacing( spacing );
      sMask->SetOrigin( this->GetOrigin() );
      sMask->SetRegions( sMask->GetBufferedRegion().GetSize() );

      typedef itk::ImageMomentsCalculator< VolumeType > CalculatorType;
      typename CalculatorType::Pointer calculator = CalculatorType::New();
      calculator->SetImage( sMask );
      calculator->Compute();

      typename CalculatorType::VectorType cog = calculator->GetCenterOfGravity();

      typedef SubImage::TransformType TransformType;
      typename TransformType::Pointer transform = TransformType::New();
      TransformType::ParametersType p( transform->GetNumberOfParameters() );
      p.Fill( 0.0 );

      for( unsigned int j = 0; j < Dimension; j++ )
        {
        p[j+3] = m_Center[j];
        }

      for( unsigned int j = 0; j < Dimension - 1; j++ )
        {
        p[j+6] = cog[j] - m_Center[j];
        }

      transform->SetParameters( p );
      subimages[k]->SetTransform( transform );

      double size = calculator->GetTotalMass();
      size *= spacing[0] * spacing[1];
      size /= 1000 * 1000;
      subimages[k]->SetTissueSize( size );

      if ( m_Verbose )
        {
        std::cout << subimages[k]->GetSpecimenTissueIndex() << ": " << size << std::endl;
        }
      
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      std::cerr << err << std::endl;
      continue;
      }
    catch( ... )
      {
      std::cerr << "SpecimenTissueIndex: " << subimages[k]->GetSpecimenTissueIndex() << std::endl;
      std::cerr << "Caught unknown exception" << std::endl;
      continue;
      }    


    }

}

/**
 * Compute transforms to align centroid of each section
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
CentroidAlignment(
MaskVolumePointer & mask,
PointType & target )
{

  unsigned int nSlices = mask->GetBufferedRegion().GetSize( 2 );
  std::vector<Transform2DPointer> transforms( nSlices );

  typedef itk::SliceImageConstIterator<VolumeType> SliceIterator;
  typedef typename SliceIterator::OutputImageType ImageType;

  SliceIterator maskIterator( mask );

  for ( unsigned int j = 0; j < nSlices && !maskIterator.IsAtEnd(); j++, ++maskIterator )
    {

    try 
      {
      Transform2DPointer t = Transform2DType::New();
      t->SetIdentity();
      transforms[j] = t;
  
      typedef itk::ImageMomentsCalculator< ImageType > CalculatorType;
      typename CalculatorType::Pointer calculator = CalculatorType::New();
      calculator->SetImage( maskIterator.GetOutput() );
      calculator->Compute();

      typename CalculatorType::VectorType cog = calculator->GetCenterOfGravity();
      Transform2DType::ParametersType p( t->GetNumberOfParameters() );
      p.Fill( 0.0 );

      for ( unsigned int k = 0; k < 2; k++ )
        {
        p[k+1] = target[k];
        p[k+3] = cog[k] - target[k];
        }

      t->SetParameters( p );

      }
     catch( itk::ExceptionObject & err )
      {
      if ( m_Verbose )
        {
        std::cerr << "sliceIndex:" << j << " " << err << std::endl;
        }
      continue;
      }
    catch( ... )
      {
      if ( m_Verbose )
        {
        std::cerr << "sliceIndex:" << j << " " << "Caught unknown exception" << std::endl;
        }
      continue;
      }    
 
    }

  this->ComposeTransforms( transforms );


}


template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
PopulateSizeVolume(
SizeVolumePointer & sizeVolume,
SizeMaskPointer & sizeMask )
{
  // assumes that ComputeVolumeMetaInformation has already been executed
  sizeVolume = SizeVolumeType::New();
  sizeMask = SizeMaskType::New();

  SizeVolumeType::SpacingType spacing;
  SizeVolumeType::PointType   origin;
  SizeVolumeType::SizeType    size;

  spacing[0] = this->GetSpacing()[2];
  spacing[1] = 1.0;

  origin[0] = this->GetOrigin()[2];
  origin[1] = 0.0;

  size[0] = this->GetSize()[2];
  size[1] = 1;

  sizeVolume->SetSpacing( spacing );
  sizeVolume->SetOrigin( origin );
  sizeVolume->SetRegions( size );
  sizeVolume->Allocate();
  sizeVolume->FillBuffer( 0.0 );

  sizeMask->SetSpacing( spacing );
  sizeMask->SetOrigin( origin );
  sizeMask->SetRegions( size );
  sizeMask->Allocate();
  sizeMask->FillBuffer( 0 );

  // iterate through the sub-images
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }

    SizeVolumeType::PointType point;
    point[0] = subimages[k]->GetSpecimenTissueIndex() * this->GetSectionThickness();
    point[1] = 0.0;

    SizeVolumeType::IndexType index;
    sizeVolume->TransformPhysicalPointToIndex( point, index );

    if ( subimages[k]->GetTissueSize() > 0 )
      {
      sizeVolume->SetPixel( index, subimages[k]->GetTissueSize() );
      sizeMask->SetPixel( index, 255 );
      }

    } // for each k

}



template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
PopulateSizeVolume(
MaskVolumePointer & mask,
SizeVolumePointer & sizeVolume,
SizeMaskPointer & sizeMask )
{

  // assumes that ComputeVolumeMetaInformation has already been executed
  sizeVolume = SizeVolumeType::New();
  sizeMask = SizeMaskType::New();

  SizeVolumeType::SpacingType spacing;
  SizeVolumeType::PointType   origin;
  SizeVolumeType::SizeType    size;

  spacing[0] = mask->GetSpacing()[2];
  spacing[1] = 1.0;

  origin[0] = mask->GetOrigin()[2];
  origin[1] = 0.0;

  size[0] = mask->GetBufferedRegion().GetSize()[2];
  size[1] = 1;

  sizeVolume->SetSpacing( spacing );
  sizeVolume->SetOrigin( origin );
  sizeVolume->SetRegions( size );
  sizeVolume->Allocate();
  sizeVolume->FillBuffer( 0.0 );

  sizeMask->SetSpacing( spacing );
  sizeMask->SetOrigin( origin );
  sizeMask->SetRegions( size );
  sizeMask->Allocate();
  sizeMask->FillBuffer( 0 );

  itk::SliceImageConstIterator< VolumeType > vit( mask );
  itk::ImageRegionIterator< SizeVolumeType > sit( sizeVolume, sizeVolume->GetBufferedRegion() );
  itk::ImageRegionIterator< SizeMaskType > mit( sizeMask, sizeMask->GetBufferedRegion() );

  vit.SetSliceIndex( 0 );
  sit.GoToBegin();
  mit.GoToBegin();

  while( !vit.IsAtEnd() )
    {

    typedef typename itk::SliceImageConstIterator< VolumeType >::OutputImageType ImageType;
    typename ImageType::ConstPointer image = vit.GetOutput();
    itk::ImageRegionConstIterator<ImageType> it( image, image->GetBufferedRegion() );
    
    double size = 0.0;
    it.GoToBegin();

    while( !it.IsAtEnd() )
      {
      if ( it.Get() )
        {
        size += 1.0;
        }
      ++it;
      }

    size *= mask->GetSpacing()[0] * mask->GetSpacing()[1];
    size /= 1000 * 1000;

    if ( size > 0.0 )
      {
      sit.Set( size );
      mit.Set( 255 );
      }

    ++vit;
    ++sit;
    ++mit;
    }

}


/**
 * Compose transform to each subimage transform
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
ComposeTransform(
const TransformType * transform )
{

  // compose input transform to each subimage transform
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }

    SubImage::TransformPointer output;
    TransformUtilities::Compose( transform, subimages[k]->GetTransform(), output );
    subimages[k]->SetTransform( output );

    }

}

/**
 * Compose transform to each subimage transform
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
ComposeTransform(
const Transform2DType * transform )
{

  // compose input transform to each subimage transform
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }

    SubImage::TransformPointer output;
    TransformUtilities::Compose( transform, subimages[k]->GetTransform(), output );
    subimages[k]->SetTransform( output );

    }

}

/**
 * Compose transform to each subimage transform
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
ComposeTransforms(
const std::vector<Transform2DPointer> & transforms )
{

  typename VolumeType::Pointer dummy = VolumeType::New();
  dummy->SetSpacing( this->GetSpacing() );
  dummy->SetOrigin( this->GetOrigin() );
  dummy->SetRegions( this->GetSize() );
  dummy->Allocate();

  // compose transform to subimage transforms
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }

    typename VolumeType::PointType point;
    typename VolumeType::IndexType index;

    point.Fill( 0.0 );
    point[2] = subimages[k]->GetSpecimenTissueIndex() * m_ImageSeries->GetSectionThickness();
    dummy->TransformPhysicalPointToIndex( point, index );

    if ( transforms[index[2]].IsNotNull() )
      {
      if ( m_Verbose )
        {
        std::cout << subimages[k]->GetSpecimenTissueIndex() << ": " << transforms[index[2]]->GetParameters() << std::endl;
        }
    
      typename SubImage::TransformPointer output;
      TransformUtilities::Compose( transforms[index[2]], subimages[k]->GetTransform() , output );
      subimages[k]->SetTransform( output );

      }

    }

}

/**
 * Compose transform to each subimage transform
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
ComposeTransforms(
const std::vector<TransformPointer> & transforms )
{

  typename VolumeType::Pointer dummy = VolumeType::New();
  dummy->SetSpacing( this->GetSpacing() );
  dummy->SetOrigin( this->GetOrigin() );
  dummy->SetRegions( this->GetSize() );
  dummy->Allocate();

  // compose transform to subimage transforms
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }

    typename VolumeType::PointType point;
    typename VolumeType::IndexType index;

    point.Fill( 0.0 );
    point[2] = subimages[k]->GetSpecimenTissueIndex() * this->GetSectionThickness();
    dummy->TransformPhysicalPointToIndex( point, index );
    
    if ( transforms[index[2]].IsNotNull() )
      {
      if ( m_Verbose )
        {
        std::cout << subimages[k]->GetSpecimenTissueIndex() << ": " << transforms[index[2]]->GetParameters() << std::endl;
        }
    
      typename SubImage::TransformPointer output;
      TransformUtilities::Compose( transforms[index[2]], subimages[k]->GetTransform() , output );
      subimages[k]->SetTransform( output );
      }

    }

}


/**
 * Backup the current transform parameters
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
BackupTransformParameters()
{

  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

   // the backup vector is sorted in the same order as the subimages
   m_BackupParameters.resize( subimages.size() );

   for( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }
    if ( !subimages[k]->GetTransform() )
      {
      m_BackupParameters[k] = ParametersType(0);
      }

    m_BackupParameters[k] = subimages[k]->GetTransform()->GetParameters();

    }

}

/**
 * revert transform parameters to the backup state
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
RevertTransformParameters()
{

  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

   for( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }
    if ( !subimages[k]->GetTransform() )
      {
      continue;
      }
    if ( subimages[k]->GetTransform()->GetNumberOfParameters() != m_BackupParameters[k].size() )
      {
      continue;
      }

    SubImage::TransformPointer t = SubImage::TransformType::New();
    t->SetParameters( m_BackupParameters[k] );
    subimages[k]->SetTransform( t );

    }

}


/**
 * write transform parameters to XML file
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
WriteTransformParametersToXML(
ImageSeriesVector & series,
const Affine3DTransformType * volumeToAtlasTransform,
const Affine3DTransformType * atlasToGridTransform,
double metric,
TiXmlNode * node )
{

  Affine3DTransformType::Pointer atlasToVolumeTransform = Affine3DTransformType::New();
  atlasToVolumeTransform->SetCenter( volumeToAtlasTransform->GetCenter() );
  volumeToAtlasTransform->GetInverse( atlasToVolumeTransform ); 
  
  atlasToVolumeTransform->Compose( atlasToGridTransform, true );
  typedef itk::idp::DecomposeAffineMatrix3D DecomposerType;
  DecomposerType::Pointer decomposer = DecomposerType::New();
  decomposer->SetInputMatrix( atlasToVolumeTransform->GetMatrix() );
  decomposer->Compute();

  for ( unsigned int s = 0; s < series.size(); s++ )
    {

    // setup a alignment3d node
    TiXmlElement * e = new TiXmlElement( "alignment3d" );

    XMLUtilities::InsertTextNode<unsigned long>( e, "image-series-id", series[s]->GetId() );

    typedef itk::AffineTransform<double,3> OutputTransformType;
    OutputTransformType::Pointer tvr = OutputTransformType::New();
    tvr->SetMatrix( volumeToAtlasTransform->GetMatrix() );
    tvr->SetOffset( volumeToAtlasTransform->GetOffset() );

    OutputTransformType::Pointer trv = OutputTransformType::New();
    tvr->GetInverse( trv );

    for ( unsigned int p = 0; p < tvr->GetParameters().GetSize(); p++ )
      {
      std::ostringstream ostr;
      ostr << "tvr-";
      ostr.width(2);
      ostr.fill( '0' );
      ostr.setf( std::ios::right );
      ostr << p;
      XMLUtilities::InsertTextNode<double>( e, ostr.str().c_str(), tvr->GetParameters()[p] );
      }
    for ( unsigned int p = 0; p < trv->GetParameters().GetSize(); p++ )
      {
      std::ostringstream ostr;
      ostr << "trv-";
      ostr.width(2);
      ostr.fill( '0' );
      ostr.setf( std::ios::right );
      ostr << p;
      XMLUtilities::InsertTextNode<double>( e, ostr.str().c_str(), trv->GetParameters()[p] );
      }

    XMLUtilities::InsertTextNode<double>( e, "metric", metric );
    XMLUtilities::InsertTextNode<double>( e, "scale-x", decomposer->GetScale()[0] );
    XMLUtilities::InsertTextNode<double>( e, "scale-y", decomposer->GetScale()[1] );
    XMLUtilities::InsertTextNode<double>( e, "scale-z", decomposer->GetScale()[2] );
    if (std::isnan(decomposer->GetRotation()[0])||std::isnan(decomposer->GetRotation()[1])||std::isnan(decomposer->GetRotation()[2])) {
        XMLUtilities::InsertTextNode<double>( e, "rotation-x", 0.0 ); //if nan, just write 0
        XMLUtilities::InsertTextNode<double>( e, "rotation-y", 0.0 );
        XMLUtilities::InsertTextNode<double>( e, "rotation-z", 0.0 );
    } else {
        XMLUtilities::InsertTextNode<double>( e, "rotation-x", decomposer->GetRotation()[0] ); 
        XMLUtilities::InsertTextNode<double>( e, "rotation-y", decomposer->GetRotation()[1] );
        XMLUtilities::InsertTextNode<double>( e, "rotation-z", decomposer->GetRotation()[2] );    
    }
    XMLUtilities::InsertTextNode<double>( e, "skew-x", decomposer->GetShear()[0] );
    XMLUtilities::InsertTextNode<double>( e, "skew-y", decomposer->GetShear()[1] );
    XMLUtilities::InsertTextNode<double>( e, "skew-z", decomposer->GetShear()[2] );

    node->LinkEndChild( e );

    }


}


/**
 * write transform parameters to XML file
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
WriteTransformParametersToXML(
ImageSeriesVector & series,
TiXmlNode * node )
{
  
  for ( unsigned int s = 0; s < series.size(); s++ )
    {
    for ( unsigned int k = 0; k < series[s]->GetSubImages().size(); k++ )
      {

      // set up alignment2d node
      TiXmlElement * e = new TiXmlElement( "alignment2d" );

      XMLUtilities::InsertTextNode<unsigned long>( e, "sub-image-id", series[s]->GetSubImages()[k]->GetId() );

      typedef itk::AffineTransform<double,2> OutputTransformType;
      OutputTransformType::Pointer tvs = OutputTransformType::New();

      OutputTransformType::MatrixType matrix;
      OutputTransformType::OutputVectorType offset;

      matrix[0][0] = series[s]->GetSubImages()[k]->GetTransform()->GetMatrix()[0][0];
      matrix[0][1] = series[s]->GetSubImages()[k]->GetTransform()->GetMatrix()[0][1];
      matrix[1][0] = series[s]->GetSubImages()[k]->GetTransform()->GetMatrix()[1][0];
      matrix[1][1] = series[s]->GetSubImages()[k]->GetTransform()->GetMatrix()[1][1];
      offset[0] = series[s]->GetSubImages()[k]->GetTransform()->GetOffset()[0];
      offset[1] = series[s]->GetSubImages()[k]->GetTransform()->GetOffset()[1];

      tvs->SetMatrix( matrix );
      tvs->SetOffset( offset );

      tvs->Scale( 1.0 / series[s]->GetSubImages()[k]->GetImageResolution(), false );

      OutputTransformType::Pointer tsv = OutputTransformType::New();
      tvs->GetInverse( tsv );

   
      for ( unsigned int p = 0; p < tvs->GetParameters().GetSize(); p++ )
        {
        std::ostringstream ostr;
        ostr << "tvs-";
        ostr.width(2);
        ostr.fill( '0' );
        ostr.setf( std::ios::right );
        ostr << p;
        XMLUtilities::InsertTextNode<double>( e, ostr.str().c_str(), tvs->GetParameters()[p] );
        }
      for ( unsigned int p = 0; p < tsv->GetParameters().GetSize(); p++ )
        {
        std::ostringstream ostr;
        ostr << "tsv-";
        ostr.width(2);
        ostr.fill( '0' );
        ostr.setf( std::ios::right );
        ostr << p;
        XMLUtilities::InsertTextNode<double>( e, ostr.str().c_str(), tsv->GetParameters()[p] );
        }

      XMLUtilities::InsertTextNode<double>( e, "metric", series[s]->GetSubImages()[k]->GetMetric() );

      node->LinkEndChild( e );

     }
    }

}

/**
 * load transform parameters from XML File
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
ReadTransformParametersFromXML(
TiXmlNode * node )
{

  TiXmlNode * n = NULL;

  while ( ( n = node->IterateChildren( "alignment2d", n ) ) )
    {
    this->ReadSubImageTransformParametersFromXML( n );
    }

}


template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
ReadSubImageTransformParametersFromXML(
TiXmlNode * node )
{

  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  TiXmlNode * n = NULL;
  TiXmlElement * e;

  ParametersType ip( 6 );
  ip.Fill( 0.0 );

  unsigned long sid = 0;

  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();

      if ( parameter.compare( "sub-image-id" ) == 0 )
        {
        sid = atoi( e->GetText() );
        }
      else if ( parameter.find( "tvs-" ) != std::string::npos )
        {
        unsigned int index = atoi( parameter.substr(4,2).c_str() );
        ip[index] = atof( e->GetText() );
        }

      } // if (e)
    } // while (n)


  // search for the right subimage
  for( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetId() == sid )
      {

      typedef itk::AffineTransform<double,2> AffineTransformType;
      AffineTransformType::Pointer it = AffineTransformType::New();
      it->SetParameters( ip );
      it->Scale( subimages[k]->GetImageResolution(), false );
      ip = it->GetParameters();

      double rotation = -1.0 * asin( ip[1] );

      ParametersType op( 9 );
      op.Fill( 0.0 );
      op[2] = rotation;
      op[6] = ip[4];
      op[7] = ip[5];

      SubImage::TransformPointer ot = SubImage::TransformType::New();
      ot->SetParameters( op );

      std::cout << subimages[k]->GetSpecimenTissueIndex() << ": " << op << std::endl;
      subimages[k]->SetTransform( ot );
      }
    }

}


/**
 * load transform parameters from XML file
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
ReadTransformParametersFromXML(
TiXmlNode * node,
Affine3DTransformPointer & transform )
{

  TiXmlNode * al3d = NULL;

  while ( ( al3d = node->IterateChildren( "alignment3d", al3d ) ) )
    {
    
    TiXmlNode * n = NULL;
    TiXmlElement * e;

    ParametersType op( 15 );
    op.Fill( 0.0 );

    unsigned long iid = 0;

    while ( ( n = al3d->IterateChildren( n ) ) )
      {

      e = n->ToElement();

      if ( e )
        {
        std::string parameter = e->Value();

        if ( parameter.compare( "image-series-id" ) == 0 )
          {
          iid = atoi( e->GetText() );
          }
        else if ( parameter.find( "tvr-" ) != std::string::npos )
          {
          unsigned int index = atoi( parameter.substr(4,2).c_str() );
          if ( index < 9 )
            {
            op[index] = atof( e->GetText() );
            }
          else
            {
            op[index+3] = atof( e->GetText() );
            }
          }

        } // if (e)
      } // while (n)


    if ( m_ImageSeries->GetId() && m_ImageSeries->GetId() != iid )
      {
      continue;
      }

    transform = Affine3DTransformType::New();
    transform->SetParameters( op );
    std::cout << iid << ": " << op << std::endl;
     
    } // iterate alignment3d

}


/**
 * Set the metric for each sub_image from input vector
 */
template<class TPixel>
void
ImageSeriesUtilities<TPixel>::
SetSubImageMetric(
const std::vector<double> & metric )
{

  typename VolumeType::Pointer dummy = VolumeType::New();
  dummy->SetSpacing( this->GetSpacing() );
  dummy->SetOrigin( this->GetOrigin() );
  dummy->SetRegions( this->GetSize() );
  dummy->Allocate();

  // compose transform to subimage transforms
  ImageSeries::SubImageArray subimages = m_ImageSeries->GetSubImages();

  for ( unsigned int k = 0; k < subimages.size(); k++ )
    {
    if ( subimages[k]->GetFailed() )
      {
      continue;
      }

    typename VolumeType::PointType point;
    typename VolumeType::IndexType index;

    point.Fill( 0.0 );
    point[2] = subimages[k]->GetSpecimenTissueIndex() * m_ImageSeries->GetSectionThickness();
    dummy->TransformPhysicalPointToIndex( point, index );

    subimages[k]->SetMetric( metric[index[2]] );

    if ( m_Verbose )
      {
        std::cout << subimages[k]->GetSpecimenTissueIndex() << ": " << metric[index[2]] << std::endl;
      }

    }

}


/**
 * Read in high resolution image
 */
template<class TPixel>
template <class TImage>
void
ImageSeriesUtilities<TPixel>::
ReadImage(
const SubImage * subimage,
typename TImage::Pointer & image,
int reduce,
int channel )
{
  std::string fname = subimage->GetJp2Filename();
  ChangePaths( fname );

  if ( !itksys::SystemTools::FileExists( fname.c_str(), true ) )
    {
    itkExceptionMacro( << "Missing file " << fname << std::endl );
    }

  //std::cout << "Reading: " << fname << std::endl;

  typedef ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( fname.c_str() );

  typename TImage::RegionType region;
  region.SetIndex( subimage->GetIndex() );
  region.SetSize( subimage->GetSize() );

  JP2ImageIO::Pointer io = JP2ImageIO::New();
  io->SetReduce( reduce );
  io->SetRegionOfInterest( region );

  if (channel >= 0)
      io->SetReadSingleComponentIndex(channel);
  
  reader->SetImageIO( io );

  try
    {
    bool ok = io->CanReadFile( fname.c_str() );

    if ( !ok )
      {
      std::cout << fname << "is not a valid jp2 file" << std::endl;
      image = 0;
      return;
      }

    reader->Update();
    image = reader->GetOutput();
    image->DisconnectPipeline();
    reader = 0;

    double s = pow( 2.0, reduce ) * subimage->GetImageResolution();
    typename TImage::SpacingType spacing;
    spacing.Fill( s );
    image->SetSpacing( spacing );
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    throw err;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    image = 0;
    return;
    }

}



/**
 * Read in high resolution expresssion image
 */
template<class TPixel>
template <class TImage>
void
ImageSeriesUtilities<TPixel>::
ReadExpressionImage(
const SubImage * subimage,
typename TImage::Pointer & image,
int reduce,
int channel )
{

  std::string fname = subimage->GetExpressionJp2Filename();
  ChangePaths( fname );

  std::cout << "Reading: " << fname << std::endl;

  typedef ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( fname.c_str() );

  typename RGBImageType::RegionType region;
  region.SetIndex( subimage->GetIndex() );
  region.SetSize( subimage->GetSize() );

  JP2ImageIO::Pointer io = JP2ImageIO::New();
  io->SetReduce( reduce );
  io->SetRegionOfInterest( region );

  if ( channel >= 0 )
      io->SetReadSingleComponentIndex(channel);
  
  reader->SetImageIO( io );

  try
    {
    bool ok = io->CanReadFile( fname.c_str() );

    if ( !ok )
      {
      std::cout << fname << "is not a valid jp2 file" << std::endl;
      image = 0;
      return;
      }

    reader->Update();
    image = reader->GetOutput();
    image->DisconnectPipeline();
    reader = 0;

    double s = pow( 2.0, reduce ) * subimage->GetImageResolution();
    typename ImageType::SpacingType spacing;
    spacing.Fill( s );
    image->SetSpacing( spacing );

    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    image = 0;
    return;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    image = 0;
    return;
    }    


}

/**
 * Read in high resolution projection mask image
 */
template<class TPixel>
template <class TImage>
void
ImageSeriesUtilities<TPixel>::
ReadProjectionMaskImage(
const SubImage * subimage,
typename TImage::Pointer & image,
int reduce,
int channel )
{
  std::string fname = subimage->GetProjectionMaskJp2Filename();

  ChangePaths( fname );

  std::cout << "Reading: " << fname << std::endl;

  typedef ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( fname.c_str() );

  typename RGBImageType::RegionType region;
  region.SetIndex( subimage->GetIndex() );
  region.SetSize( subimage->GetSize() );

  JP2ImageIO::Pointer io = JP2ImageIO::New();
  io->SetReduce( reduce );
  io->SetRegionOfInterest( region );

  if ( channel >= 0 )
      io->SetReadSingleComponentIndex(channel);
  
  reader->SetImageIO( io );

  try
    {
    bool ok = io->CanReadFile( fname.c_str() );

    if ( !ok )
      {
      std::cout << fname << "is not a valid jp2 file" << std::endl;
      image = 0;
      return;
      }

    reader->Update();
    image = reader->GetOutput();
    image->DisconnectPipeline();
    reader = 0;

    double s = pow( 2.0, reduce ) * subimage->GetImageResolution();
    typename ImageType::SpacingType spacing;
    spacing.Fill( s );
    image->SetSpacing( spacing );

    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    image = 0;
    return;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    image = 0;
    return;
    }    


}


/**
 * Resample input image into volume using current transform
 */
template<class TPixel>
template <class TOutputPixel>
void
ImageSeriesUtilities<TPixel>::
ResampleIntoVolume(
const SubImage * subimage,
const Image<TOutputPixel,2> * image,
Image<TOutputPixel,3> * volume,
const char * interpolationType )
{

  // wrap input as 3D volume with 1 slice
  typedef ImportImageFilter< TOutputPixel, 3 > ImporterType;
  typedef typename ImporterType::OutputImageType  VolumeType;

  typename ImporterType::Pointer importer = ImporterType::New();

  typename VolumeType::SpacingType spacing;
  typename VolumeType::PointType   origin;
  typename VolumeType::IndexType   index;
  typename VolumeType::SizeType    size;

  for ( int j = 0; j < 2; j++ )
    {
    spacing[j] = image->GetSpacing()[j];
    origin[j] = image->GetOrigin()[j];
    index[j] = image->GetBufferedRegion().GetIndex()[j];
    size[j] = image->GetBufferedRegion().GetSize()[j];
    }

  spacing[2] = m_Spacing[2];
  origin[2] = subimage->GetSpecimenTissueIndex() * this->GetImageSeries()->GetSectionThickness();
  index[2] = 0;
  size[2] = 1;
 
  typename VolumeType::RegionType region;
  region.SetIndex( index );
  region.SetSize( size );
  
  importer->SetSpacing( spacing );
  importer->SetOrigin( origin );
  importer->SetRegion( region );
  importer->SetImportPointer( const_cast<TOutputPixel *>( image->GetBufferPointer() ), region.GetNumberOfPixels(), false );
  
  // compute reference output size
  spacing = volume->GetSpacing();
  origin = volume->GetOrigin();
  origin[2] = subimage->GetSpecimenTissueIndex() * this->GetImageSeries()->GetSectionThickness();
  index = volume->GetBufferedRegion().GetIndex();
  index[2] = 0;
  size = volume->GetBufferedRegion().GetSize();
  size[2] = 1;
  region.SetIndex( index );
  region.SetSize( size );

  typename VolumeType::Pointer ref = VolumeType::New();
  ref->SetSpacing( spacing );
  ref->SetOrigin( origin );
  ref->SetRegions( region );
  ref->Allocate();
  ref->FillBuffer( 0 );

  // resample image
  typename VolumeType::Pointer input = importer->GetOutput();
  typename VolumeType::Pointer output;
  ResampleImage<VolumeType>( input, ref, subimage->GetTransform(), output, 0, interpolationType );

  typename VolumeType::PointType d;
  d.Fill( 0.0 );
  d[2] = subimage->GetSpecimenTissueIndex() * this->GetImageSeries()->GetSectionThickness();

  volume->TransformPhysicalPointToIndex( d, index );
  
  region.SetSize( size );
  region.SetIndex( index );

  // copy into volume
  typename VolumeType::Pointer vptr = volume;
  CopyRegion<VolumeType>( output, output->GetBufferedRegion(), vptr, region );

}


template<class TPixel>
template <class TImage>
void
ImageSeriesUtilities<TPixel>::
ReadImageFromFile(
    const std::string& filename, 
    typename TImage::Pointer & image )
{
  typedef ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );  

  try
    {
    reader->Update();
    image = reader->GetOutput();
    image->DisconnectPipeline();
    reader = 0;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    image = 0;
    return;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    image = 0;
    return;
    }    
}


} // end namespace idp
} //end namespace itk

#endif
