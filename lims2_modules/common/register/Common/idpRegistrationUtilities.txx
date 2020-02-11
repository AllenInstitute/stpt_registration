/*=========================================================================

  idpRegistrationUtilities.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpRegistrationUtilities_txx
#define __idpRegistrationUtilities_txx

#include "idpRegistrationUtilities.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkTranslationTransform.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImportImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkMeanProjectionImageFilter.h"
#include "itkMedianProjectionImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkSumProjectionImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkChangePixelSpacingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkMaskImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkNotImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkLog10ImageFilter.h"
#include "itkLogImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkConstantPadImageFilter.h"


namespace itk
{
namespace idp
{

/*--------------------------
 * Read image in from file
 * --------------------------
 */
template <typename ImageType>
void
ReadImage( const char * filename,
           typename ImageType::Pointer & image ) throw (ExceptionObject)
{
  typedef ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename );

  try
    {
    reader->Update();
    image = reader->GetOutput();
    image->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

};

/*--------------------------
 * Write image in from file
 * --------------------------
 */
template <typename ImageType>
void
WriteImage( const char * filename,
            const ImageType * image,
            bool compression ) throw (ExceptionObject)
{
  typedef ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename );

  if (compression)
    {
    writer->UseCompressionOn();
    }

  writer->SetInput( image );

  try
    {
    writer->Update();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

};

/*---------------------------
 * Resample image using input transform
 * --------------------------
 */
template <typename ImageType>
void
ResampleImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & ref,
const Transform< double, ImageType::ImageDimension, ImageType::ImageDimension> * trans,
typename ImageType::Pointer & output,
typename ImageType::PixelType pad,
const char * interpolatorType
 )
throw( ExceptionObject )
{

  typedef ResampleImageFilter< ImageType, ImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetTransform( trans );
  filter->SetInput( input );
  filter->SetSize( ref->GetBufferedRegion().GetSize() );
  filter->SetOutputOrigin( ref->GetOrigin() );
  filter->SetOutputSpacing( ref->GetSpacing() );
  filter->SetDefaultPixelValue( pad );

  typedef InterpolateImageFunction<ImageType,double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;

  std::string str = interpolatorType;

#define _CreateInterpolator( type ) \
  ( str.compare( #type ) == 0 ) \
    { \
    typedef type##InterpolateImageFunction<ImageType,double> IType; \
    typename IType::Pointer iptr = IType::New(); \
    interpolator = iptr; \
    }

  if ( !str.empty() )
    {
    if _CreateInterpolator( Linear )
    else if _CreateInterpolator( NearestNeighbor )
    else if _CreateInterpolator( BSpline )
    }

#undef _CreateInterpolator

  if ( interpolator.IsNotNull() )
    {
    filter->SetInterpolator( interpolator );
    }

  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

/*---------------------------
 * Resample vectorimage using input transform
 * --------------------------
 */
template <typename ImageType>
void
VectorResampleImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & ref,
const Transform< double, ImageType::ImageDimension, ImageType::ImageDimension> * trans,
typename ImageType::Pointer & output,
typename ImageType::PixelType pad,
const char * interpolatorType
 )
throw( ExceptionObject )
{

  typedef VectorResampleImageFilter< ImageType, ImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetTransform( trans );
  filter->SetInput( input );
  filter->SetSize( ref->GetBufferedRegion().GetSize() );
  filter->SetOutputOrigin( ref->GetOrigin() );
  filter->SetOutputSpacing( ref->GetSpacing() );
  filter->SetDefaultPixelValue( pad );

  typedef VectorInterpolateImageFunction<ImageType,double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;

  std::string str = interpolatorType;

#define _CreateInterpolator( type ) \
  ( str.compare( #type ) == 0 ) \
    { \
    typedef type##InterpolateImageFunction<ImageType,double> IType; \
    typename IType::Pointer iptr = IType::New(); \
    interpolator = iptr; \
    }

  if ( !str.empty() )
    {
    if _CreateInterpolator( VectorLinear )
    else if _CreateInterpolator( VectorNearestNeighbor )    
    }

#undef _CreateInterpolator

  if ( interpolator.IsNotNull() )
    {
    filter->SetInterpolator( interpolator );
    }

  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

/*--------------------------
 * Compute correlation between two images
 * --------------------------
 */
template <typename ImageType, typename MaskImageType>
void
ComputeCorrelation( const ImageType * image1,
                    const ImageType * image2,
                    const MaskImageType * mask1,
                    const MaskImageType * mask2,
                    double & value ) throw (ExceptionObject)
{
  typedef NormalizedCorrelationImageToImageMetric< ImageType, ImageType > MetricType;
  typename MetricType::Pointer metric = MetricType::New();

  metric->SetFixedImage( image1 );
  metric->SetMovingImage( image2 );

  typedef TranslationTransform< double, ImageType::ImageDimension > TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  metric->SetTransform( transform );
  typename TransformType::ParametersType parameters( metric->GetNumberOfParameters() );
  parameters.Fill( 0.0 );

  typedef LinearInterpolateImageFunction< ImageType > InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  metric->SetInterpolator( interpolator );

  typedef ImageMaskSpatialObject<ImageType::ImageDimension> MaskType;
  typename MaskType::Pointer maskObject1 = MaskType::New();
  typename MaskType::Pointer maskObject2 = MaskType::New();

  if ( mask1 )
    {
    maskObject1->SetImage( mask1 );
    metric->SetFixedImageMask( maskObject1 );
    }

  if ( mask2 )
    {
    maskObject2->SetImage( mask2 );
    metric->SetMovingImageMask( maskObject2 );
    }

  try
    {
    metric->ComputeGradientOff();
    metric->SubtractMeanOn();
    metric->SetFixedImageRegion( image1->GetBufferedRegion() );
    metric->Initialize();

    value = metric->GetValue( parameters );
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }
 
}


/*---------------------------
 * Copy region from one image to another
 * --------------------------
 */
template <typename ImageType>
void
CopyRegion(
typename ImageType::Pointer & fromImage,
const typename ImageType::RegionType & fromRegion,
typename ImageType::Pointer & toImage,
const typename ImageType::RegionType & toRegion
)
throw( ExceptionObject )
{
  if ( fromRegion.GetSize() != toRegion.GetSize() )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "From and to regions are not the same size", ITK_LOCATION );

    }

  try
    {
    ImageRegionConstIterator<ImageType> fi( fromImage, fromRegion );
    ImageRegionIterator<ImageType> ti( toImage,toRegion );

    fi.GoToBegin();
    ti.GoToBegin();

    while( !fi.IsAtEnd() )
      {
      ti.Set( fi.Get() );
      ++fi;
      ++ti;
      }
   
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}


/*---------------------------
 * Compute minimum projection image
 * --------------------------
 */
template <typename ImageType>
void
MinimumProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int dim ) throw (ExceptionObject)
{
  typedef MinimumProjectionImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetProjectionDimension( dim );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}


/*---------------------------
 * Compute minimum projection image
 * --------------------------
 */
template <typename ImageType>
void
MaximumProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int dim ) throw (ExceptionObject)
{
  typedef MaximumProjectionImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetProjectionDimension( dim );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

/*---------------------------
 * Compute average projection image
 * --------------------------
 */
template <typename ImageType>
void
AverageProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int dim ) throw (ExceptionObject)
{
  typedef MeanProjectionImageFilter<ImageType,ImageType,double> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetProjectionDimension( dim );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

/*---------------------------
 * Compute average projection image
 * --------------------------
 */
template <typename ImageType>
void
MedianProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int dim ) throw (ExceptionObject)
{
  typedef MedianProjectionImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetProjectionDimension( dim );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to compute the average projection image
template <typename ImageType, typename MaskImageType>
void
MaskedAverageProjection(
typename ImageType::Pointer & input,
typename MaskImageType::Pointer & mask,
typename ImageType::Pointer & output,
unsigned int dim ) throw( ExceptionObject )
{

  typedef Image<double, ImageType::ImageDimension> DoubleImageType;
  typedef BinaryThresholdImageFilter<MaskImageType,ImageType> ThresholdType;
  typedef SumProjectionImageFilter<ImageType,DoubleImageType> FilterType;
  typedef ShiftScaleImageFilter<DoubleImageType,ImageType> ScalerType;
  typedef MinimumMaximumImageFilter<DoubleImageType> CalculatorType;

  typename FilterType::Pointer filter1 = FilterType::New();
  filter1->SetInput( input );
  filter1->SetProjectionDimension( dim );

  typename ThresholdType::Pointer threshold = ThresholdType::New();
  threshold->SetInput( mask );
  threshold->SetLowerThreshold( 1 );
  threshold->SetUpperThreshold( NumericTraits<typename ImageType::PixelType>::max() );
  threshold->SetInsideValue( 1 );
  threshold->SetOutsideValue( 0 );

  typename FilterType::Pointer filter2 = FilterType::New();
  filter2->SetInput( threshold->GetOutput() );
  filter2->SetProjectionDimension( dim );

  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetInput( filter2->GetOutput() );

  typename ScalerType::Pointer scaler = ScalerType::New();
  scaler->SetInput( filter1->GetOutput() );
   
  try
    {
    calculator->Update();
    //std::cout << calculator->GetMaximum() << std::endl;
    scaler->SetScale( 1.0 / calculator->GetMaximum() );
    scaler->Update();
    output = scaler->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

/*---------------------------
 * Invert image intensity
 * --------------------------
 */
template <typename ImageType>
void
InvertIntensity(
const ImageType  * input,
typename ImageType::Pointer & output ) throw (ExceptionObject)
{
  typedef InvertIntensityImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

/*---------------------------
 * Extract one slice from a volume
 * --------------------------
 */
template <typename InputImageType, typename OutputImageType>
void
ExtractSlice(
typename InputImageType::Pointer & input,
typename OutputImageType::Pointer & output,
unsigned int slice,
unsigned int dimension ) throw( ExceptionObject )
{
  typedef ExtractImageFilter<InputImageType,OutputImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  
  typename InputImageType::RegionType region;
  region = input->GetBufferedRegion();
  region.SetIndex( dimension, slice );
  region.SetSize( dimension, 0 );

  filter->SetInput( input );
  filter->SetExtractionRegion( region );

  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to pad an image
template <typename ImageType>
void
PadImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
typename ImageType::SizeType  upperExtendSize,
typename ImageType::SizeType  lowerExtendSize,
typename ImageType::PixelType  padValue        
)
throw( ExceptionObject )
{
  //typename ImageType::SizeType lowerExtendRegion;
  //lowerExtendRegion.Fill( lowerExtendSize );

  //typename ImageType::SizeType upperExtendRegion;
  //upperExtendRegion.Fill( upperExtendSize );

  typedef ConstantPadImageFilter< ImageType, ImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetPadLowerBound(lowerExtendSize);
  filter->SetPadUpperBound(upperExtendSize);
  filter->SetConstant( padValue );
  
  filter->Update();
  
  output = filter->GetOutput();
    
}


/*---------------------------
 * Extract region from input volume
 * --------------------------
 */
template <typename ImageType>
void
ExtractRegion(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
typename ImageType::RegionType region ) throw( ExceptionObject )
{
  typedef ExtractImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  region.Crop( input->GetBufferedRegion() );
  
  filter->SetInput( input );
  filter->SetExtractionRegion( region );

  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}


/*---------------------------
 * Change the image pixel spacing
 * --------------------------
 */
template <typename ImageType>
void
ChangePixelSpacing(
typename ImageType::Pointer & input,
const FixedArray<double,ImageType::ImageDimension> & spacing,
const FixedArray<double,ImageType::ImageDimension> & variance,
typename ImageType::Pointer & output,
const char * interpolatorType,
typename ImageType::PixelType pad
) throw( ExceptionObject )
{
  typedef ChangePixelSpacingImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  typedef typename FilterType::InternalImageType InternalImageType;
  typedef InterpolateImageFunction<InternalImageType,double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;

  std::string str = interpolatorType;

#define _CreateInterpolator( type ) \
  ( str.compare( #type ) == 0 ) \
    { \
    typedef type##InterpolateImageFunction<InternalImageType,double> IType; \
    typename IType::Pointer iptr = IType::New(); \
    interpolator = iptr; \
    }

  if ( !str.empty() )
    {
    if _CreateInterpolator( Linear )
    else if _CreateInterpolator( NearestNeighbor )
    else if _CreateInterpolator( BSpline )
    }

#undef _CreateInterpolator

  if ( interpolator.IsNotNull() )
    {
    filter->SetInterpolator( interpolator );
    }

  filter->SetInput( input );
  filter->SetSpacing( spacing );
  filter->SetSmoothingVariance( variance );
  filter->SetDefaultPixelValue( pad );

  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}


// helper function to threshold the image
template <typename ImageType>
void
BinaryThreshold(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
typename ImageType::PixelType lower,
typename ImageType::PixelType upper,
typename ImageType::PixelType inside,
typename ImageType::PixelType outside
 ) throw( ExceptionObject )
{  
  typedef BinaryThresholdImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetLowerThreshold( lower );
  filter->SetUpperThreshold( upper );
  filter->SetInsideValue( inside );
  filter->SetOutsideValue( outside );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }
}

// helper function to cast image type to another
template <typename InputImageType, typename OutputImageType>
void
CastImage(
typename InputImageType::Pointer & input,
typename OutputImageType::Pointer & output
) throw( ExceptionObject )
{
  typedef CastImageFilter<InputImageType,OutputImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}


// helper function to write transform to file
template< typename TransformType >
void WriteTransform(
const char * filename,
typename TransformType::Pointer & transform
)throw( ExceptionObject )
{
  typedef TransformFileWriter WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( transform );
  writer->SetFileName( filename );

  try
    {
    writer->Update();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to read transform to file
template< typename TransformType >
void ReadTransform(
const char * filename,
typename TransformType::Pointer & transform
)throw( ExceptionObject )
{
  typedef TransformFileReader ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename );

  try
    {
    reader->Update();

    typedef ReaderType::TransformListType TransformListType;
    TransformListType * tlist = reader->GetTransformList();
    TransformListType::const_iterator it = tlist->begin();

    transform = static_cast<TransformType*>((*it).GetPointer());

    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}


// helper function to fill a specifiy region
template <typename ImageType>
void
FillRegion(
typename ImageType::Pointer & ref,
typename ImageType::Pointer & output,
typename ImageType::RegionType region,
typename ImageType::PixelType foreground,
typename ImageType::PixelType background
)
throw( ExceptionObject )
{
 
  try
    {
    output = ImageType::New();
    output->SetRegions( ref->GetBufferedRegion() );
    output->SetSpacing( ref->GetSpacing() );
    output->SetOrigin( ref->GetOrigin() );
    output->Allocate();
    output->FillBuffer( background );
    region.Crop( output->GetBufferedRegion() );

    typedef ImageRegionIterator<ImageType> Iterator;
    Iterator it( output, region );
  
    while( !it.IsAtEnd() )
      {
      it.Set( foreground );
      ++it;
      }

    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function make a default image
template <typename ImageType, typename RefType >
void
MakeImage(
typename RefType::Pointer & ref,
typename ImageType::Pointer & output,
typename ImageType::PixelType background
)
throw( ExceptionObject )
{
 
  try
    {
    output = ImageType::New();
    output->SetRegions( ref->GetBufferedRegion() );
    output->SetSpacing( ref->GetSpacing() );
    output->SetOrigin( ref->GetOrigin() );
    output->Allocate();
    output->FillBuffer( background );
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

template <typename ImageType, typename MaskType>
void
MaskImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & mask,
typename ImageType::Pointer & output,
typename ImageType::PixelType outsideValue
) throw( ExceptionObject )
{
  typedef MaskImageFilter<ImageType,MaskType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( input );
  filter->SetInput2( mask );
  filter->SetOutsideValue( outsideValue );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }
}

template <typename ImageType>
void
GaussianSmoothImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
const double * variance
) throw( ExceptionObject )
{

  typedef DiscreteGaussianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  // variance is specified in spacing units.  find the maximum variance 
  // in any dimension in pixel units.
  float max_pixel_variance = 0;
  typename FilterType::ArrayType v;

  for ( int j = 0; j < ImageType::ImageDimension; j++ )
    {
    v[j] = variance[j];

    float pixel_variance = v[j] / input->GetSpacing()[j];
    if (pixel_variance > max_pixel_variance)
        max_pixel_variance = pixel_variance;
    }

  // the blurring kernel shouldn't need more than 6 sigma worth of pixels
  int kernel_width = 6.0 * sqrtf(max_pixel_variance);
  if (kernel_width > 32) {
      filter->SetMaximumKernelWidth(kernel_width);
  }

  filter->SetVariance( v );
  filter->SetUseImageSpacingOn();
  filter->SetInput( input );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }


}

// helper function to divide one image by another
template <typename ImageType>
void
DivideImage(
typename ImageType::Pointer & input1,
typename ImageType::Pointer & input2,
typename ImageType::Pointer & output
)
throw( ExceptionObject )
{

  typedef DivideImageFilter<ImageType,ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput1( input1 );
  filter->SetInput2( input2 );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to add one image by another
template <typename ImageType>
void
AddImage(
typename ImageType::Pointer & input1,
typename ImageType::Pointer & input2,
typename ImageType::Pointer & output,
bool inPlace
)
throw( ExceptionObject )
{

  typedef AddImageFilter<ImageType,ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput1( input1 );
  filter->SetInput2( input2 );
  filter->SetInPlace( inPlace );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to compute image statistics
template <typename ImageType>
void
ImageStatistics(
typename ImageType::Pointer & input,
double & minimum,
double & maximum,
double & sum,
double & mean,
double & variance,
double & sigma
)
throw( ExceptionObject )
{

  typedef StatisticsImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  
  try
    {
    filter->Update();
    minimum = filter->GetMinimum();
    maximum = filter->GetMaximum();
    sum = filter->GetSum();
    mean = filter->GetMean();
    variance = filter->GetVariance();
    sigma = filter->GetSigma();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to not a image
template <typename ImageType>
void
NotImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output
)
throw( ExceptionObject )
{

  typedef NotImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to shift and scale an image
template <typename InputImageType, typename OutputImageType>
void
ShiftScale(
typename InputImageType::Pointer & input,
typename OutputImageType::Pointer & output,
double shift,
double scale
)throw( ExceptionObject )
{

  typedef ShiftScaleImageFilter<InputImageType,OutputImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  filter->SetShift( shift );
  filter->SetScale( scale );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }


}

// select on component out of a vector
template <typename InputImageType, typename OutputImageType>
void
VectorIndexSelection( 
typename InputImageType::Pointer & input,
unsigned int index,
typename OutputImageType::Pointer & output
) throw( ExceptionObject )
{
  typedef VectorIndexSelectionCastImageFilter<InputImageType,OutputImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetIndex( index );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to log 10 an image
template <typename ImageType>
void
Log10Image(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output
)
throw( ExceptionObject )
{

  typedef Log10ImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to log 10 an image
template <typename ImageType>
void
LogImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output
)
throw( ExceptionObject )
{

  typedef LogImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  
  try
    {
    filter->Update();
    output = filter->GetOutput();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

// helper function to threshold below
template <typename ImageType>
void
ThresholdBelow(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
double thresh,
double outsideValue
)
throw( ExceptionObject )
{

  typedef ThresholdImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  filter->ThresholdBelow( thresh );
  filter->SetOutsideValue( outsideValue );

  try
    {
    filter->Update();
    output = filter->GetOutput();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }


}

// helper function to threshold above
template <typename ImageType>
void
ThresholdAbove(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
double thresh,
double outsideValue
)
throw( ExceptionObject )
{

  typedef ThresholdImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  filter->ThresholdAbove( thresh );
  filter->SetOutsideValue( outsideValue );

  try
    {
    filter->Update();
    output = filter->GetOutput();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }

}

/*---------------------------
 * Resample image using deformation field
 * --------------------------
 */
template <typename ImageType, typename RefImageType, typename DeformationFieldType>
void ResampleImageByDfmfld(
typename ImageType::Pointer & input,
typename RefImageType::Pointer & ref,
typename DeformationFieldType::Pointer & dfmfld,
typename ImageType::Pointer & output,
const char * interpolatorType
 )
throw( ExceptionObject )
{  
  // warp the image 
  typedef itk::WarpImageFilter< ImageType, ImageType, DeformationFieldType  >  WarpImageFilterType;
 
  typename WarpImageFilterType::Pointer warpImageFilter =  WarpImageFilterType::New();
 
  typedef InterpolateImageFunction<ImageType,double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;

  std::string str = interpolatorType;

#define _CreateInterpolator( type ) \
  ( str.compare( #type ) == 0 ) \
    { \
    typedef type##InterpolateImageFunction<ImageType,double> IType; \
    typename IType::Pointer iptr = IType::New(); \
    interpolator = iptr; \
    }

  if ( !str.empty() )
    {
    if _CreateInterpolator( Linear )
    else if _CreateInterpolator( NearestNeighbor )
    else if _CreateInterpolator( BSpline )
    }

#undef _CreateInterpolator

  if ( interpolator.IsNotNull() )
  {
    warpImageFilter->SetInterpolator( interpolator );
	warpImageFilter->SetOutputSpacing( ref->GetSpacing() );
    warpImageFilter->SetOutputOrigin(  ref->GetOrigin() );
    warpImageFilter->SetOutputSize( ref->GetBufferedRegion().GetSize() );    
    warpImageFilter->SetDeformationField( dfmfld );
    warpImageFilter->SetInput( input );
    warpImageFilter->Update();
    
  }

  try
    {
    warpImageFilter->Update();
    output = warpImageFilter->GetOutput();
    output->DisconnectPipeline();
    }
  catch( ExceptionObject & excp )
    {
    throw excp;
    }
  catch( ... )
    {
    ExceptionObject e( __FILE__, __LINE__, 
                              "Caught unknown exception", ITK_LOCATION );
    throw e;
    }		
 
}

} // end namespace idp
} //end namespace itk

#endif

  
