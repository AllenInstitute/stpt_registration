/*=========================================================================

  itkChangePixelSpacingImageFilter.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkChangePixelSpacingImageFilter_txx
#define _itkChangePixelSpacingImageFilter_txx

#include "itkChangePixelSpacingImageFilter.h"
#include "itkProgressAccumulator.h"

#include "itkCastImageFilter.h"

#include "itkGaussianOperator.h"
#include "itkNeighborhoodOperatorImageFilter.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"


namespace itk
{

/**
 * Default constructor
 */
template <class TInputImage, class TOutputImage>
ChangePixelSpacingImageFilter<TInputImage, TOutputImage>
::ChangePixelSpacingImageFilter()
{
  // Default output spacing
  m_Spacing.Fill( 1.0 );

  // Default Gaussian kernel parameters
  m_SmoothingVariance.Fill( 1.0 );
  m_MaximumKernelError = 0.1;
  m_MaximumKernelWidth = 5;
  
  // default pixel
  m_DefaultPixelValue = 0;

  // Default interpolation method
  typedef BSplineInterpolateImageFunction<InternalImageType,double> BSplineInterpolatorType;
  typename BSplineInterpolatorType::Pointer interp = BSplineInterpolatorType::New();

  interp->SetSplineOrder( 3 );

  m_Interpolator = interp;

}


/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
ChangePixelSpacingImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Spacing: " << m_Spacing << std::endl;
  os << indent << "SmoothingVariance: " << m_SmoothingVariance << std::endl;
  os << indent << "MaximumKernelError: " << m_MaximumKernelError << std::endl;
  os << indent << "MaximumKernelWidth: " << m_MaximumKernelWidth << std::endl;
  os << indent << "DefaultPixelValue: " << m_DefaultPixelValue << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
}


/**
 * ThreadedGenerateData
 */
template <class TInputImage, class TOutputImage>
void
ChangePixelSpacingImageFilter<TInputImage, TOutputImage>
::GenerateData()
{

  // Get the input and output pointers
  OutputImagePointer outputPtr = this->GetOutput();
  InputImagePointer  inputPtr  = this->GetInput();

  if( !inputPtr || !outputPtr )
    {
    return;
    }

  // Cast input to internal image type
  typedef CastImageFilter<InputImageType,InternalImageType> InputCasterType;
  typename InputCasterType::Pointer inputCaster = InputCasterType::New();

  inputCaster->SetInput( inputPtr );
  inputCaster->UpdateLargestPossibleRegion();

  // Smooth the image
  typedef GaussianOperator<float,ImageDimension> OperatorType;
  typedef NeighborhoodOperatorImageFilter<InternalImageType,InternalImageType> SmootherType;

  typename InternalImageType::Pointer smoothedImage = inputCaster->GetOutput();

  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    if ( m_SmoothingVariance[j] == 0.0 )
      {
      continue;
      }
    std::cout << "Smoothing direction " << j << std::endl;
    typename SmootherType::Pointer filter = SmootherType::New();
    OperatorType * oper = new OperatorType;
    oper->SetDirection( j );
    oper->SetVariance( m_SmoothingVariance[j] );
    oper->SetMaximumKernelWidth( m_MaximumKernelWidth );
    oper->SetMaximumError( m_MaximumKernelError );
    oper->CreateDirectional();
    filter->SetOperator( *oper );
    filter->SetInput( smoothedImage );
    filter->UpdateLargestPossibleRegion();
    smoothedImage = filter->GetOutput();
    smoothedImage->DisconnectPipeline();
    delete oper;
    if ( j == 0 )
      {
      inputCaster = NULL;
      }
    }


  // Resample the image
  typedef IdentityTransform<double,ImageDimension> TransformType;
  typename TransformType::Pointer transform = TransformType::New();

  typedef ResampleImageFilter<InternalImageType,OutputImageType> ResamplerType;
  typename ResamplerType::Pointer resampler = ResamplerType::New();

  resampler->SetInput( smoothedImage );
  resampler->SetTransform( transform );
  resampler->SetInterpolator( m_Interpolator );
  resampler->SetOutputOrigin( outputPtr->GetOrigin() );
  resampler->SetOutputSpacing( outputPtr->GetSpacing() );
  resampler->SetDefaultPixelValue( m_DefaultPixelValue );
  resampler->SetSize( outputPtr->GetLargestPossibleRegion().GetSize() );
  resampler->UpdateLargestPossibleRegion();
 
  this->GraftOutput( resampler->GetOutput() );

}


/**
 * EnlargeOutputRequestedRegion
 */
template <class TInputImage, class TOutputImage>
void
ChangePixelSpacingImageFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion( DataObject * data )
{

  // Call the superclass' implementation of this method
  Superclass::EnlargeOutputRequestedRegion( data );
  if ( data )
    {
    data->SetRequestedRegionToLargestPossibleRegion();
    }
}

/**
 * GenerateOutputInformation
 */
template <class TInputImage, class TOutputImage>
void
ChangePixelSpacingImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{

  // Call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // Get the input and output pointers
  OutputImagePointer outputPtr = this->GetOutput();
  InputImagePointer  inputPtr  = this->GetInput();

  if( !inputPtr || !outputPtr )
    {
    return;
    }

  // Set the output spacing
  outputPtr->SetSpacing( m_Spacing.GetDataPointer() );

  // Compute the start and end points of the input
  typedef typename InputImageType::IndexType IndexType;
  typedef Point<double,ImageDimension> PointType;

  Offset<ImageDimension>       offset;
  offset.Fill( -1 );

  IndexType inputStartIndex = inputPtr->GetLargestPossibleRegion().GetIndex();
  IndexType inputEndIndex   = inputStartIndex +
    inputPtr->GetLargestPossibleRegion().GetSize() + offset;

  PointType inputStartPoint;
  PointType inputEndPoint;

  inputPtr->TransformIndexToPhysicalPoint( inputStartIndex, inputStartPoint );
  inputPtr->TransformIndexToPhysicalPoint( inputEndIndex,   inputEndPoint   );

  // Set the output origin
  outputPtr->SetOrigin( inputStartPoint.GetDataPointer() );

  typedef ContinuousIndex<double,ImageDimension> ContinuousIndexType;
  ContinuousIndexType outputEndContIndex;
  
  outputPtr->TransformPhysicalPointToContinuousIndex( inputEndPoint,
    outputEndContIndex );

  typedef typename OutputImageType::SizeType SizeType;
  SizeType size;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    size[i] = static_cast<long>( outputEndContIndex[i] + 1 + 0.5 );
    }

  typedef typename OutputImageType::RegionType RegionType;
  RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( size );

  outputPtr->SetLargestPossibleRegion( largestPossibleRegion );

}


/**
 * GenerateInputRequestedRegion
 */
template <class TInputImage, class TOutputImage>
void
ChangePixelSpacingImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{

  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  InputImageType * input = 
        const_cast< InputImageType * >( this->GetInput() );
  if( input )
    {
    input->SetRequestedRegionToLargestPossibleRegion();
    }

}


} // end namespace itk

#endif
