/*=========================================================================

  itkChangePixelSpacingImageFilter.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkChangePixelSpacingImageFilter_h
#define __itkChangePixelSpacingImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkInterpolateImageFunction.h"

namespace itk
{

/** \class ChangePixelSpacingImageFilter
 * \brief Reslice an image to user specified pixel spacing.
 *
 * ChangePixelSpacingImageFilter reslices an image to a user specified
 * pixel spacing. The input image is first separably smoothed using a Gaussian
 * kernel of user specified variance. The smoothed image is then resampled to
 * to produce an image of the required pixel spacing. By default a B-spline interpolator
 * is used for the resample. The interpolation can be replaced using the 
 * SetInterpolator method.
 *
 */
template < 
  class TInputImage,
  class TOutputImage >
class ITK_EXPORT ChangePixelSpacingImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ChangePixelSpacingImageFilter         Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(ChangePixelSpacingImageFilter, ImageToImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Inherit some types from superclass. */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename OutputImageType::PixelType  OutputPixelType;
  typedef typename InputImageType::ConstPointer  InputImagePointer;
  typedef typename OutputImageType::Pointer  OutputImagePointer;

  /** Typedef for input parameter arrays. */
  typedef FixedArray<double,ImageDimension> ArrayType;

  /** Typedef for the internal image type. */
  typedef Image<float,ImageDimension> InternalImageType;

  /** Typedef of interpolator. */
  typedef InterpolateImageFunction<InternalImageType,double> InterpolatorType;

  /** Set/Get the interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );
  itkGetObjectMacro( Interpolator, InterpolatorType );


  /** Set/Get the output image spacing. */
  itkSetMacro( Spacing, ArrayType );
  itkGetConstMacro( Spacing, ArrayType );

  /** Set/Get the variance of the Gaussian smoothing kernel for each dimension. */
  itkSetMacro( SmoothingVariance, ArrayType );
  itkGetConstMacro( SmoothingVariance, ArrayType );

  /** Set/Get the maximum difference between the continuous Gaussian kernel
   * and the discrete approximation. */
  itkSetMacro( MaximumKernelError, double );
  itkGetConstMacro( MaximumKernelError, double );

  /** Set/Get the maximum kernel size allowed */
  itkSetMacro( MaximumKernelWidth, unsigned int );
  itkGetConstMacro( MaximumKernelWidth, unsigned int );
  
  /** Set/Get the default outtput pixel. */
  itkSetMacro( DefaultPixelValue, OutputPixelType );
  itkGetConstMacro( DefaultPixelValue, OutputPixelType );

  /** This filter produces an output which is of a different size to the
   * input. As such, it needs to provide an implementation for
   * GenerateOutputInformation.
   * The original documentation of this method is below.
   * \sa ProcessObject::GenerateOutputInformation() */
  virtual void GenerateOutputInformation();

  /** This filter requires all of the output image
   * to be produced at once. As such, it needs to provide an implemenation for
   * EnlargeOutputRequestedRegion.
   * The original documentation of this method is below.
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  virtual void EnlargeOutputRequestedRegion(DataObject *);

  /** This filter requires all of the input image to be present
   * to compute the output. As such, it needs to provide and implementation for
   * GenerateInputRequestedRegion.
   * The original documentation of this method is below.
   * \sa GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion();


protected:
  ChangePixelSpacingImageFilter();
  ~ChangePixelSpacingImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This method causes the filter to generate its output. */
  virtual void GenerateData();
  
private:
  ChangePixelSpacingImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Pixel spacing for the output images. */
  ArrayType       m_Spacing;
 
  /** Variance of gaussian smoothing kernels for each dimension. */
  ArrayType       m_SmoothingVariance;

  double          m_MaximumKernelError;
  unsigned int    m_MaximumKernelWidth;

  typename InterpolatorType::Pointer   m_Interpolator;
  
  OutputPixelType m_DefaultPixelValue;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkChangePixelSpacingImageFilter.txx"
#endif

#endif
