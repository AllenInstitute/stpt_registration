/*=========================================================================

  itkRegistrationImageCoordinator.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkRegistrationImageCoordinator_h
#define _itkRegistrationImageCoordinator_h

#include "itkObject.h"
#include "itkImageToImageMetric.h"
#include "itkImageMaskSpatialObject.h"


namespace itk
{

/** \class RegistrationImageCoordinator
 * \brief Helper class to coordinate images for IterativeRegistrationMethod
 */
template <typename TFixedImage, typename TMovingImage>
class ITK_EXPORT RegistrationImageCoordinator : public Object 
{
public:

  /** Standard class typedefs. */
  typedef RegistrationImageCoordinator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(RegistrationImageCoordinator, Object);

  /**  Type of the Fixed image. */
  typedef          TFixedImage                      FixedImageType;
 
  /**  Type of the Moving image. */
  typedef          TMovingImage                     MovingImageType;

  /** Type transform. */
  typedef     ImageToImageMetric< FixedImageType,
                                  MovingImageType>    MetricType;
  typedef  typename MetricType::TransformType         TransformType;
  typedef  typename TransformType::InputPointType     InputPointType;
  typedef  typename TransformType::OutputVectorType   OutputVectorType;

  /** Set/Get the number of registration levels. */
  itkSetClampMacro( NumberOfLevels, unsigned long, 1,
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfLevels, unsigned long );

   /** Initialize. This method must be called before using any of the
     Get image methods and get transformation center and translation 
     methods.
   */
  virtual void Initialize() throw (ExceptionObject) {};

  /** Get the fixed image for a particular level. */
  virtual FixedImageType * GetFixedImage( unsigned long level ) const = 0;

  /** Get the moving image for a particular level. */
  virtual MovingImageType * GetMovingImage( unsigned long level ) const = 0;

  /** Constants for the image dimensions */
  itkStaticConstMacro(FixedImageDimension, unsigned int,
                      TFixedImage::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int,
                      TMovingImage::ImageDimension);

  /** Type of the Fixed image mask. */
  typedef ImageMaskSpatialObject<FixedImageDimension>   FixedImageMaskType;
  typedef typename FixedImageMaskType::Pointer          FixedImageMaskPointer;

  /** Type of the Moving image mask. */
  typedef ImageMaskSpatialObject<MovingImageDimension>    MovingImageMaskType;
  typedef typename MovingImageMaskType::Pointer           MovingImageMaskPointer;


  /** Get the mask spatial object associated with the fixed image. */
  virtual FixedImageMaskType * GetFixedImageMask() const
    { return NULL; }

  /** Get the mask spatial object associated with the moving image. */
  virtual MovingImageMaskType * GetMovingImageMask() const
    { return NULL; }

  /** Get the transform center. Typically this is the center of the fixed image. */
  itkGetMacro( TransformCenter, InputPointType );

  /** Get the initial translation. Typically this is the distance between
   * the fixed image center and the moving image center. */
  itkGetMacro( InitialTranslation, OutputVectorType );

protected:
  RegistrationImageCoordinator ()
    { 
    m_NumberOfLevels = 0;
    m_TransformCenter.Fill( 0.0 );
    m_InitialTranslation.Fill( 0.0 );
    };

  virtual ~RegistrationImageCoordinator () {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    this->Superclass::PrintSelf( os, indent );
    os << indent << "NumberOfLevels: " << m_NumberOfLevels << std::endl;
    os << indent << "TransformCenter: " << m_TransformCenter << std::endl;
    os << indent << "InitialTranslation: " << m_InitialTranslation << std::endl;
    }

  InputPointType               m_TransformCenter;
  OutputVectorType             m_InitialTranslation;

private:
  RegistrationImageCoordinator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned long                m_NumberOfLevels;
  
};

} // namespace itk

#endif
