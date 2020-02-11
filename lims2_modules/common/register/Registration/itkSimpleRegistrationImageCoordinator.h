/*=========================================================================

  itkSimpleRegistrationImageCoordinator.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkSimpleRegistrationImageCoordinator_h
#define _itkSimpleRegistrationImageCoordinator_h

#include "itkRegistrationImageCoordinator.h"
#include "itkFixedArray.h"

namespace itk
{

/** \class SimpleRegistrationImageCoordinator
 * \brief Helper class to coordinate images for IterativeRegistrationMethod
 */
template < typename TFixedImage, typename TMovingImage >
class ITK_EXPORT SimpleRegistrationImageCoordinator : 
  public RegistrationImageCoordinator<TFixedImage,TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef SimpleRegistrationImageCoordinator Self;
  typedef RegistrationImageCoordinator<TFixedImage,TMovingImage> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(SimpleRegistrationImageCoordinator, RegistrationImageCoordinator);

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Constants for the image dimensions */
  itkStaticConstMacro(FixedImageDimension, unsigned int,
                      TFixedImage::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int,
                      TMovingImage::ImageDimension);

  /**  Type of the Fixed image. */
  typedef typename Superclass::FixedImageType     FixedImageType;
  typedef typename FixedImageType::Pointer        FixedImagePointer;
  typedef typename FixedImageType::RegionType     RegionType;

  /** Type of the Fixed image mask. */
  typedef ImageMaskSpatialObject<FixedImageDimension>   FixedImageMaskType;
  typedef typename FixedImageMaskType::Pointer          FixedImageMaskPointer;
  typedef typename FixedImageMaskType::ImageType        FixedMaskImageType;
  typedef typename FixedMaskImageType::Pointer          FixedMaskImagePointer;

  /**  Type of the Moving image. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename MovingImageType::Pointer       MovingImagePointer;

  /** Type of the Moving image mask. */
  typedef ImageMaskSpatialObject<MovingImageDimension>    MovingImageMaskType;
  typedef typename MovingImageMaskType::Pointer           MovingImageMaskPointer;
  typedef typename MovingImageMaskType::ImageType         MovingMaskImageType;
  typedef typename MovingMaskImageType::Pointer           MovingMaskImagePointer;

  /** Type of the registration class. */
  typedef  typename Superclass::TransformType      TransformType;
  typedef  typename Superclass::InputPointType     InputPointType;
  typedef  typename Superclass::OutputVectorType   OutputVectorType;

  /** Set/Get the fixed image. */
  itkSetConstObjectMacro( FixedImage, FixedImageType );
  itkGetConstObjectMacro( FixedImage, FixedImageType );

  /** Set/Get the fixed image mask. */
  itkSetConstObjectMacro( FixedMaskImage, FixedMaskImageType );
  itkGetConstObjectMacro( FixedMaskImage, FixedMaskImageType );

  /** Set/Get the moving image. */
  itkSetConstObjectMacro( MovingImage, MovingImageType );
  itkGetConstObjectMacro( MovingImage, MovingImageType );

  /** Set/Get the moving image mask. */
  itkSetConstObjectMacro( MovingMaskImage, MovingMaskImageType );
  itkGetConstObjectMacro( MovingMaskImage, MovingMaskImageType );

  /** Type of defining the shrink factors. */
  typedef FixedArray<unsigned int, itkGetStaticConstMacro( FixedImageDimension ) > FactorsType;

  /** Set/Get the fixed image starting shrink factors. */
  itkSetMacro( FixedImageStartingFactors, FactorsType );
  itkGetMacro( FixedImageStartingFactors, FactorsType );

  /** Set/Get the moving image starting shrink factors. */
  itkSetMacro( MovingImageStartingFactors, FactorsType );
  itkGetMacro( MovingImageStartingFactors, FactorsType );

  /** Set/Get verbose flag. */
  itkSetMacro( Verbose, bool );
  itkGetConstMacro( Verbose, bool );

   /** Initialize. This method must be called before using any of the
     Get image methods and get transformation center and translation 
     methods.
   */
  virtual void Initialize() throw (ExceptionObject);

  /** Get the fixed image for a particular level. */
  virtual FixedImageType * GetFixedImage( unsigned long level ) const;

  /** Get the mask spatial object associated with the fixed image. */
  virtual FixedImageMaskType * GetFixedImageMask() const;

  /** Get the moving image for a particular level. */
  virtual MovingImageType * GetMovingImage( unsigned long level ) const;

  /** Get the mask spatial object associated with the moving image. */
  virtual MovingImageMaskType * GetMovingImageMask() const;

protected:
  SimpleRegistrationImageCoordinator ();
  virtual ~SimpleRegistrationImageCoordinator () {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  SimpleRegistrationImageCoordinator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::vector<FixedImagePointer>                  m_FixedImages;
  typename FixedImageType::ConstPointer           m_FixedImage;
  FixedImageMaskPointer                           m_FixedImageMask;
  typename FixedMaskImageType::ConstPointer       m_FixedMaskImage;
  FactorsType                                     m_FixedImageStartingFactors;

  std::vector<MovingImagePointer>                 m_MovingImages;
  typename MovingImageType::ConstPointer          m_MovingImage;
  MovingImageMaskPointer                          m_MovingImageMask;
  typename MovingMaskImageType::ConstPointer      m_MovingMaskImage;
  FactorsType                                     m_MovingImageStartingFactors;

  bool                                            m_Verbose;
    
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSimpleRegistrationImageCoordinator.txx"
#endif

#endif
