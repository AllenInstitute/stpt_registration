/*=========================================================================

  idpMutualInfoVersorRigidVolumeRegistration.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpMutualInfoVersorRigidVolumeRegistration_h
#define __idpMutualInfoVersorRigidVolumeRegistration_h

#include "itkIterativeImageRegistrationMethod.h"
#include "itkVersorRigid3DTransform.h"
#include "itkMattesMutualInformationImageToImageMetric.h"

#include "itkSimpleRegistrationImageCoordinator.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegistrationImageCoordinator.h"

namespace itk
{
    namespace idp
    {

/** \class MutualInfoVersorRigidVolumeRegistration
 * \brief Versor rigid 3D registration using mutual information and gradient descent optimization.
 *
 * The registration process is initiated by method StartRegistration().
 * The user must set the parameters of each component before calling
 * this method.
 *
 * The number of level to process can be set via
 * SetNumberOfLevels(). At each level, the user specified 
 * registration components are used to register the current input
 * images by computing the transform parameters that will map one image onto 
 * the other image.
 *
 * Before each level an IterationEvent is invoked providing an
 * opportunity for a user interface to change any of the components,
 * change component parameters, or stop the registration.
 *
 * This class is templated over the fixed image type and the moving image
 * type.
 *
 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */
template <typename TFixedImage, typename TMovingImage>
class ITK_EXPORT MutualInfoVersorRigidVolumeRegistration : 
  public IterativeImageRegistrationMethod<TFixedImage, TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef MutualInfoVersorRigidVolumeRegistration  Self;
  typedef IterativeImageRegistrationMethod<TFixedImage,TMovingImage>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MutualInfoVersorRigidVolumeRegistration, IterativeImageRegistrationMethod);

  /**  Type of the Fixed image. */
  typedef typename Superclass::FixedImageType      FixedImageType;

  /**  Type of the Moving image. */
  typedef typename Superclass::MovingImageType     MovingImageType;

  /** Constants for the image dimensions */
  itkStaticConstMacro(FixedImageDimension, unsigned int, TFixedImage::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int, TMovingImage::ImageDimension);

  /**  Type of the metric. */
  typedef typename Superclass::MetricType                      MetricType;
  typedef typename MetricType::CoordinateRepresentationType    CoordinateRepresentationType;
  typedef MattesMutualInformationImageToImageMetric<FixedImageType,MovingImageType>  MutualInfoMetricType;

  /**  Type of the Transform . */
  typedef VersorRigid3DTransform<CoordinateRepresentationType> VersorTransformType;

  /**  Type of the Interpolator. */
  typedef LinearInterpolateImageFunction<MovingImageType,CoordinateRepresentationType> LinearInterpolatorType;
  typedef typename LinearInterpolatorType::Pointer LinearInterpolatorPointer;
  
  /**  Type of the optimizer. */
  typedef VersorRigid3DTransformOptimizer OptimizerType;
  
  /**  Type of the image-coordinator. */
  typedef SimpleRegistrationImageCoordinator< FixedImageType,MovingImageType >    CoordinatorType;

  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef  typename Superclass::ParametersType    ParametersType;

  /** Get the metric. */
  virtual MutualInfoMetricType * GetMutualInfoMetric();

  /** Get the affine transform. */
  virtual VersorTransformType * GetVersorTransform();

  /** Get the linear interpolator. */
  virtual LinearInterpolatorType * GetLinearInterpolator();

  /** Get the gradient descent optimizer. */
  virtual VersorRigid3DTransformOptimizer * GetGradientDescentOptimizer();

  /** Set/Get the cooridinator. */
  itkSetObjectMacro( Coordinator, CoordinatorType );
  itkGetObjectMacro( Coordinator, CoordinatorType );
  
  /** Set/Get the input Fixed image. */
  itkSetConstObjectMacro( FixedInputImage, FixedImageType );
  itkGetConstObjectMacro( FixedInputImage, FixedImageType ); 

  /** Set/Get the input Moving image. */
  itkSetConstObjectMacro( MovingInputImage, MovingImageType );
  itkGetConstObjectMacro( MovingInputImage, MovingImageType );

  /** Optimizer iteration handler. */
  void OptimizerIterationHandler();
  
  /** Load registration parameters from XML file. */
  void LoadParametersFromXML( const char * fileName );

  //itkSetMacro( NumberOfLevels, int );
  //itkGetConstMacro( NumberOfLevels, int );
  itkSetMacro( NumberOfIterations, int );
  itkGetConstMacro( NumberOfIterations, int );
  itkSetMacro( MaximumStepLength, double );
  itkGetConstMacro( MaximumStepLength, double );
  itkSetMacro( MinimumStepLength, double );
  itkGetConstMacro( MinimumStepLength, double );
  itkSetMacro( MaximumStepLengthRate, double );
  itkGetConstMacro( MaximumStepLengthRate, double );
  itkSetMacro( MinimumStepLengthRate, double );
  itkGetConstMacro( MinimumStepLengthRate, double );
  itkSetMacro( AngleScale, double );
  itkGetConstMacro( AngleScale, double );
  itkSetMacro( TranslationScale, double );
  itkGetConstMacro( TranslationScale, double );
  itkSetMacro( NumberOfHistogramBins, int );
  itkGetConstMacro( NumberOfHistogramBins, int );
  itkSetMacro( NumberOfSpatialSamples, int );
  itkGetConstMacro( NumberOfSpatialSamples, int );
  itkSetMacro( RelaxationFactor, double );
  itkGetConstMacro( RelaxationFactor, double );
  itkSetMacro( RelaxationFactorRate, double );
  itkGetConstMacro( RelaxationFactorRate, double );

  /** Set/Get the verbosity flag. */
  itkSetMacro( Verbose, bool );
  itkGetMacro( Verbose, bool );
  itkBooleanMacro( Verbose );
  
  /** Get output affine transformation */
  itkGetConstObjectMacro( OutputTransform, VersorTransformType );

protected:
  MutualInfoVersorRigidVolumeRegistration();
  virtual ~MutualInfoVersorRigidVolumeRegistration();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This method is called before starting registration at each level. */
  virtual void BeforeIteration();

  /** This method is called after registration at each level. */
  virtual void AfterIteration();

  /** This method is called before starting registration. */
  virtual void BeforeRegistration();

  /** This method is called after registration has finished or is stopped. */
  virtual void AfterRegistration();


private:
  MutualInfoVersorRigidVolumeRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double                    m_RelaxationFactorRate;
  int                       m_NumberOfIterations;
  double                    m_MaximumStepLength;
  double                    m_MinimumStepLength;
  double                    m_MaximumStepLengthRate;
  double                    m_MinimumStepLengthRate;
  double                    m_AngleScale;
  double                    m_TranslationScale;
  int                       m_NumberOfHistogramBins;
  int                       m_NumberOfSpatialSamples;
  double                    m_RelaxationFactor;

  unsigned long             m_ObserverTag;
  bool                      m_Verbose;
    
  typename CoordinatorType::Pointer            m_Coordinator;
  typename FixedImageType::ConstPointer        m_FixedInputImage; //the input image at fullsize and used as input to m_Coordinator.
  typename MovingImageType::ConstPointer       m_MovingInputImage;//the downsampled output from m_Coordinator is the actrual input to superclass
  
  typename VersorTransformType::Pointer        m_OutputTransform;
};

} // end namespace idp

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "idpMutualInfoVersorRigidVolumeRegistration.txx"
#endif

#endif
