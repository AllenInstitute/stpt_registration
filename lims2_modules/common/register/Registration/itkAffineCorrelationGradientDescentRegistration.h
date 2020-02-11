/*=========================================================================

  itkAffineCorrelationGradientDescentRegistration.h
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkAffineCorrelationGradientDescentRegistration_h
#define __itkAffineCorrelationGradientDescentRegistration_h

#include "itkIterativeImageRegistrationMethod.h"
#include "itkCenteredAffineTransform.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkScaledDerivativeImageToImageMetricAdaptor.h"
#include "itkGradientDescentOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegistrationImageCoordinator.h"


namespace itk
{

/** \class AffineCorrelationGradientDescentRegistration
 * \brief Affine registration using correlation and gradient descent optimization.
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
class ITK_EXPORT AffineCorrelationGradientDescentRegistration : 
  public IterativeImageRegistrationMethod<TFixedImage, TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef AffineCorrelationGradientDescentRegistration  Self;
  typedef IterativeImageRegistrationMethod<TFixedImage,TMovingImage>  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(AffineCorrelationGradientDescentRegistration, IterativeImageRegistrationMethod);

  /**  Type of the Fixed image. */
  typedef typename Superclass::FixedImageType      FixedImageType;

  /**  Type of the Moving image. */
  typedef typename Superclass::MovingImageType     MovingImageType;

  /** Constants for the image dimensions */
  itkStaticConstMacro(FixedImageDimension, unsigned int,
                      TFixedImage::ImageDimension);
  itkStaticConstMacro(MovingImageDimension, unsigned int,
                      TMovingImage::ImageDimension);

  /**  Type of the metric. */
  typedef typename Superclass::MetricType          MetricType;
  typedef typename MetricType::CoordinateRepresentationType 
                                             CoordinateRepresentationType;
  typedef NormalizedCorrelationImageToImageMetric< 
                            FixedImageType,
                            MovingImageType>  CorrelationMetricType;

  typedef ScaledDerivativeImageToImageMetricAdaptor<
                            FixedImageType,
                            MovingImageType,
                            CorrelationMetricType > MetricAdaptorType;

  /**  Type of the Transform . */
  typedef CenteredAffineTransform<CoordinateRepresentationType, FixedImageDimension> AffineTransformType;

  /**  Type of the Interpolator. */
  typedef LinearInterpolateImageFunction<
                  MovingImageType,
                  CoordinateRepresentationType> LinearInterpolatorType;
  typedef typename LinearInterpolatorType::Pointer LinearInterpolatorPointer;


  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef  typename Superclass::ParametersType    ParametersType;

  /** Type of the image coordinator. */
  typedef RegistrationImageCoordinator<
                                   FixedImageType,
                                   MovingImageType >  CoordinatorType;
  typedef typename CoordinatorType::Pointer           CoordinatorPointer;


  /** Get the correlation metric. */
  virtual CorrelationMetricType * GetCorrelationMetric();
  virtual MetricAdaptorType *  GetMetricAdaptor();

  /** Get the affine transform. */
  virtual AffineTransformType * GetAffineTransform();

  /** Get the linear interpolator. */
  virtual LinearInterpolatorType * GetLinearInterpolator();

  /** Get the gradient descent optimizer. */
  virtual GradientDescentOptimizer * GetGradientDescentOptimizer();

  /** Set/Get the cooridinator. */
  itkSetObjectMacro( Coordinator, CoordinatorType );
  itkGetObjectMacro( Coordinator, CoordinatorType );

  /** Optimizer iteration handler. */
  void OptimizerIterationHandler();

  /** Set/Get the matrix parameter scale. */
  itkSetMacro( MatrixScale, double );
  itkGetMacro( MatrixScale, double );

  /** Set/Get the angle parameter scale rate. */
  itkSetMacro( MatrixScaleRate, double );
  itkGetMacro( MatrixScaleRate, double );

  /** Set/Get the matrix parameter scale. */
  itkSetMacro( CenterScale, double );
  itkGetMacro( CenterScale, double );

  /** Set/Get the step size. */
  itkSetMacro( StepSize, double );
  itkGetMacro( StepSize, double );

  /** Set/Get the step size rate. */
  itkSetMacro( StepSizeRate, double );
  itkGetMacro( StepSizeRate, double );

  /** Set/Get the number of iterations. */
  itkSetMacro( NumberOfIterations, unsigned long );
  itkGetMacro( NumberOfIterations, unsigned long );

  /** Set/Get the iteration rate. */
  itkSetMacro( IterationRate, double );
  itkGetMacro( IterationRate, double );

  /** Set/Get the verbosity flag. */
  itkSetMacro( Verbose, bool );
  itkGetMacro( Verbose, bool );
  itkBooleanMacro( Verbose );

protected:
  AffineCorrelationGradientDescentRegistration();
  virtual ~AffineCorrelationGradientDescentRegistration();
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
  AffineCorrelationGradientDescentRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  CoordinatorPointer         m_Coordinator;  

  double                     m_MatrixScale;
  double                     m_MatrixScaleRate;
  double                     m_CenterScale;

  double                     m_StepSize;
  double                     m_StepSizeRate;
  unsigned long              m_NumberOfIterations;
  double                     m_IterationRate;

  unsigned long              m_ObserverTag;
  bool                       m_Verbose;
    
};


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAffineCorrelationGradientDescentRegistration.txx"
#endif

#endif



