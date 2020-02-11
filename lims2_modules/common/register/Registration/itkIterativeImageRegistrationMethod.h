/*=========================================================================

  itkIterativeImageRegistrationMethod.h
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkIterativeImageRegistrationMethod_h
#define __itkIterativeImageRegistrationMethod_h

#include "itkProcessObject.h"
#include "itkImageToImageMetric.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class IterativeImageRegistrationMethod
 * \brief Base class for iterative image registration methods
 *
 * This class provides a generic interface for iterative 
 * registration using components of the registration framework.
 * See documentation for ImageRegistrationMethod for a description
 * of the registration framework components.
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
class ITK_EXPORT IterativeImageRegistrationMethod : public ProcessObject 
{
public:
  /** Standard class typedefs. */
  typedef IterativeImageRegistrationMethod  Self;
  typedef ProcessObject  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(IterativeImageRegistrationMethod, ProcessObject);

  /**  Type of the Fixed image. */
  typedef          TFixedImage                     FixedImageType;
  typedef typename FixedImageType::Pointer         FixedImagePointer;
  typedef typename FixedImageType::ConstPointer    FixedImageConstPointer;
  typedef typename FixedImageType::RegionType      FixedImageRegionType;

  /**  Type of the Moving image. */
  typedef          TMovingImage                    MovingImageType;
  typedef typename MovingImageType::ConstPointer   MovingImageConstPointer;

  /**  Type of the metric. */
  typedef ImageToImageMetric< FixedImageType,
                              MovingImageType >    MetricType;
  typedef typename MetricType::Pointer             MetricPointer;

  /**  Type of the Transform . */
  typedef typename MetricType::TransformType       TransformType;
  typedef typename TransformType::Pointer          TransformPointer;

  /**  Type of the Interpolator. */
  typedef typename MetricType::InterpolatorType    InterpolatorType;
  typedef typename InterpolatorType::Pointer       InterpolatorPointer;

  /**  Type of the optimizer. */
  typedef SingleValuedNonLinearOptimizer           OptimizerType;

  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef  typename MetricType::TransformParametersType    ParametersType;

  /** Method that initiates the registration. */
  void StartRegistration();

  /** Method to stop the registration. */
  void StopRegistration();

  /** Set/Get the Fixed image. */
  itkSetConstObjectMacro( FixedImage, FixedImageType );
  itkGetConstObjectMacro( FixedImage, FixedImageType ); 

  /** Set/Get the Moving image. */
  itkSetConstObjectMacro( MovingImage, MovingImageType );
  itkGetConstObjectMacro( MovingImage, MovingImageType );

  /** Set/Get the Optimizer. */
  itkSetObjectMacro( Optimizer,  OptimizerType );
  itkGetObjectMacro( Optimizer,  OptimizerType );

  /** Set/Get the Metric. */
  itkSetObjectMacro( Metric, MetricType );
  itkGetObjectMacro( Metric, MetricType );

  /** Set/Get the Metric. */
  itkSetMacro( FixedImageRegion, FixedImageRegionType );
  itkGetConstReferenceMacro( FixedImageRegion, FixedImageRegionType );

  /** Set/Get the Transfrom. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetObjectMacro( Transform, TransformType );

  /** Set/Get the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );
  itkGetObjectMacro( Interpolator, InterpolatorType );

  /** Set/Get the number of registration levels. */
  itkSetClampMacro( NumberOfLevels, unsigned long, 1,
                    NumericTraits<unsigned long>::max() );
  itkGetMacro( NumberOfLevels, unsigned long );

  /** Get the current resolution level being processed. */
  itkGetMacro( CurrentLevel, unsigned long );

  /** Set/Get the initial transformation parameters. */
  itkSetMacro( InitialTransformParameters, ParametersType );
  itkGetConstReferenceMacro( InitialTransformParameters, ParametersType );

  /** Set/Get the initial transformation parameters of the next resolution
   level to be processed. The default is the last set of parameters of
   the last resolution level. */
  itkSetMacro( InitialTransformParametersOfNextLevel, ParametersType );
  itkGetConstReferenceMacro( InitialTransformParametersOfNextLevel, ParametersType );

  /** Get the last transformation parameters visited by 
   * the optimizer. */
  itkGetConstReferenceMacro( LastTransformParameters, ParametersType );  

protected:
  IterativeImageRegistrationMethod();
  virtual ~IterativeImageRegistrationMethod() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize by setting the interconnects between the components.
      This method is executed at every level of the pyramid with the
      values corresponding to this resolution
   */
  void Initialize() throw (ExceptionObject);

  /** Set the current level to be processed */  
  itkSetMacro( CurrentLevel, unsigned long );

  /** This method is called before starting registration at each level. */
  virtual void BeforeIteration() {};

  /** This method is called after registration at each level. */
  virtual void AfterIteration() {};

  /** This method is called before starting registration. */
  virtual void BeforeRegistration() {};

  /** This method is called after registration has finished or is stopped. */
  virtual void AfterRegistration() {};


private:
  IterativeImageRegistrationMethod(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  MetricPointer                    m_Metric;
  OptimizerType::Pointer           m_Optimizer;

  MovingImageConstPointer          m_MovingImage;
  FixedImageConstPointer           m_FixedImage;

  TransformPointer                 m_Transform;
  InterpolatorPointer              m_Interpolator;

  ParametersType                   m_InitialTransformParameters;
  ParametersType                   m_InitialTransformParametersOfNextLevel;
  ParametersType                   m_LastTransformParameters;

  FixedImageRegionType             m_FixedImageRegion;

  unsigned long                    m_NumberOfLevels;
  unsigned long                    m_CurrentLevel;

  bool                             m_Stop;
  
};


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIterativeImageRegistrationMethod.txx"
#endif

#endif



