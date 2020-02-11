/*=========================================================================

  itkScaledDerivativeImageToImageMetricAdaptor.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkScaledDerivativeImageToImageMetricAdaptor_h
#define __itkScaledDerivativeImageToImageMetricAdaptor_h

#include "itkImageToImageMetric.h"


namespace itk
{
  
/** \class ScaledDerivativeImageToImageMetricAdaptor
 * \brief This class is an adaptor that rescaled the 
 * derivative of the underlying metric by user defined scaling factors.
 * 
 * \ingroup Numerics Optimizers
 */
template <typename TFixedImage, typename TMovingImage, typename TMetric> 
class ScaledDerivativeImageToImageMetricAdaptor: 
    public ImageToImageMetric<TFixedImage,TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef ScaledDerivativeImageToImageMetricAdaptor     Self;
  typedef ImageToImageMetric<TFixedImage,TMovingImage>   Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ScaledDerivativeImageToImageMetricAdaptor, ImageToImageMetric );

  /** Underlying metric type. */
  typedef TMetric                                   MetricType;
  typedef typename MetricType::Pointer              MetricPointer;

  /** Scales typedef */
  typedef Array<double>                             ScalesType;

  /** Inherit typedef from superclass. */
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::FixedImageType       FixedImageType;
  typedef typename Superclass::MovingImageType      MovingImageType;
  typedef typename Superclass::FixedImageMaskType   FixedImageMaskType;
  typedef typename Superclass::MovingImageMaskType  MovingImageMaskType;
  typedef typename Superclass::FixedImageRegionType FixedImageRegionType;
  typedef typename Superclass::TransformType        TransformType;
  typedef typename Superclass::InterpolatorType     InterpolatorType;
  typedef typename Superclass::MeasureType          MeasureType;
  typedef typename Superclass::DerivativeType       DerivativeType;

  /** Set the underlying metric */
  itkSetObjectMacro( Metric, MetricType );
    
  /** Get the underlying metric */
  itkGetObjectMacro( Metric, MetricType );
    
  /**  Delegate computation of the value to the underlying Metric. */
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const;
    
  /**  Delegate computation of the gradient to the costFunction.  */
  virtual void GetValueAndDerivative( const ParametersType & parameters,
                                      MeasureType & value,
                                      DerivativeType & derivative ) const;

  /**  Delegate computation of the value to the costFunction.  */
  virtual MeasureType GetValue( const ParametersType & parameters ) const;
    
  /** Set/Get current parameters scaling. */
  virtual void SetScales(const ScalesType & scales);
  itkGetMacro( Scales, ScalesType );

  /** Get the number of pixels considered in the computation. */
  virtual const unsigned long & GetNumberOfPixelsCounted() const;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject );

protected:
  ScaledDerivativeImageToImageMetricAdaptor();
  virtual ~ScaledDerivativeImageToImageMetricAdaptor() {};

private:
  ScaledDerivativeImageToImageMetricAdaptor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MetricPointer           m_Metric;
  bool                    m_ScalesInitialized;
  ScalesType              m_Scales;

};  // end of Class CostFunction

    
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScaledDerivativeImageToImageMetricAdaptor.txx"
#endif

#endif



