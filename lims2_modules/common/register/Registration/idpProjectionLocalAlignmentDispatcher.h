/*=========================================================================

  idpProjectionAlignmentDispatcher.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpProjectionLocalAlignmentDispatcher_h
#define __idpProjectionLocalAlignmentDispatcher_h
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkRecursiveMultiResolutionPyramidImageFilter.h>
#include "itkLinearInterpolateImageFunction.h"
#include "itkImage.h"
#include "idpRegistrationDispatcher.h"

namespace itk
{
namespace idp
{

/** \class ProjectionLocalAlignmentDispatcher
 *
 */
template <class TPixel>
class ProjectionLocalAlignmentDispatcher : public RegistrationDispatcher<TPixel>
{
public:

  /** Standard typedefs. */
  typedef ProjectionLocalAlignmentDispatcher  Self;
  typedef RegistrationDispatcher<TPixel>  Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(ProjectionLocalAlignmentDispatcher, RegistrationDispatcher);
  
  /**  Inherited types. */
  typedef TPixel PixelType;
  typedef typename Superclass::ImageSeriesPointer   ImageSeriesPointer;
  typedef typename Superclass::ImageSeriesVector    ImageSeriesVector;
  typedef typename Superclass::SpecimenPointer      SpecimenPointer;
  
  typedef typename Superclass::Rigid2DTransformType     Rigid2DTransformType;
  typedef typename Superclass::Affine3DTransformType    Affine3DTransformType;
  typedef typename Superclass::Affine2DTransformType    Affine2DTransformType;
  typedef typename Superclass::ParametersType           ParametersType;
  typedef typename Superclass::Affine3DParametersType   Affine3DParametersType;
  
  typedef typename Superclass::ImageSeriesUtilitiesType     ImageSeriesUtilitiesType;
  typedef typename Superclass::VolumeType                   VolumeType;
  typedef typename Superclass::MaskVolumeType               MaskVolumeType;
  typedef typename Superclass::ImageType                    ImageType;
  typedef typename Superclass::MaskImageType                MaskImageType;

  typedef typename Superclass::PointType              PointType;
  typedef typename Superclass::RegionOfInterestType   RegionOfInterestType;
  
  typedef typename Superclass::ArrayType              ArrayType;  
  
  /**  new types. */
  typedef unsigned short IOPixelType;
  typedef itk::Image<float, 3>  RegistrationImageType;
  typedef itk::Image<IOPixelType, 3>  IOImageType;
  
  itkGetConstMacro( SampleNum, int );
  itkSetMacro( SampleNum, int );
  
  itkGetConstMacro( ImageSize0, int );
  itkSetMacro( ImageSize0, int );

  itkGetConstMacro( ImageSize1, int );
  itkSetMacro( ImageSize1, int );

  itkGetConstMacro( ImageSize2, float );
  itkSetMacro( ImageSize2, float );

  itkGetConstMacro( MetricValue0, float );
  itkSetMacro( MetricValue0, float );

  itkGetConstMacro( MetricValue1, float );
  itkSetMacro( MetricValue1, float );

  itkGetConstMacro( NumBin, unsigned int );
  itkSetMacro( NumBin, unsigned int );

  itkGetConstMacro( MaxFixed, float );
  itkSetMacro( MaxFixed, float );

  itkGetConstMacro( MaxMoving, float );
  itkSetMacro( MaxMoving, float );

protected:
  ProjectionLocalAlignmentDispatcher();
  ~ProjectionLocalAlignmentDispatcher();
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  void minimumInd(float *numbers,float &value,int &index,int length);

  void writeOutput(float* data,char* name,int length);
  
  void DfmfldModulation(float* im1,int* ordered,int* parents,int step1);

  void energy(float *val,int* ind,int len,float offset,int k,int* v,float* z,float* f,int* ind1);

  void regulate(float* r,int* indr,int rl,float dx,float dy,float dz);

  void interpolate3d(float* interp,float* input,float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2,bool flag);

  void readRaw(char str2[], float*& pixels,int SZ);

  void InterpDeformationFld(float* u1,float* v1,float* w1,float* u0,float* v0,float* w0,int m,int n,int o,int m2,int n2,int o2);

  void fldinv(float* ui,float* vi,float* wi,float* u,float* v,float* w,int m,int n,int o);

  void combineDeformation(float* u3,float* v3,float* w3,float* u1,float* v1,float* w1,float* u2,float* v2,float* w2,int m,int n,int o);

  void smoothDfm(float* u1,float* v1,float *w1,int m,int n,int o,int expsteps,int factor);

  void symmetricMapping(float* u,float* v,float* w,float* u2,float* v2,float* w2,int m,int n,int o,int factor);

  void genDfm(float* u1, float* v1, float* w1, float* u0, float* v0, float* w0, float* costall, float alpha, int hw, int step1, float quant, int* ordered, int* parents);

  void optimization(float* im1, float* im1b, float* costall, float alpha, int hw, float step1, float quant);

  void updateImage(float* warped,float* im1,float* im1b,float* u1,float* v1,float* w1);
  
public:

  void deformableReg(std::string modelDirectory, std::string outputDirectory, int debugLevel);

  void ResampleDeformationFld(char* controlGrid, std::string modelDirectory, const RegistrationImageType::Pointer fldX, const RegistrationImageType::Pointer fldY, const RegistrationImageType::Pointer fldZ);



private:

  ProjectionLocalAlignmentDispatcher(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  int m_SampleNum; 
  
  int m_ImageSize0; 
  int m_ImageSize1; 
  int m_ImageSize2; 

  float m_MetricValue0;
  float m_MetricValue1;

  unsigned int m_NumBin; //# of bins to compute entropy
  float m_MaxFixed; // max intensity of fixed 
  float m_MaxMoving; // max intensity of moving 
  

};

} // end namespace idp
} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "idpProjectionLocalAlignmentDispatcher.txx"
#endif

#endif
