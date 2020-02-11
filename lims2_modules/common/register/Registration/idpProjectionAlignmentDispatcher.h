/*=========================================================================

  idpProjectionAlignmentDispatcher.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpProjectionAlignmentDispatcher_h
#define __idpProjectionAlignmentDispatcher_h
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

/** \class ProjectionAlignmentDispatcher
 *
 */
template <class TPixel>
class ProjectionAlignmentDispatcher : public RegistrationDispatcher<TPixel>
{
public:

  /** Standard typedefs. */
  typedef ProjectionAlignmentDispatcher  Self;
  typedef RegistrationDispatcher<TPixel>  Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(ProjectionAlignmentDispatcher, RegistrationDispatcher);
  
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

  
  /** Volumetric affine correlation registration. */
  void AffineCorrelationRegistration( const char * fname = "affineCorrelationParameters-1.xml" );

  /** Volumetric rigid versor registration, usually use as initialization */
  void RigidVersor3DRegistration( const char * fname = "rigidVersor3DParameters-1.xml" );

  /** Volumetric affine registration. */
  void CenteredAffineVolumeRegistration( const char * fname = "CenteredAffineMutualInfoVolumeRegistration-1.xml" );
  
  /** detect the software pause which causes shift between sections */
  void DetectSoftwarePause();
  
  /** compute the mean intensity of each subimage */
  void subimageMeanIntensity(std::vector<double>& subimageMeanIntensity);
  
  /** align one subimage to the other one */
  void subimage2subimageRigid2DRegistration(unsigned int fixed, unsigned int moving, double similarityDrop);
  
  /** compare the similarity of two subimages using cross-correlation 
      Note: subimages are extracted from m_Volume directly*/
  double compareSubimages(unsigned int fixedIdx, unsigned int movingIdx);
  
  /** compare the similarity of two subimages using mutual information 
      Note: subimages are extracted from m_Volume directly*/
  double compareSubimagesMI(unsigned int fixedIdx, unsigned int movingIdx, typename Superclass::ImageSeriesUtilitiesType::Transform2DType::Pointer& tran);
  
  /** compare the similarity of two subimages using mutral information 
      Note: subimages are extracted from m_Volume directly*/
  double compareSubimagesMI(unsigned int fixedIdx, unsigned int movingIdx);
  
  /** set the global alignment model directory*/
  void SetGlobalAlignmentModelDirectory(void);
  
protected:
  ProjectionAlignmentDispatcher();
  ~ProjectionAlignmentDispatcher();
  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  ProjectionAlignmentDispatcher(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace idp
} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "idpProjectionAlignmentDispatcher.txx"
#endif

#endif
