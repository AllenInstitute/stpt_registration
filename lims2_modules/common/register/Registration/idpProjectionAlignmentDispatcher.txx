/*=========================================================================

  idpProjectionAlignmentDispatcher.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpProjectionAlignmentDispatcher_txx
#define __idpProjectionAlignmentDispatcher_txx

#include "idpProjectionAlignmentDispatcher.h"
#include "tinyxml.h"
#include "idpRegistrationUtilities.h"
#include "itkSliceImageConstIterator.h"
#include "idpXMLUtilities.h"
#include "idpImageSeries.h"
#include "idpAffineCorrelationVolumeRegistration.h"
#include "idpMutualInfoVersorRigidVolumeRegistration.h"
#include "itkIdentityTransform.h"
#include "idpMutualInfoCenteredAffineRegularGradientDescentRegistration.h"
#include <sstream>

namespace itk {
namespace idp {

/**
    * Constructor
    */
template<class TPixel>
ProjectionAlignmentDispatcher<TPixel>::
ProjectionAlignmentDispatcher() {


}

/**
    * Destructor
    */
template<class TPixel>
ProjectionAlignmentDispatcher<TPixel>::
~ProjectionAlignmentDispatcher() {

}

/**
    * PrintSelf
    */
template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>::
PrintSelf(std::ostream& os, Indent indent) const {
    Superclass::PrintSelf(os, indent);

}

/**
    * Volumetric affine correlation registration
    */
template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>::
AffineCorrelationRegistration(const char * pfilename) {
    if (this->m_Verbose) {
        std::cout << "AffineCorrelationRegistration ..." << pfilename << std::endl;
    }

    // do registration in 8-bits
    typedef itk::Image<unsigned char, 3 > UCharVolumeType;
    UCharVolumeType::Pointer ucharVolume, ucharAtlas;

    double shift = 0.0;
    double scale = 0.2;
    ShiftScale<VolumeType, UCharVolumeType > (this->m_Volume, ucharVolume, shift, scale);
    ShiftScale<VolumeType, UCharVolumeType > (this->m_AtlasVolume, ucharAtlas, shift, scale);

    typedef AffineCorrelationVolumeRegistration<UCharVolumeType> RegistrationType;
    RegistrationType::Pointer reg = RegistrationType::New();

    reg->SetFixedVolume(ucharVolume);
    reg->SetFixedMask(this->m_Mask);
    reg->SetMovingVolume(ucharAtlas);
    reg->SetMovingMask(this->m_AtlasMask);

    std::string fname = this->GetModelDirectory();
    fname += "/";
    fname += pfilename;
    reg->LoadParametersFromXML(fname.c_str());

    reg->SetInputTransform(this->GetVolumeToAtlasTransform());
    reg->Compute();

    ParametersType q = reg->GetOutputTransform()->GetParameters();
    this->m_VolumeToAtlasTransform->SetParameters(q);

    typedef DecomposeAffineMatrix3D DecomposerType;
    DecomposerType::Pointer decomposer = DecomposerType::New();
    decomposer->SetInputMatrix(this->m_VolumeToAtlasTransform->GetMatrix());
    decomposer->SetUseShear(false);
    decomposer->Compute();

    this->m_VolumeToAtlasTransform->SetMatrix(decomposer->GetOutputMatrix());
    this->m_AtlasToVolumeTransform->SetCenter(this->m_VolumeToAtlasTransform->GetCenter());
    this->m_VolumeToAtlasTransform->GetInverse(this->m_AtlasToVolumeTransform);

    this->m_TransformModificationCount++;

    if (this->m_Verbose) {
        std::cout << " Transform: " << this->m_VolumeToAtlasTransform->GetParameters();
        std::cout << " [done]" << std::endl << std::endl;
    }
}

/**
    * Volumetric rigid versor registration.
    */
template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>::
RigidVersor3DRegistration(const char * pfilename) {
    std::cout << "RigidVersor3DRegistration ..." << pfilename << std::endl;

    //prepare registration
    typedef itk::Image<TPixel, 3 > VolumeType;
    typedef itk::idp::MutualInfoVersorRigidVolumeRegistration< VolumeType,VolumeType > RegistrationType;
    typename RegistrationType::Pointer reg = RegistrationType::New();

    std::string fname = this->GetModelDirectory();
    fname += "/";
    fname += pfilename;
    reg->LoadParametersFromXML(fname.c_str());   

    reg->SetFixedInputImage(this->m_Volume);
    reg->SetMovingInputImage(this->m_AtlasVolume);
        
    // do registration
    if (this->m_Verbose) {
            reg->SetVerbose(true);
        }        
    reg->StartRegistration();

    //after registration
    this->m_VolumeToAtlasTransform->SetCenter(reg->GetOutputTransform()->GetCenter());
    this->m_VolumeToAtlasTransform->SetTranslation(reg->GetOutputTransform()->GetTranslation());
    this->m_VolumeToAtlasTransform->SetMatrix(reg->GetOutputTransform()->GetMatrix());

    this->m_AtlasToVolumeTransform->SetCenter(this->m_VolumeToAtlasTransform->GetCenter());
    this->m_VolumeToAtlasTransform->GetInverse(this->m_AtlasToVolumeTransform);

    this->m_TransformModificationCount++;

    std::cout << " Transform: " << this->m_VolumeToAtlasTransform->GetParameters();
    std::cout << " [done]" << std::endl << std::endl;
}


template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>::
CenteredAffineVolumeRegistration(const char * pfilename) 
{
    std::cout << "CenteredAffineMutualInfoVolumeRegistration ..." << pfilename << std::endl;

    //prepare registration
    typedef itk::Image<TPixel, 3 > VolumeType;
    typedef itk::idp::MutualInfoCenteredAffineRegularGradientDescentRegistration< VolumeType, VolumeType > RegistrationType;
    typename RegistrationType::Pointer reg = RegistrationType::New();

    std::string fname = this->GetModelDirectory();
    fname += "/";
    fname += pfilename;
    reg->LoadParametersFromXML( fname.c_str() );

    reg->SetFixedInputImage(this->m_Volume);
    reg->SetMovingInputImage(this->m_AtlasVolume);
    //reg->SetInitialTransformParameters( this->GetVolumeToAtlasTransform()->GetParameters() );
    reg->SetInputTransform( this->GetVolumeToAtlasTransform() );

    //do registration
    if (this->m_Verbose) {
            reg->SetVerbose(true);
        }        
    reg->StartRegistration();

    //after registration
    Affine3DParametersType q = reg->GetOutputTransform()->GetParameters();
    this->m_VolumeToAtlasTransform->SetParameters(q);

    this->m_VolumeToAtlasTransform->SetCenter(reg->GetOutputTransform()->GetCenter());
    this->m_VolumeToAtlasTransform->SetTranslation(reg->GetOutputTransform()->GetTranslation());
    this->m_VolumeToAtlasTransform->SetMatrix(reg->GetOutputTransform()->GetMatrix());

    this->m_AtlasToVolumeTransform->SetCenter(this->m_VolumeToAtlasTransform->GetCenter());
    this->m_VolumeToAtlasTransform->GetInverse(this->m_AtlasToVolumeTransform);

    this->m_TransformModificationCount++;

    std::cout << " Transform: " << this->m_VolumeToAtlasTransform->GetParameters();
    std::cout << " [done]" << std::endl << std::endl;
}

/** detect the software pause which causes shift between sections */
template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>:: 
DetectSoftwarePause()
{
	std::cout<<"Detecting the software pause ................................................."<<std::endl;
	
	// for debug only --------------------------
	//typedef ImageFileWriter<typename	Superclass::ImageSeriesUtilitiesType::ImageType > WriterType;
	//typename WriterType::Pointer writer = WriterType::New();
	
	//some parameters need to tweek
	double leastDropOfSimilarity = 0.15; //the threshold indicating a suddenly change of similarity
	unsigned searchRange = 3;            //range within which the similarity recover again.
	if (this->m_Specimen->GetId()==714144)
	{
		searchRange = 5; //714144 is bigger than usual. So needs special treatment
		std::cout<<"software pause with larger gap detected!"<<std::endl;
	}

	double leastGainOfSimilarity = 0.8; //the threshold (in percentage) of the least recovery of similarity
	
	ImageSeries::SubImageArray subimages = Superclass::m_Series->GetSubImages();
	std::vector<double> similarityVect (subimages.size(),0);   //vector contains the similarity between each neighboring subimage
	std::vector<double> imgEnergy (subimages.size(),0);   //vector contains energy of each subimage
	
	std::cout<<"There are "<<subimages.size()<<" subimages."<<std::endl;
	
	// get the mean intensity of each subimage
	this->subimageMeanIntensity(imgEnergy);
	
	for ( unsigned int k = 0; k < subimages.size()-1; k++ )
	//for ( unsigned int k = 70; k < 80; k++ )
	{
        similarityVect[k]=compareSubimages(k, k+1);
		std::cout<<k<<": "<<similarityVect[k]<<" "<<imgEnergy[k]<<std::endl;
	}
	
	// search for the software pause
	double biggestSimilarityDecrease = 0.0;
	double biggestIntensityDecrease = 0.0;
	double highestSimilarity = 0.0;
	double highestIntensity = 0.0;
	unsigned int biggestSimilarityDecreaseIndex = 0;
	unsigned int biggestIntensityDecreaseIndex = 0;
	for ( unsigned int k = 1; k < subimages.size()-1; k++ )
	{
		if ((similarityVect[k-1]-similarityVect[k])>biggestSimilarityDecrease)
		{
			biggestSimilarityDecrease = similarityVect[k-1]-similarityVect[k];
			biggestSimilarityDecreaseIndex = k;
		}
		
		if (similarityVect[k]>highestSimilarity)
		{
			highestSimilarity = similarityVect[k];
		}
		
		if (imgEnergy[k]>highestIntensity)
		{
			highestIntensity = imgEnergy[k];
		}
		
		if ((imgEnergy[k-1]-imgEnergy[k])>biggestIntensityDecrease)
		{
			biggestIntensityDecrease = imgEnergy[k-1]-imgEnergy[k];
			biggestIntensityDecreaseIndex = k;
		}
	}
	std::cout<<"The largest drop of similarity measure on section "<<biggestSimilarityDecreaseIndex<<": "<<biggestSimilarityDecrease<<" Highest: "<<highestSimilarity<<std::endl;
	std::cout<<"The largest drop of intensity measure on section "<<biggestIntensityDecreaseIndex<<": "<<biggestIntensityDecrease<<" Highest: "<<highestIntensity<<std::endl;
	
	//	the drop is not big enough
	if (biggestSimilarityDecrease<leastDropOfSimilarity)
	{
		std::cout<<"The drop is not significant enough. Probably, it's not a software pause."<<std::endl;
		return;
	}
	
	// check if within a few steps the similarity measure recovers
	double biggestSimilarityIncrease = 0.0;
	unsigned int biggestSimilarityIncreaseIndex = 0;
	for ( unsigned int k = biggestSimilarityDecreaseIndex; (k < biggestSimilarityDecreaseIndex+searchRange)&&(k < subimages.size()-1); k++ )
	{
		if ((similarityVect[k+1]-similarityVect[k])>biggestSimilarityIncrease)
		{
			biggestSimilarityIncrease = similarityVect[k+1]-similarityVect[k];
			biggestSimilarityIncreaseIndex = k+1;
		}
	}

	if (biggestSimilarityIncrease/biggestSimilarityDecrease>leastGainOfSimilarity)
	{
		std::cout<<"The drop of similarity re-gained at section "<<biggestSimilarityIncreaseIndex<<". Looks like a software pause."<<std::endl;
	}
	else
	{
		std::cout<<"The drop of similarity is not recovered shortly. Probably, it's not a software pause."<<std::endl;
		return;
	}
	
	// perform 2D registratio to correct the software pause error
	for (unsigned int i=biggestSimilarityDecreaseIndex;i<biggestSimilarityIncreaseIndex;i++)
	{
		if (imgEnergy[i+1]>0.0)
		{
		std::cout<<"Register section "<<i+1<<" to section "<<biggestSimilarityDecreaseIndex<<"......."<<std::endl;
		this->subimage2subimageRigid2DRegistration(biggestSimilarityDecreaseIndex, i+1, biggestSimilarityDecrease);
		}
	}
}

/** align one subimage to the other one */
template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>::
subimage2subimageRigid2DRegistration(unsigned int fixed, unsigned int moving, double similarityDrop)
{
  // do 2d registration 	
  typedef float                                      InternalPixelType;
  typedef Image< InternalPixelType, 2 >              InternalImageType;
  
  typedef typename Superclass::ImageSeriesUtilitiesType::Transform2DType          SubImgTransformType;
  typedef RegularStepGradientDescentOptimizer                                           OptimizerType;  
  typedef LinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
  
  typedef NormalizedCorrelationImageToImageMetric< InternalImageType, InternalImageType >    MetricType;
  typedef OptimizerType::ScalesType       OptimizerScalesType;

  typedef MultiResolutionImageRegistrationMethod<  InternalImageType,  InternalImageType    > RegistrationType;
  typedef RecursiveMultiResolutionPyramidImageFilter< InternalImageType, InternalImageType  > ImagePyramidType;

  typename OptimizerType::Pointer      optimizer     = OptimizerType::New();
  typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  typename RegistrationType::Pointer   registration  = RegistrationType::New();
  typename MetricType::Pointer         metric        = MetricType::New();

  // prepare the input images
  ImageSeries::SubImageArray subimages = Superclass::m_Series->GetSubImages();  
  typename SubImgTransformType::Pointer   rigidTran  = SubImgTransformType::New();
  registration->SetTransform( rigidTran );
  
  typename ImagePyramidType::Pointer fixedImagePyramid = ImagePyramidType::New();
  typename ImagePyramidType::Pointer movingImagePyramid = ImagePyramidType::New();
  
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  registration->SetMetric( metric );
  registration->SetNumberOfLevels( 3 );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );

  typename InternalImageType::Pointer fixedImg, movingImg;
  ExtractSlice<typename Superclass::VolumeType, InternalImageType>( this->m_Volume, fixedImg, fixed, 2 );
  ExtractSlice<typename Superclass::VolumeType, InternalImageType>( this->m_Volume, movingImg, moving, 2 );
  
  // compute the similarity before registration
  double similarityBefore = compareSubimages(fixed, moving);
  double similarityBeforeMI = compareSubimagesMI(fixed, moving);
  
  registration->SetFixedImage( fixedImg );
  registration->SetMovingImage( movingImg );
  registration->SetFixedImageRegion( fixedImg->GetBufferedRegion() );
  
  typedef CenteredTransformInitializer< typename Superclass::ImageSeriesUtilitiesType::Transform2DType,  InternalImageType, InternalImageType >  TransformInitializerType;
  typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();

  initializer->SetTransform( rigidTran );
  initializer->SetFixedImage( fixedImg );
  initializer->SetMovingImage( movingImg );
  initializer->MomentsOn();
  initializer->GeometryOn();
  initializer->InitializeTransform();
  rigidTran->SetAngle( 0.0 );
  registration->SetInitialTransformParameters( rigidTran->GetParameters() );
  
  OptimizerScalesType optimizerScales( rigidTran->GetNumberOfParameters() );
  const float sizeRotation = 1;
  const unsigned long sizeTranslation = 1000000;
  
  optimizerScales[0] = 1.0 / sizeRotation; //size
  optimizerScales[1] = 1.0 / sizeTranslation; //rotation
  optimizerScales[2] = 1.0 / sizeTranslation; //center
  optimizerScales[3] = 1.0 / sizeTranslation; 
  optimizerScales[4] = 1.0 / sizeTranslation; //translation
  
  optimizer->SetScales( optimizerScales );

  //metric->SetNumberOfHistogramBins( 128 );
  //metric->SetNumberOfSpatialSamples( 25000 );
  //metric->ReinitializeSeed( 76926294 );
 
	unsigned long NumberOfIterations=180;
	float RelaxationFactor=0.8;
	float MaxStepLength=3.0;
	float MinStepLength=0.001;
	
	if (this->m_Specimen->GetId()==714144)
	{
		NumberOfIterations=300.0; //714144 is bigger than usual. So needs special treatment
		MinStepLength=0.0001;
		std::cout<<"software pause with larger gap detected!"<<std::endl;
	}
 
  optimizer->SetNumberOfIterations(  NumberOfIterations  );
  optimizer->SetRelaxationFactor( RelaxationFactor );
  optimizer->SetMaximumStepLength( MaxStepLength );  
  optimizer->SetMinimumStepLength( MinStepLength );

  // Setup an optimizer observer
  //CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  //optimizer->AddObserver( itk::IterationEvent(), observer );
  
  try 
    { 
    registration->StartRegistration(); 
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return;
    } 
  
  // get the similarity measure after registration
  double similarityAfter = -optimizer->GetValue();
  double similarityAfterMI = compareSubimagesMI(fixed, moving, rigidTran);
  
  std::cout<<"Similarity before:after --- "<<similarityBefore<<":"<<similarityAfter<<std::endl;
  std::cout<<"SimilarityMI before:after --- "<<similarityBeforeMI<<":"<<similarityAfterMI<<std::endl;
  //if ((((similarityAfter-similarityBefore)/similarityDrop>0.9)&&((moving-fixed)<2))||(((moving-fixed)>=2)&&(similarityAfter/similarityBefore>1.2)))
  float thd1=1.1;
  float thd2=0.1;
	if (this->m_Specimen->GetId()==714144)
	{
		thd2=0.07; //714144 is bigger than usual. So needs special treatment
		std::cout<<"software pause with larger gap detected!"<<std::endl;
	}
	
  if ((similarityAfter/similarityBefore>thd1)&&(similarityAfterMI-similarityBeforeMI>thd2))
  {
  std::cout<<"The sudden drop of similarity CAN be corrected by regitration. It should be a software pause. Updating the subimage transform ..."<<std::endl;
  // update the transform of each 2d subimage
  typename SubImgTransformType::ParametersType finalParameters = registration->GetLastTransformParameters();
  typename SubImgTransformType::Pointer finalTransform = SubImgTransformType::New();
  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( rigidTran->GetFixedParameters() );  
  for ( unsigned int k = fixed+1; k < subimages.size(); k++ )
  {
	typename Superclass::ImageSeriesUtilitiesType::TransformType::Pointer output;
    //TransformUtilities::Compose( finalTransform, subimages[k]->GetTransform(), output );
	TransformUtilities::Compose( subimages[k]->GetTransform(), finalTransform, output );
    subimages[k]->SetTransform( output );
  }
  }
  else
  {
	std::cout<<"The sudden drop of similarity can not be corrected by regitration. It should NOT be a software pause."<<std::endl;
  }
}

/** compute the mean intensity of each subimage (to see if it is totally/almost dark) */
template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>::
subimageMeanIntensity(std::vector<double>& subimageMeanIntensity)
{
	ImageSeries::SubImageArray subimages = Superclass::m_Series->GetSubImages(); //to know how many subimages we have
	
	typedef ImageRegionConstIterator< typename Superclass::VolumeType > ConstIteratorType;
	typename Superclass::VolumeType::RegionType region;
	typename Superclass::VolumeType::RegionType::IndexType start;
	typename Superclass::VolumeType::RegionType::SizeType size;
	
	typename Superclass::VolumeType::RegionType::SizeType volumeSize = this->m_Volume->GetLargestPossibleRegion().GetSize();
	
	size[0] = volumeSize[0];
	size[1] = volumeSize[1];
	size[2] = 1; // just take one section a time
	
	for (unsigned int i=0;i<subimages.size();i++)
	{
		start[0] = 0;
		start[1] = 0;
		start[2] = i;

		region.SetSize( size );
		region.SetIndex( start );
		
		ConstIteratorType It( this->m_Volume, region );
		double sum = 0.0;
		for ( It.GoToBegin(); !It.IsAtEnd();++It)
		{
			sum+=It.Get();
		}
		subimageMeanIntensity[i] = sum/volumeSize[0]/volumeSize[1];
	}
}

/** compare the similarity of two subimages using cross-correlation */
template<class TPixel>
double
ProjectionAlignmentDispatcher<TPixel>::
compareSubimages(unsigned int fixedIdx, unsigned int movingIdx)
{
	// for each two adjacent subimages, compute the similarity metric
	typedef NormalizedCorrelationImageToImageMetric<typename Superclass::ImageSeriesUtilitiesType::ImageType,typename Superclass::ImageSeriesUtilitiesType::ImageType > MMIMetricType;
	typename MMIMetricType::Pointer metricMMI = MMIMetricType::New();
	
	typedef LinearInterpolateImageFunction<typename Superclass::ImageSeriesUtilitiesType::ImageType,double > InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	
	typename Superclass::ImageSeriesUtilitiesType::ImageType::Pointer fixed, moving;
	ExtractSlice<typename Superclass::VolumeType,typename Superclass::ImageSeriesUtilitiesType::ImageType>( this->m_Volume, fixed, fixedIdx, 2 );
	ExtractSlice<typename Superclass::VolumeType,typename Superclass::ImageSeriesUtilitiesType::ImageType>( this->m_Volume, moving, movingIdx, 2 );
        
	metricMMI->SetFixedImage(fixed);
	metricMMI->SetMovingImage(moving);
	interpolator->SetInputImage(moving);
		
	typename Superclass::ImageSeriesUtilitiesType::Transform2DType::Pointer initTran = Superclass::ImageSeriesUtilitiesType::Transform2DType::New();
	typename Superclass::ImageSeriesUtilitiesType::Transform2DType::ParametersType p( initTran->GetNumberOfParameters() );
	p.Fill(0.0);	
	initTran->SetParameters( p );
		
	metricMMI->SetTransform(initTran);
	metricMMI->SetFixedImageRegion(fixed->GetLargestPossibleRegion());
	metricMMI->SetInterpolator(interpolator);
	metricMMI->Initialize();
	return -metricMMI->GetValue(initTran->GetParameters());
}

/** compare the similarity of two subimages using mutral information */
template<class TPixel>
double
ProjectionAlignmentDispatcher<TPixel>::
compareSubimagesMI(unsigned int fixedIdx, unsigned int movingIdx)
{
	// for each two adjacent subimages, compute the similarity metric
	typedef MattesMutualInformationImageToImageMetric<typename Superclass::ImageSeriesUtilitiesType::ImageType,typename Superclass::ImageSeriesUtilitiesType::ImageType > MMIMetricType;
	typename MMIMetricType::Pointer metricMMI = MMIMetricType::New();
	
	typedef LinearInterpolateImageFunction<typename Superclass::ImageSeriesUtilitiesType::ImageType,double > InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	
	typename Superclass::ImageSeriesUtilitiesType::ImageType::Pointer fixed, moving;
	ExtractSlice<typename Superclass::VolumeType,typename Superclass::ImageSeriesUtilitiesType::ImageType>( this->m_Volume, fixed, fixedIdx, 2 );
	ExtractSlice<typename Superclass::VolumeType,typename Superclass::ImageSeriesUtilitiesType::ImageType>( this->m_Volume, moving, movingIdx, 2 );
        
	metricMMI->SetFixedImage(fixed);
	metricMMI->SetMovingImage(moving);
	interpolator->SetInputImage(moving);
		
	typename Superclass::ImageSeriesUtilitiesType::Transform2DType::Pointer initTran = Superclass::ImageSeriesUtilitiesType::Transform2DType::New();
	typename Superclass::ImageSeriesUtilitiesType::Transform2DType::ParametersType p( initTran->GetNumberOfParameters() );
	p.Fill(0.0);	
	initTran->SetParameters( p );
	
	metricMMI->ReinitializeSeed( 76926294 );
	metricMMI->SetNumberOfHistogramBins(64);
	metricMMI->SetNumberOfSpatialSamples(100000);
		
	metricMMI->SetTransform(initTran);
	metricMMI->SetFixedImageRegion(fixed->GetLargestPossibleRegion());
	metricMMI->SetInterpolator(interpolator);
	metricMMI->Initialize();
	return -metricMMI->GetValue(initTran->GetParameters());
}

/** compare the similarity of two subimages using cross-correlation 
      Note: subimages are extracted from m_Volume directly*/
template<class TPixel>
double
ProjectionAlignmentDispatcher<TPixel>::
compareSubimagesMI(unsigned int fixedIdx, unsigned int movingIdx, typename Superclass::ImageSeriesUtilitiesType::Transform2DType::Pointer& tran)
{
	// for each two adjacent subimages, compute the similarity metric
	typedef MattesMutualInformationImageToImageMetric<typename Superclass::ImageSeriesUtilitiesType::ImageType,typename Superclass::ImageSeriesUtilitiesType::ImageType > MMIMetricType;
	typename MMIMetricType::Pointer metricMMI = MMIMetricType::New();
	
	typedef LinearInterpolateImageFunction<typename Superclass::ImageSeriesUtilitiesType::ImageType,double > InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	
	typename Superclass::ImageSeriesUtilitiesType::ImageType::Pointer fixed, moving;
	ExtractSlice<typename Superclass::VolumeType,typename Superclass::ImageSeriesUtilitiesType::ImageType>( this->m_Volume, fixed, fixedIdx, 2 );
	ExtractSlice<typename Superclass::VolumeType,typename Superclass::ImageSeriesUtilitiesType::ImageType>( this->m_Volume, moving, movingIdx, 2 );
        
	metricMMI->SetFixedImage(fixed);
	metricMMI->SetMovingImage(moving);
	interpolator->SetInputImage(moving);
	
	metricMMI->ReinitializeSeed( 76926294 );
	metricMMI->SetNumberOfHistogramBins(64);
	metricMMI->SetNumberOfSpatialSamples(100000);
		
	metricMMI->SetTransform(tran);
	metricMMI->SetFixedImageRegion(fixed->GetLargestPossibleRegion());
	metricMMI->SetInterpolator(interpolator);
	metricMMI->Initialize();
	return -metricMMI->GetValue(tran->GetParameters());
}

/** Set global alignment model directory*/
template<class TPixel>
void
ProjectionAlignmentDispatcher<TPixel>::
SetGlobalAlignmentModelDirectory(void)
{
	this->m_ModelDirectory += "/alignment/";
}

} // end namespace idp
} //end namespace itk

#endif
