/*=========================================================================

  idpRegistrationDispatcher.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpRegistrationDispatcher_h
#define __idpRegistrationDispatcher_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "idpImageSeries.h"
#include "idpImageSeriesUtilities.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredAffineTransform.h"
#include "idpReferenceSpaceUtilities.h"
#include "itkFixedArray.h"
#include <vector>
#include <string>
#include <sstream>

namespace itk
{
namespace idp
{

/** \class RegistrationDispatcher
 *
 */
template <class TPixel>
class RegistrationDispatcher : public Object
{
public:

  /** Standard typedefs. */
  typedef RegistrationDispatcher      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(RegistrationDispatcher, Object);
  
  /** Image series related types. */
  typedef TPixel PixelType;
  typedef ImageSeries::Pointer            ImageSeriesPointer;
  typedef std::vector<ImageSeriesPointer> ImageSeriesVector;
  typedef Specimen::Pointer               SpecimenPointer;

  /** Transform types. */
  typedef CenteredRigid2DTransform<double>        Rigid2DTransformType;
  typedef CenteredAffineTransform<double,3>       Affine3DTransformType;
  typedef CenteredAffineTransform<double,2>       Affine2DTransformType;
  typedef Rigid2DTransformType::ParametersType    ParametersType;
  typedef Affine3DTransformType::ParametersType   Affine3DParametersType;

  /** Volume types. */
  typedef ImageSeriesUtilities<PixelType>                     ImageSeriesUtilitiesType;
  typedef typename ImageSeriesUtilitiesType::VolumeType       VolumeType;
  typedef typename ImageSeriesUtilitiesType::MaskVolumeType   MaskVolumeType;
  typedef Image<PixelType,2>     ImageType;
  typedef Image<unsigned char,2> MaskImageType;

  /** Fiducials related types. */
  typedef ReferenceSpaceUtilities::PointType              PointType;
  typedef ReferenceSpaceUtilities::RegionOfInterestType   RegionOfInterestType;

  /** Other types. */
  typedef FixedArray<double,3> ArrayType;


  /** Get the number of times volume has been modified. */
  itkGetConstMacro( VolumeModificationCount, unsigned long );
 
  /** Get the number of times volume-atlas transform been modified. */
  itkGetConstMacro( TransformModificationCount, unsigned long );

  /** Volume is the 3D reconstruction of sections at downsampled 16 resolution. */
  itkGetObjectMacro( Volume, VolumeType );

  /** Mask is the corresponding 3D mask of Volume. */
  itkGetObjectMacro( Mask, MaskVolumeType );

  /** SliceMask indicates valid sections within the reconstruction. */
  itkGetObjectMacro( SliceMask, MaskVolumeType );

  /** ZProjection of the current reconstruction. */
  itkGetObjectMacro( ZProjection, ImageType );

  /** ZProjectionMask is the 2D mask associated with the z projection. */
  itkGetObjectMacro( ZProjectionMask, MaskImageType );

  /** AtlasVolume is the 3D atlas at downsampled 16 resolution. */
  itkGetObjectMacro( AtlasVolume, VolumeType );

  /** AtlasMask is the corresponding 3D mask of AtlasVolume. */
  itkGetObjectMacro( AtlasMask, MaskVolumeType );

  /** AtlasHeadMask is the mask restricted to the head region. */
  itkGetObjectMacro( AtlasHeadMask, MaskVolumeType );

  /** AtlasZProjection is the z projection of the atlas volume. */
  itkGetObjectMacro( AtlasZProjection, ImageType );
  itkGetObjectMacro( ResampledAtlasZProjection, ImageType );

  /** AtlasZProjectionMask is the 2D mask associated with the z projection. */
  itkGetObjectMacro( AtlasZProjectionMask, MaskImageType );

  /** Volume to atlas transform. */
  itkGetObjectMacro( VolumeToAtlasTransform, Affine3DTransformType );

  /** Atlas to volume transform. */
  itkGetObjectMacro( AtlasToVolumeTransform, Affine3DTransformType );

  /** Zprojection volume to atlas transform. */
  itkGetObjectMacro( ZProjectionTransform, Affine2DTransformType );

  /** Atlas to grid transform. */
  itkGetObjectMacro( GridToAtlasTransform, Affine3DTransformType );
  itkGetObjectMacro( AtlasToGridTransform, Affine3DTransformType );

  /** Set/Get the Zprojection metric. */
  itkSetMacro( ZMetric, double );
  itkGetConstMacro( ZMetric, double );
  itkGetConstMacro( ZMetricBackup, double );

  /** Set/Get the registration metric. */
  itkSetMacro( Metric, double );
  itkGetConstMacro( Metric, double );
  itkGetConstMacro( MetricBackup, double );

  /** Set/Get verbose flag. */
  itkSetMacro( Verbose, bool );
  itkGetConstMacro( Verbose, bool );

  /** Get combined series pointer. */
  itkGetObjectMacro( Series, ImageSeries );
  itkGetConstObjectMacro( Series, ImageSeries );

  /** Get image series vector. */
  ImageSeriesVector & GetSeriesVector()
    {
    return m_SeriesVector;
    }

  /** Get specimen pointer. */
  itkGetObjectMacro( Specimen, Specimen );
  itkGetConstObjectMacro( Specimen, Specimen );

  
  /** Get image series utility. */
  itkGetObjectMacro( Utility, ImageSeriesUtilitiesType );
  itkGetConstObjectMacro( Utility, ImageSeriesUtilitiesType );

  /** Get reference space fiducials. */
  itkGetObjectMacro( Fiducials, ReferenceSpaceUtilities );

  /** Set model directory. */
  virtual void SetModelDirectory( const char * dir );
  itkGetStringMacro( ModelDirectory );

  /** Load image series information from xml file. */
  virtual void LoadImageSeriesFromXML( const char * fileName );

  /** Write transform data to xml file. */
  virtual void WriteTransformsToXML( const char * fileName );

  /** Read transform data from xml file. */
  virtual void LoadTransformsFromXML( const char * fileName );
  virtual void LoadEmbeddedTransformsFromXML( const char * fileName ); 

  /** Initialize data before registration. */
  virtual void Initialize();

  /** Backup and rollback transform parameters. */
  virtual void BackupImageTransforms();
  virtual void RollbackImageTransforms();

  /** Populate volume from current parameters. */
  virtual void PopulateVolume();
  virtual void PopulateVolume(double redScale, double greenScale, double blueScale);

  /** Write volume to output directory. */
  virtual void WriteVolume( const char * dir, bool useCount = false );

  /** Read in atlas data. */
  virtual void ReadAtlasData();
  virtual void ReadAtlasZProjectionMask();

  /** Resample atlas z projection. */
  virtual void ResampleAtlasZProjection();

  /** Resample atlas data. */
  virtual void ResampleAtlas();

  /** Write current resampled atlas data to output directory. */
  virtual void WriteResampledAtlas( const char * dir, bool useCount = false );
  virtual void WriteResampledVolume( const char * dir, bool useCount = false );
  
  /** Compute z metric. */
  virtual void ComputeZMetric();

  /** Compute registration metric. */
  virtual void ComputeMetric();

  /** Compute section based metric. */
  virtual void ComputeSectionBasedMetric();

  /** Registration checkpoint and rollback. */
  virtual void ZProjectionCheckpoint();
  virtual bool ZProjectionRollback(bool force = false, bool checkpoint = true, double threshold = 0.0 );

  /** Registration checkpoint and rollback. */
  virtual void Checkpoint();
  virtual bool Rollback( bool force = false, bool checkpoint = true, double threshold = 0.0 );

  /** Resample volume to atlas space. */
  virtual void ResampleVolume();

  /** Extract view slice from resampled volume. */
  virtual void ExtractViewSlice( typename VolumeType::Pointer & volume, unsigned int dim, typename ImageType::Pointer & image );
  virtual void ExtractViewSlice( unsigned int dim, typename ImageType::Pointer & image );

  /** Set volume to atlas transform. */
  virtual void SetVolumeToAtlasTransform( Affine3DTransformType * t )
    {
    m_VolumeToAtlasTransform = t;
    m_AtlasToVolumeTransform->SetCenter( m_VolumeToAtlasTransform->GetCenter() );
    m_VolumeToAtlasTransform->GetInverse( m_AtlasToVolumeTransform );
    m_TransformModificationCount++;
    }

  /** Compose transform with z-axis flip. */
  virtual void ComposeWithZAxisFlip();

  /** Set/Get the red channel scale. */
  itkSetMacro( RedScale, double );
  itkGetConstMacro( RedScale, double );

  /** Set/Get the green channel scale. */
  itkSetMacro( GreenScale, double );
  itkGetConstMacro( GreenScale, double );

  /** Set/Get the blue channel scale. */
  itkSetMacro( BlueScale, double );
  itkGetConstMacro( BlueScale, double );

  /** Set/Get the image intensity background. */
  itkSetMacro( ImageBackground, PixelType );
  itkGetConstMacro( ImageBackground, PixelType );

  /** Set/Get the intensity inverse flag */
  itkSetMacro( InvertIntensity, bool );
  itkGetConstMacro( InvertIntensity, bool );

  /** Read color scale from xml file. */
  virtual void LoadColorScalesFromXML( const char * fileName );

  /** Compose with atlas-to-grid transform. */
  virtual void ComposeWithAtlasToGridTransform();
  virtual void ComposeWithGridToAtlasTransform();


protected:
  RegistrationDispatcher();
  ~RegistrationDispatcher();
  void PrintSelf(std::ostream& os, Indent indent) const;


  bool                             m_Verbose;
  ImageSeriesVector                m_SeriesVector;
  ImageSeriesPointer               m_Series;
  SpecimenPointer                  m_Specimen;
  std::string                      m_ModelDirectory;

  typename ImageSeriesUtilitiesType::Pointer    m_Utility;
  ReferenceSpaceUtilities::Pointer              m_Fiducials;

  unsigned long                    m_VolumeModificationCount;
  unsigned long                    m_TransformModificationCount;

  typename VolumeType::Pointer              m_Volume; 
  typename MaskVolumeType::Pointer              m_Mask;
  typename MaskVolumeType::Pointer              m_SliceMask;
  typename ImageType::Pointer               m_ZProjection;
  typename MaskImageType::Pointer               m_ZProjectionMask;

  Affine3DTransformType::Pointer   m_VolumeToAtlasTransform;
  Affine3DTransformType::Pointer   m_AtlasToVolumeTransform;

  Affine3DTransformType::Pointer   m_GridToAtlasTransform;
  Affine3DTransformType::Pointer   m_AtlasToGridTransform;

  Affine2DTransformType::Pointer   m_ZProjectionTransform;

  typename VolumeType::Pointer              m_AtlasVolume;
  typename MaskVolumeType::Pointer              m_AtlasMask;
  typename MaskVolumeType::Pointer              m_AtlasHeadMask;
  typename ImageType::Pointer               m_AtlasZProjection;
  typename MaskImageType::Pointer               m_AtlasZProjectionMask;

  typename VolumeType::Pointer              m_ResampledAtlasVolume;
  typename MaskVolumeType::Pointer              m_ResampledAtlasMask;
  typename MaskVolumeType::Pointer              m_ResampledAtlasHeadMask;
  typename ImageType::Pointer               m_ResampledAtlasZProjection;
  typename MaskImageType::Pointer               m_ResampledAtlasZProjectionMask;

  typename VolumeType::Pointer              m_ResampledVolume;
  typename MaskVolumeType::Pointer              m_ResampledMask;

  double                           m_ZMetric;
  double                           m_ZMetricBackup;
  Affine2DTransformType::Pointer   m_ZProjectionTransformBackup;

  double                           m_Metric;
  double                           m_MetricBackup;
  Affine3DTransformType::Pointer   m_TransformBackup;

  double                           m_RedScale;
  double                           m_GreenScale;
  double                           m_BlueScale;
  PixelType                        m_ImageBackground;
  bool                             m_InvertIntensity;


private:
  RegistrationDispatcher(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


};

} // end namespace idp
} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "idpRegistrationDispatcher.txx"
#endif

#endif
