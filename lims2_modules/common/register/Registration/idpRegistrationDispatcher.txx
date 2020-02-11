/*=========================================================================

  idpRegistrationDispatcher.txx
  
  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpRegistrationDispatcher_txx
#define __idpRegistrationDispatcher_txx

#include "idpRegistrationDispatcher.h"
#include "tinyxml.h"
#include "idpRegistrationUtilities.h"
#include "itkSliceImageConstIterator.h"
#include "idpXMLUtilities.h"
#include "idpImageSeries.h"

#include <sstream>

namespace itk
{
namespace idp
{

/**
 * Constructor
 */
template<class TPixel>
RegistrationDispatcher<TPixel>::
RegistrationDispatcher()
{

  m_Verbose = false;
  m_SeriesVector.resize( 0 );
  m_Series = 0;
  m_Specimen = 0;
  m_ModelDirectory = "";
  m_Utility = 0;
  m_Fiducials = 0;

  m_VolumeModificationCount = 0;
  m_TransformModificationCount = 0;

  m_VolumeToAtlasTransform = Affine3DTransformType::New();
  m_AtlasToVolumeTransform = Affine3DTransformType::New();
  m_TransformBackup        = Affine3DTransformType::New();

  m_VolumeToAtlasTransform->SetIdentity();
  m_AtlasToVolumeTransform->SetIdentity();
  m_TransformBackup->SetIdentity();

  m_ZProjectionTransform = Affine2DTransformType::New();
  m_ZProjectionTransform->SetIdentity();
  m_ZProjectionTransformBackup = Affine2DTransformType::New();
  m_ZProjectionTransformBackup->SetParameters( m_ZProjectionTransform->GetParameters() );

  m_ZMetric = 0.0;
  m_Metric = 0.0;

  // default for HP yellow counterstain
  m_RedScale = 0.0;
  m_GreenScale = 0.0;
  m_BlueScale = 1.0;
  m_ImageBackground = 255;
  m_InvertIntensity = true;

  // allow the final grid to be different to atlas space
  m_AtlasToGridTransform = Affine3DTransformType::New();
  m_GridToAtlasTransform = Affine3DTransformType::New();

  m_AtlasToGridTransform->SetIdentity();
  m_GridToAtlasTransform->SetIdentity();
 
}

/**
 * Destructor
 */
template<class TPixel>
RegistrationDispatcher<TPixel>::
~RegistrationDispatcher()
{

}

/**
 * PrintSelf
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << "Verbose: " << this->GetVerbose() << std::endl;
  os << "#Series: " << m_SeriesVector.size() << std::endl;
  os << "Series: " << this->GetSeries() << std::endl;
  os << "ModelDirectory: " << this->GetModelDirectory() << std::endl;
  os << "Utility: " << this->GetUtility() << std::endl;
  os << "VolumeModificationCount: " << this->GetVolumeModificationCount() << std::endl;
  os << "VolumeToAtlasTransform: " << m_VolumeToAtlasTransform << std::endl;
  os << "AtlasToVolumeTransform: " << m_AtlasToVolumeTransform << std::endl;
  os << "TransformModificationCount: " << this->GetTransformModificationCount() << std::endl;
  os << "ZMetric: " << this->GetZMetric() << std::endl;
  os << "Metric: " << this->GetMetric() << std::endl;
  os << "RedScale: " << this->GetRedScale() << std::endl;
  os << "GreenScale: " << this->GetGreenScale() << std::endl;
  os << "BlueScale: " << this->GetBlueScale() << std::endl;
  os << "GridToAtlasTransform: " << m_GridToAtlasTransform << std::endl;

}


/**
 * Load image series from xml file
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
LoadImageSeriesFromXML( const char * xmlFile )
{

  /*
   * Read image series information from XML file using tinyxml
   */
  m_Specimen = 0;
  m_Series = 0;
  m_SeriesVector.resize( 0 );

  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile << std::endl );
    }


  TiXmlNode * node;
  bool specimenModel = false;
  if ( ( node = doc.FirstChild( "specimen" ) ) )
    {
    // first try to find a specimen node  
    m_Specimen = Specimen::New();
    m_Specimen->LoadFromXML( node );
    if ( m_Verbose )
      {
      std::cout << "Loaded specimen: " << m_Specimen->GetId() << std::endl;
      std::cout << "#Series: " << m_Specimen->GetImageSeries().size() << std::endl;
      }
    specimenModel = true;
    }
  else if ( ( node = doc.FirstChild( "image-series" ) ) )
    {
    // else try an image-series node
    m_Series = ImageSeries::New();
    m_Series->LoadFromXML( node );
    if ( m_Verbose )
      {
      std::cout << "Loaded imageseries: " << m_Series->GetId() << std::endl;
      }
    m_SeriesVector.push_back( m_Series );
    }
  else if ( ( node = doc.FirstChild( "specimen-image" ) ) )
    {
    // else try an image-series node
    m_Series = ImageSeries::New();
    m_Series->LoadFromXML( node );
    if ( m_Verbose )
      {
      std::cout << "Loaded specimen-image: " << m_Series->GetId() << std::endl;
      }
    m_SeriesVector.push_back( m_Series );
    }    
   else
    {
    itkExceptionMacro( << "no specimen or image-series node found" << std::endl );
    }


  /*
   * Combine multiple image series into one
   */
  if ( specimenModel )
    {
    m_Series = ImageSeries::New();
    m_SeriesVector = m_Specimen->GetImageSeries();
    m_Series->CopyInformation( m_SeriesVector[0] );
    }

  if ( specimenModel && m_SeriesVector.size() > 1 )
    {
    m_Series->SetId( 0 );
    m_Series->ResetTreatments();
    }

  for ( unsigned int s = 0; specimenModel && s < m_SeriesVector.size(); s++ )
    {
    if ( m_SeriesVector[s]->HasTreatment( "Annotated" ) )
      {
      // skip annotated image series
      continue;
      }
    
    std::string wfs = m_SeriesVector[s]->GetWorkflowState();
    if ( wfs.compare( "processing" ) == 0 )
      {
      // skip if series is currently being processed
      continue;
      }
    
    for( unsigned int k = 0; k < m_SeriesVector[s]->GetSubImages().size(); k++ )
      {
      m_Series->InsertSubImage( m_SeriesVector[s]->GetSubImages()[k] );
      }

    if ( m_SeriesVector[s]->GetImageResolution() > m_Series->GetImageResolution() )
      {
      m_Series->SetImageResolution( m_SeriesVector[s]->GetImageResolution() );
      }

    }

  m_Series->SortSubImages();


  if ( m_Verbose )
    {
    std::cout << "#Sub images: " << m_Series->GetSubImages().size() << std::endl;
    std::cout << "Interval: " << m_Series->GetInterSectionInterval() << std::endl;
    }

  m_Utility = ImageSeriesUtilitiesType::New();
  m_Utility->SetImageSeries( m_Series );

}

/**
 * Setup model directory
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
SetModelDirectory( const char * dir )
{

  m_ModelDirectory = dir;

  if ( m_Series->GetReferenceSpace() != 0 )
    {
    m_ModelDirectory += "/";
    m_ModelDirectory += m_Series->GetReferenceSpace()->GetAge();
    m_ModelDirectory += "/";
    m_ModelDirectory += m_Series->GetPlaneOfSection();
    }
    
  m_Utility->SetModelDirectory( m_ModelDirectory.c_str() );

  if ( m_Verbose )
    {
    std::cout << "ModelDirectory: " << this->GetModelDirectory() << std::endl;
    }

}



/**
 * Initialize before registration
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
Initialize()
{
  m_Volume = 0;
  m_Mask = 0;
  m_SliceMask = 0;
  m_ZProjection = 0;
  m_ZProjectionMask = 0;
  m_VolumeModificationCount = 0;

  m_TransformModificationCount = 0;

  m_Utility->SetModelDirectory( m_ModelDirectory );
  m_Utility->ComputeVolumeMetaInformation();

  m_Fiducials = ReferenceSpaceUtilities::New();
  std::string fname = m_ModelDirectory;
  fname += "/fiducials.xml";
  if ( itksys::SystemTools::FileExists( fname.c_str(), true ) )
    {
    m_Fiducials->InsertFromXML( fname.c_str() );
    }

  fname = m_ModelDirectory;
  fname += "/colorScales.xml";
  if ( itksys::SystemTools::FileExists( fname.c_str(), true ) )
    {
    this->LoadColorScalesFromXML( fname.c_str() );
    }
  

  m_AtlasVolume = 0;
  m_AtlasMask = 0;
  m_AtlasHeadMask = 0;
  m_AtlasZProjection = 0;
  m_AtlasZProjectionMask = 0;

  m_ResampledAtlasVolume = 0;
  m_ResampledAtlasMask = 0;
  m_ResampledAtlasHeadMask = 0;
  m_ResampledAtlasZProjection = 0;
  m_ResampledAtlasZProjectionMask = 0;

  m_ResampledVolume = 0;
  m_ResampledMask = 0;

  m_VolumeToAtlasTransform->SetIdentity();
  m_AtlasToVolumeTransform->SetIdentity();
  m_TransformBackup->SetIdentity();
  m_Metric = 0.0;

  m_ZProjectionTransform->SetIdentity();
  typename ImageType::PointType c;
  c[0] = m_Utility->GetCenter()[0];
  c[1] = m_Utility->GetCenter()[1];
  m_ZProjectionTransform->SetCenter( c );
  m_ZProjectionTransformBackup->SetParameters( m_ZProjectionTransform->GetParameters() );

  m_ZMetric = 0.0;

  fname = m_ModelDirectory;
  fname += "/atlasVolumeToGrid.meta";
  if ( itksys::SystemTools::FileExists( fname.c_str(), true ) )
    {
    ReadTransform<Affine3DTransformType>( fname.c_str(), m_GridToAtlasTransform );
    }
  m_AtlasToGridTransform->SetCenter( m_GridToAtlasTransform->GetCenter() );
  m_GridToAtlasTransform->GetInverse( m_AtlasToGridTransform );

}


/**
 * Populate volume using current transform
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
PopulateVolume()
{

  if ( m_Verbose )
    {
    std::cout << "PopulateVolume ...";
    }

  m_Utility->PopulateVolume( m_Volume, m_Mask, m_SliceMask, m_RedScale, m_GreenScale, m_BlueScale, m_ImageBackground );

  if ( m_InvertIntensity )
    {
    InvertIntensity<VolumeType>( m_Volume, m_Volume );
    }

  typename VolumeType::Pointer dummy;

  MaskedAverageProjection<VolumeType,MaskVolumeType>( m_Volume, m_Mask, dummy );
  ExtractSlice<VolumeType,ImageType>( dummy, m_ZProjection, 0, 2 );

  m_VolumeModificationCount++;  

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}


/**
 * Populate volume using current transform
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
PopulateVolume(
double redScale,
double greenScale,
double blueScale )
{

  if ( m_Verbose )
    {
    std::cout << "PopulateVolume ...";
    }

  m_Utility->PopulateVolume( m_Volume, m_Mask, m_SliceMask, redScale, greenScale, blueScale, m_ImageBackground );

  if ( m_InvertIntensity )
    {
    InvertIntensity<VolumeType>( m_Volume, m_Volume );
    }

  typename VolumeType::Pointer dummy;
  MaskedAverageProjection<VolumeType,MaskVolumeType>( m_Volume, m_Mask, dummy );
  ExtractSlice<VolumeType,ImageType>( dummy, m_ZProjection, 0, 2 );

  m_VolumeModificationCount++;  

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}


/**
 * Write out current volume to output directory
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
WriteVolume(
const char * outputDirectory,
bool useCount
)
{

  if ( m_Verbose )
    {
    std::cout << "WriteVolume ...";
    }

  std::string fname;

  std::ostringstream ostr;
  ostr << m_VolumeModificationCount;

  // write out volume
  fname = outputDirectory;
  fname += "/volume";
  if ( useCount )
    {
    fname += "-";
    fname += ostr.str();
    }
  fname += ".mhd";

  WriteImage<VolumeType>( fname.c_str(), m_Volume );

  // write out mask
  fname = outputDirectory;
  fname += "/mask";
  if ( useCount )
    {
    fname += "-";
    fname += ostr.str();
    }
  fname += ".mhd";

  WriteImage<MaskVolumeType>( fname.c_str(), m_Mask );


  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}

/**
 * Backup current image transform parameters
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
BackupImageTransforms()
{
  if ( m_Verbose )
    {
    std::cout << "BackupImageTransforms ...";
    }

  m_Utility->BackupTransformParameters();

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }
}

/**
 *  rollback transform parameters
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
RollbackImageTransforms()
{
  if ( m_Verbose )
    {
    std::cout << "RevertImageTransforms ...";
    }

  m_Utility->RevertTransformParameters();

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }
}



/**
 * Read atlas data
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ReadAtlasData()
{
  if ( m_Verbose )
    {
    std::cout << "ReadAtlasData ...";
    }

  std::string fname;

  fname = m_ModelDirectory;
  fname += "/atlasVolume.mhd";
  ReadImage<VolumeType>( fname.c_str(), m_AtlasVolume );

  fname = m_ModelDirectory;
  fname += "/atlasMask.mhd";
  ReadImage<MaskVolumeType>( fname.c_str(), m_AtlasMask );

  fname = m_ModelDirectory;
  fname += "/atlasHeadMask.mhd";
  if ( itksys::SystemTools::FileExists( fname.c_str(), true ) )
    {
    ReadImage<MaskVolumeType>( fname.c_str(), m_AtlasHeadMask );
    }
  else
    {
    CastImage<MaskVolumeType,MaskVolumeType>( m_AtlasMask, m_AtlasHeadMask );
    }

  fname = m_ModelDirectory;
  fname += "/avgProj.png";
  if ( itksys::SystemTools::FileExists( fname.c_str(), true ) )
    {
    ReadImage<ImageType>( fname.c_str(), m_AtlasZProjection );
    InvertIntensity<ImageType>( m_AtlasZProjection, m_AtlasZProjection );
    }


  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }
}


/**
 * Read atlas z projection mask
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ReadAtlasZProjectionMask()
{
  if ( m_Verbose )
    {
    std::cout << "ReadAtlasZProjectionMask ...";
    }

  std::string fname;

  fname = m_ModelDirectory;
  fname += "/avgProjMask.png";

  if ( itksys::SystemTools::FileExists( fname.c_str() ) )
    {
    ReadImage<MaskImageType>( fname.c_str(), m_AtlasZProjectionMask );
    }
  else
    {
    m_AtlasZProjectionMask = 0;
    }

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }
}


/**
 * Resample the atlas z projection
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ResampleAtlasZProjection()
{

  if ( m_Verbose )
    {
    std::cout << "ResampleAtlasZProjection ...";
    }

  ResampleImage<ImageType>( m_AtlasZProjection, m_ZProjection, m_ZProjectionTransform,
                            m_ResampledAtlasZProjection, 0, "Linear" );

  if ( m_AtlasZProjectionMask.IsNotNull() )
    {
    MaskImageType::Pointer ref;
    MakeImage<MaskImageType,ImageType>( m_ZProjection, ref, 0 );
    ResampleImage<MaskImageType>( m_AtlasZProjectionMask, ref, m_ZProjectionTransform,
                                  m_ResampledAtlasZProjectionMask, 0, "NearestNeighbor" );

    }
  else
    {
    m_ResampledAtlasZProjectionMask = 0;
    }


  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}

/**
 * Compute Z metric
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ComputeZMetric()
{

  if ( m_Verbose )
    {
    std::cout << "ComputeZMetric ...";
    }

  if ( m_ResampledAtlasZProjection.IsNotNull() &&
       m_ZProjection.IsNotNull() )
    {
    ComputeCorrelation<ImageType,MaskImageType>( m_ResampledAtlasZProjection, m_ZProjection,
                                                 m_ResampledAtlasZProjectionMask, m_ZProjectionMask,
                                                 m_ZMetric );
    }
  else
    {
    m_ZMetric = 0.0;
    }

  if ( m_Verbose )
    {
    std::cout << "Metric = " << m_ZMetric;
    }

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}

/**
 * Z Projection checkpoint
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ZProjectionCheckpoint()
{

  if ( m_Verbose )
    {
    std::cout << "ZProjectionCheckpoint ...";
    }

  m_ZMetricBackup = m_ZMetric;
  m_ZProjectionTransformBackup->SetParameters( m_ZProjectionTransform->GetParameters() );

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}


/**
 * Z Projection rollback
 */
template<class TPixel>
bool
RegistrationDispatcher<TPixel>::
ZProjectionRollback(
bool force,
bool checkpoint,
double threshold )
{

  if ( m_Verbose )
    {
    std::cout << "ZProjectionRollback ...";
    }

  bool reverted = false;

  if ( force || m_ZMetric > ( m_ZMetricBackup + threshold ) )
    {
    if ( m_Verbose )
      {
      std::cout << " Reverting ";
      }
    m_ZMetric = m_ZMetricBackup;
    m_ZProjectionTransform->SetParameters( m_ZProjectionTransformBackup->GetParameters() );
    reverted = true;
    }
  else
    {
    if ( m_Verbose )
      {
      std::cout << " No Rollback ";
      }
    }

  if ( !reverted && checkpoint )
    {
    this->ZProjectionCheckpoint();
    }

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

  return reverted;

}


/**
 * registration checkpoint
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
Checkpoint()
{

  if ( m_Verbose )
    {
    std::cout << "Checkpoint ...";
    }

  m_MetricBackup = m_Metric;
  m_TransformBackup->SetParameters( m_VolumeToAtlasTransform->GetParameters() );

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}


/**
 * registration rollback
 */
template<class TPixel>
bool
RegistrationDispatcher<TPixel>::
Rollback(
bool force,
bool checkpoint,
double threshold )
{

  if ( m_Verbose )
    {
    std::cout << "Rollback ...";
    }

  std::cout << " " << (m_MetricBackup+threshold) << " " << m_Metric << std::endl;

  bool reverted = false;

  if ( force ||  ( m_Metric > ( m_MetricBackup + threshold ) ) )
    {
    if ( m_Verbose )
      {
      std::cout << " Reverting" ;
      }

    m_Metric = m_MetricBackup;
    m_VolumeToAtlasTransform->SetParameters( m_TransformBackup->GetParameters() );
    m_AtlasToVolumeTransform->SetCenter( m_VolumeToAtlasTransform->GetCenter() );
    m_VolumeToAtlasTransform->GetInverse( m_AtlasToVolumeTransform );
   
    reverted = true;
    }
  else
    {
    if ( m_Verbose )
      {
      std::cout << " No Rollback ";
      }
    }

  if ( !reverted && checkpoint )
    {
    this->Checkpoint();
    }

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

  return reverted;

}

/**
 * Resample the atlas data
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ResampleAtlas()
{

  if ( m_Verbose )
    {
    std::cout << "ResampleAtlas ...";
    }
	
	
  // resample the atlas volume
  ResampleImage<VolumeType>( m_AtlasVolume, m_Volume, m_VolumeToAtlasTransform, 
                             m_ResampledAtlasVolume, 0, "Linear" );

  // resample the atlas mask
  ResampleImage<MaskVolumeType>( m_AtlasMask, m_Mask, m_VolumeToAtlasTransform, 
                                 m_ResampledAtlasMask, 0, "NearestNeighbor" );

  // resample the atlas head mask
  ResampleImage<MaskVolumeType>( m_AtlasHeadMask, m_Mask, m_VolumeToAtlasTransform, 
                                 m_ResampledAtlasHeadMask, 0, "NearestNeighbor" );


  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}

/**
 * Resample the volume to atlas space
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ResampleVolume()
{

  if ( m_Verbose )
    {
    std::cout << "ResampleVolume ...";
    }

  // resample the volume
  ResampleImage<VolumeType>( m_Volume, m_AtlasVolume, m_AtlasToVolumeTransform, 
                             m_ResampledVolume, 0, "Linear" );

  // resample the mask
  ResampleImage<MaskVolumeType>( m_Mask, m_AtlasMask, m_AtlasToVolumeTransform, 
                                 m_ResampledMask, 0, "NearestNeighbor" );

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}


/**
 * Compute metric
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ComputeMetric()
{

  if ( m_Verbose )
    {
    std::cout << "ComputeMetric ...";
    }

  ComputeCorrelation<VolumeType,MaskVolumeType>( m_ResampledAtlasVolume, m_Volume, 
                                                 m_ResampledAtlasHeadMask, m_Mask,
                                                 m_Metric );

  if ( m_Verbose )
    {
    std::cout << "Metric = " << m_Metric;
    }

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}

/**
 * Compute section based metric
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ComputeSectionBasedMetric()
{

  if ( m_Verbose )
    {
    std::cout << "ComputeSectionBasedMetric ..." << std::endl;
    }

  unsigned int nSlices = m_Volume->GetBufferedRegion().GetSize( 2 );
  std::vector<double> metric( nSlices );
  
  typedef SliceImageConstIterator<VolumeType> SliceIterator;
  typedef typename SliceIterator::OutputImageType ImageType;
  typedef SliceImageConstIterator<MaskVolumeType> MaskSliceIterator;

  SliceIterator vIterator( m_Volume );
  MaskSliceIterator mIterator( m_Mask );
  SliceIterator ravIterator( m_ResampledAtlasVolume );
  MaskSliceIterator rahmIterator( m_ResampledAtlasHeadMask );

  vIterator.GoToBegin();
  mIterator.GoToBegin();
  ravIterator.GoToBegin();
  rahmIterator.GoToBegin();

  for ( unsigned int j = 0; j < nSlices; j++ )
    {

    double m;
    ComputeCorrelation<ImageType,MaskImageType>( ravIterator.GetOutput(), vIterator.GetOutput(),
                                   rahmIterator.GetOutput(), mIterator.GetOutput(), m );
    metric[j] = m;

    if ( m_Verbose )
      {
      std::cout << j << ": " << m << std::endl;
      }

    ++vIterator;
    ++mIterator;
    ++ravIterator;
    ++rahmIterator;

    }


  // delegate upload to the utility
  m_Utility->SetSubImageMetric( metric );

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}


/**
 * Write out resampled atlas data to output directory
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
WriteResampledAtlas(
const char * outputDirectory,
bool useCount
)
{

  if ( m_Verbose )
    {
    std::cout << "WriteResampledAtlas ...";
    }

  std::string fname;

  std::ostringstream ostr;
  ostr << m_TransformModificationCount;

  // write out volume
  fname = outputDirectory;
  fname += "/atlasVolume";
  if ( useCount )
    {
    fname += "-";
    fname += ostr.str();
    }
  fname += ".mhd";

  WriteImage<VolumeType>( fname.c_str(), m_ResampledAtlasVolume );

  // write out mask
  fname = outputDirectory;
  fname += "/atlasMask";
  if ( useCount )
    {
    fname += "-";
    fname += ostr.str();
    }
  fname += ".mhd";

  WriteImage<MaskVolumeType>( fname.c_str(), m_ResampledAtlasMask );

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}

/**
 * Write out resampled volume data to output directory
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
WriteResampledVolume(
const char * outputDirectory,
bool useCount
)
{

  if ( m_Verbose )
    {
    std::cout << "WriteResampledVolume ...";
    }

  std::string fname;

  std::ostringstream ostr;
  ostr << m_TransformModificationCount;

  // write out volume
  fname = outputDirectory;
  fname += "/volume";
  if ( useCount )
    {
    fname += "-";
    fname += ostr.str();
    }
  fname += ".mhd";

  WriteImage<VolumeType>( fname.c_str(), m_ResampledVolume );

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}
/**
 * Extract view slice from volume in atlas space
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ExtractViewSlice(
unsigned int dim,
typename ImageType::Pointer & slice
)
{
  this->ExtractViewSlice( m_ResampledVolume, dim, slice );
}

template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ExtractViewSlice(
typename VolumeType::Pointer & volume,
unsigned int dim,
typename ImageType::Pointer & slice
)
{

  if ( m_Verbose )
    {
    std::cout << "ExtractViewSlice ...";
    }

  if ( m_Fiducials.IsNull() ||
       !m_Fiducials->PointLabelExists( "SliceView" ) )
    {
    if ( m_Verbose )
      {
      std::cout << "[failed]" << std::endl << std::endl;
      }
    return;
    }

  typename VolumeType::PointType point = m_Fiducials->GetPoint( "SliceView" );
  typename VolumeType::IndexType index; 
  m_AtlasVolume->TransformPhysicalPointToIndex( point, index );

  ExtractSlice<VolumeType,ImageType>( volume, slice, index[dim], dim );

  if ( m_Verbose )
    {
    std::cout << " [done]" << std::endl << std::endl;
    }

}


/**
 * Write transforms to output xml file
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
WriteTransformsToXML(
const char * xmlFile )
{

  // open the doc
  TiXmlDocument doc;
  TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "UTF-8", "" );
  doc.LinkEndChild( decl );

  // set up root (alignment) node
  TiXmlElement * alignment = new TiXmlElement( "alignment" );
  doc.LinkEndChild( alignment );

  // set up alignment3ds node
  TiXmlElement * alignment3ds = new TiXmlElement( "alignment3ds" );
  alignment->LinkEndChild( alignment3ds );
  ImageSeriesUtilitiesType::WriteTransformParametersToXML( m_SeriesVector, m_VolumeToAtlasTransform, m_AtlasToGridTransform,
                                                       m_Metric, alignment3ds );

  // set up alignment2ds node
  TiXmlElement * alignment2ds = new TiXmlElement( "alignment2ds" );
  alignment->LinkEndChild( alignment2ds );
  ImageSeriesUtilitiesType::WriteTransformParametersToXML( m_SeriesVector, alignment2ds );

  // close file
  doc.SaveFile( xmlFile );

}


/**
 * Load image series from xml file
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
LoadTransformsFromXML( const char * xmlFile )
{

  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile << std::endl );
    }
    
  TiXmlNode * node = doc.FirstChild( "alignment" );
  if ( !node )
    {
    itkExceptionMacro( << "Could not find alignment node in file " << xmlFile << std::endl );
    }

  TiXmlNode * alignment3ds = node->FirstChild( "alignment3ds" );
  if ( !alignment3ds )
    {
    itkExceptionMacro( << "Could not find alignment3ds node in file " << xmlFile << std::endl );
    }

  Affine3DTransformType::Pointer transform = Affine3DTransformType::New();
  m_Utility->ReadTransformParametersFromXML( alignment3ds, transform );
  this->SetVolumeToAtlasTransform( transform );

  TiXmlNode * alignment2ds = node->FirstChild( "alignment2ds" );
  if ( !alignment2ds )
    {
    itkExceptionMacro( << "Could not find alignment3ds node in file " << xmlFile << std::endl );
    }

  m_Utility->ReadTransformParametersFromXML( alignment2ds );

}

/**
 * Load image series from xml file
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
LoadEmbeddedTransformsFromXML( const char * xmlFile )
{

  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile << std::endl );
    }
    
   TiXmlNode * node;
   if ( ( node = doc.FirstChild( "specimen-image" ) ) || ( node = doc.FirstChild( "image-series" ) ) )
    {

    // read in 3D transform    
    TiXmlNode * meta_data = node->FirstChild( "meta-data" );
    if ( !meta_data )
      {
      itkExceptionMacro( << "Could not find meta-data node in file " << xmlFile << std::endl );
      }
      
    Affine3DTransformType::Pointer transform = Affine3DTransformType::New();
    m_Utility->ReadTransformParametersFromXML( meta_data, transform );
    this->SetVolumeToAtlasTransform( transform );
    
    // std::cout << transform << std::endl;
    
    // read in 2D transform
    TiXmlNode * sub_images = node->FirstChild( "sub-images" );
    if ( !sub_images )
      {
      itkExceptionMacro( << "Could not find sub-images node in file " << xmlFile << std::endl );     
      }
      
    TiXmlNode * sub_image = NULL;
    
    while ( ( sub_image = sub_images->IterateChildren( "sub-image", sub_image ) ) )
      {
      //std::cout << "sub image " << std::endl;
      
      TiXmlNode * alignment2d = sub_image->FirstChild( "alignment2d" );
      if ( !alignment2d )
        {
        m_Utility->ReadSubImageTransformParametersFromXML( alignment2d );
        }
      }
   
    }
   else
    {
    itkExceptionMacro( << "this file content is not currently supported" << std::endl );
    }

}

/**
 * Load color scales from xml
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
LoadColorScalesFromXML( const char * xmlFile )
{

  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile << std::endl );
    }
    
  TiXmlNode * node = doc.FirstChild( "scales" );
  if ( !node )
    {
    itkExceptionMacro( << "Could not find scales node in file " << xmlFile << std::endl );
    }

  TiXmlNode *n;
  TiXmlElement *e;
  n = NULL;

  while( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {

      std::string parameter = e->Value();
      if ( parameter.compare( "red" ) == 0 )
        {
        this->SetRedScale( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "green" ) == 0 )
        {
        this->SetGreenScale( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "blue" ) == 0 )
        {
        this->SetBlueScale( atof( e->GetText() ) );
        }

      } // end if (e)

    } // end while

}


/**
 * Compose with atlas-to-grid transform
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ComposeWithAtlasToGridTransform()
{
  m_VolumeToAtlasTransform->Compose( m_AtlasToGridTransform, false );
  m_AtlasToVolumeTransform->SetCenter( m_VolumeToAtlasTransform->GetCenter() );
  m_VolumeToAtlasTransform->GetInverse( m_AtlasToVolumeTransform );   
}

/**
 * Compose with grid to atlas transform
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ComposeWithGridToAtlasTransform()
{
  m_VolumeToAtlasTransform->Compose( m_GridToAtlasTransform, false );
  m_AtlasToVolumeTransform->SetCenter( m_VolumeToAtlasTransform->GetCenter() );
  m_VolumeToAtlasTransform->GetInverse( m_AtlasToVolumeTransform );   
}

/**
 * Compose with z-axis flip
 */
template<class TPixel>
void
RegistrationDispatcher<TPixel>::
ComposeWithZAxisFlip()
{
  Affine3DTransformType::Pointer a3t = Affine3DTransformType::New();
  ReferenceSpaceUtilities::PointType bc;
  bc.Fill( 0.0 );

  if ( !m_Fiducials.IsNull() && m_Fiducials->PointLabelExists( "BrainCenter" ) )
    {
    bc = m_Fiducials->GetPoint( "BrainCenter" );
    }

  a3t->SetIdentity();
  a3t->SetCenter( bc );
  Affine3DTransformType::ParametersType p = a3t->GetParameters();
  p[8] = -1.0;
  a3t->SetParameters( p );
  std::cout << a3t->GetParameters() << std::endl;

  m_VolumeToAtlasTransform->Compose( a3t, true );
  m_AtlasToVolumeTransform->SetCenter( m_VolumeToAtlasTransform->GetCenter() );
  m_VolumeToAtlasTransform->GetInverse( m_AtlasToVolumeTransform ); 

  m_TransformModificationCount++;

}


} // end namespace idp
} //end namespace itk

#endif
