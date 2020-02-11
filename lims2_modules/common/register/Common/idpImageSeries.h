/*=========================================================================

  idpImageSeries.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpImageSeries_h
#define __idpImageSeries_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkSize.h"
#include "itkIndex.h"
#include "itkImageRegion.h"
#include <string>
#include <algorithm>
#include <map>

#include "tinyxml.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkPolyLineParametricPath.h"

namespace itk
{
namespace idp
{

/** \class ReferenceSpace
 *  \brief Encapsulate lims2 reference_space objet
 */
class ReferenceSpace : public Object
{
public:
  /** Standard typedefs. */
  typedef ReferenceSpace      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(ReferenceSpace, Object);

  /** Set/Get the id. */
  itkSetMacro( Id, unsigned long );
  itkGetConstMacro( Id, unsigned long );

  /** Set/Get the plane of section. */
  itkSetStringMacro( Age );
  itkGetStringMacro( Age );

  /** Set/Get the specimen name. */
  itkSetStringMacro( Name );
  itkGetStringMacro( Name );

  /** Load from XML node. */
  void LoadFromXML( TiXmlNode * node );

protected:
  ReferenceSpace();
  ~ReferenceSpace();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  ReferenceSpace(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned long         m_Id;
  std::string           m_Name;
  std::string           m_Age;

};

/** \class SubImage
 *  \brief Encapsulate lims2 sub_image object
 *
 */
class SubImage : public Object
{
public:
  /** Standard typedefs. */
  typedef SubImage      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(SubImage, Object);

  /** Other types. */
  typedef itk::Size<2>  SizeType;
  typedef itk::Index<2> IndexType;
  typedef itk::ImageRegion<2> RegionType;
  typedef itk::CenteredEuler3DTransform<double> TransformType;
  typedef TransformType::Pointer         TransformPointer;
  typedef itk::CenteredRigid2DTransform<double> Transform2DType;
  typedef Transform2DType::Pointer       Transform2DPointer;
  typedef PolyLineParametricPath<2> PolyLineType;
  typedef PolyLineType::ContinuousIndexType VertexType;

  typedef std::vector<std::string> StringArray;
  typedef std::vector<RegionType> RegionArray;
  typedef std::vector<PolyLineType::VertexListType::Pointer> PolyLineArray;
  
  typedef std::map< std::string, PolyLineArray > GraphicLayerMap;

  /** Set/Get the specimen tissue index. */
  itkSetMacro( SpecimenTissueIndex, long );
  itkGetConstMacro( SpecimenTissueIndex, long );

  /** Set/Get the id. */
  itkSetMacro( Id, unsigned long );
  itkGetConstMacro( Id, unsigned long );

  /** Set/Get the bounding box corner. */
  itkSetMacro( Index, IndexType );
  itkGetConstMacro( Index, IndexType );

  /** Set/Get the bounding box size. */
  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );

  /** Set/Get image resolution. */
  itkSetMacro( ImageResolution, double );
  itkGetConstMacro( ImageResolution, double );

  /** Set/Get the jp2 filename. */
  itkSetStringMacro( Jp2Filename );
  itkGetStringMacro( Jp2Filename );

  /** Set/Get the expression jp2 filename. */
  itkSetStringMacro( ExpressionJp2Filename );
  itkGetStringMacro( ExpressionJp2Filename );

  /** Set/Get the projection mask jp2 filename. */
  itkSetStringMacro( ProjectionMaskJp2Filename );
  itkGetStringMacro( ProjectionMaskJp2Filename );

  /** Set/Get the downsample filename. */
  itkSetStringMacro( DownsampleFilename );
  itkGetStringMacro( DownsampleFilename );

  /** Set/Get the mask filename. */
  itkSetStringMacro( MaskFilename );
  itkGetStringMacro( MaskFilename );

    /** Set/Get the csv filename. */
  itkSetStringMacro( CSVFilename );
  itkGetStringMacro( CSVFilename );
  
  /** Set/Get the downsample filename. */
  itkSetStringMacro( TifFilename );
  itkGetStringMacro( TifFilename );  

  /** Set/Get the failed flag. */
  itkSetMacro( Failed, bool );
  itkGetConstMacro( Failed, bool );

  /** Load from XML node. */
  void LoadFromXML( TiXmlNode * node );

  /** Load additional information from XML node. */
  void LoadAdditionalDataFromXML( TiXmlNode * node );

  /** Set/Get Transform. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetConstObjectMacro( Transform, TransformType );
  
  /** Set/Get 2D Transform. */
  itkSetObjectMacro( Transform2D, Transform2DType );
  itkGetConstObjectMacro( Transform2D, Transform2DType );

  /** Set/Get tissue size in squared mm. */
  itkSetMacro( TissueSize, double );
  itkGetConstMacro( TissueSize, double );

  /** Set/Get the goodness of fit metric. */
  itkSetMacro( Metric, double );
  itkGetConstMacro( Metric, double );

  /** Get the exclusion regions */
  const RegionArray& GetExclusionRegions() const
    { 
    return m_ExclusionRegions; 
    }

  /** Get the graphic layers */
  const GraphicLayerMap& GetGraphicLayers() const
    {
    return m_GraphicLayers;
    }

  /** Get the number of zoom levels in the sub image */
  itkSetMacro( TierCount, unsigned int );
  itkGetConstMacro( TierCount, unsigned int );
  
   /** Set/Get the closest Nissl. */
  itkSetObjectMacro( ClosestNissl, Self );
  itkGetConstObjectMacro( ClosestNissl, Self );
  itkGetObjectMacro( ClosestNissl, Self );
  
protected:
  SubImage();
  ~SubImage();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void ParsePolygonsXML(TiXmlNode *n, PolyLineArray& polygonArray);
  void ParseGraphicLayers( TiXmlNode * node );

private:
  SubImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned long         m_Id;
  long                  m_SpecimenTissueIndex;
  IndexType             m_Index;
  SizeType              m_Size;
  unsigned long         m_SlideImageId;
  std::string           m_StorageDirectory;
  std::string           m_Jp2Filename;
  std::string           m_ExpressionJp2Filename;
  std::string           m_ProjectionMaskJp2Filename;
  std::string           m_DownsampleFilename;
  std::string           m_MaskFilename;
  std::string           m_CSVFilename;
  std::string           m_TifFilename;
  bool                  m_Failed;
  TransformPointer      m_Transform;
  Transform2DPointer    m_Transform2D;
  double                m_TissueSize;
  StringArray           m_Treatments;
  double                m_Metric;
  RegionArray           m_ExclusionRegions;
  unsigned int          m_TierCount;
  double                m_ImageResolution;
  Pointer               m_ClosestNissl;
  GraphicLayerMap       m_GraphicLayers;

};

bool SubImageLessThan (
  SubImage::Pointer s1,
  SubImage::Pointer s2 );



/** \class ImageSeries
 *  \brief Encapsulate lims image_series object
 *
 */
class ImageSeries : public Object
{
public:
  /** Standard typedefs. */
  typedef ImageSeries      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageSeries, Object);

  /** Other types. */
  typedef itk::Size<2>  SizeType;
  typedef std::vector<SubImage::Pointer> SubImageArray;
  typedef std::vector<std::string> StringArray;

  /** Set/Get the image id. */
  itkSetMacro( Id, unsigned long );
  itkGetConstMacro( Id, unsigned long );

  /** Set/Get the plane of section. */
  itkSetStringMacro( PlaneOfSection );
  itkGetStringMacro( PlaneOfSection );

  /** Set/Get the plane of section. */
  itkSetStringMacro( Age );
  itkGetStringMacro( Age );

  /** Set/Get the specimen name. */
  itkSetStringMacro( SpecimenName );
  itkGetStringMacro( SpecimenName );

  /** Set/Get the intersection interval. */
  itkSetMacro( InterSectionInterval, unsigned long );
  itkGetConstMacro( InterSectionInterval, unsigned long );

  /** Set/Get image resolution. */
  itkSetMacro( ImageResolution, double );
  itkGetConstMacro( ImageResolution, double );

  /** Set/Get section thickness. */  
  itkSetMacro( SectionThickness, double );
  itkGetConstMacro( SectionThickness, double );

  /** Set/get reference space. */
  itkSetObjectMacro( ReferenceSpace, ReferenceSpace );
  itkGetConstObjectMacro( ReferenceSpace, ReferenceSpace );
  itkGetObjectMacro( ReferenceSpace, ReferenceSpace );

  /** Set/Get the storage directory. */
  itkSetStringMacro( StorageDirectory );
  itkGetStringMacro( StorageDirectory );

  /** Set/Get the workflow state. */
  itkSetStringMacro( WorkflowState );
  itkGetStringMacro( WorkflowState );
  
  /** Set/Get the hemisphere. */
  itkSetStringMacro( Hemisphere );
  itkGetStringMacro( Hemisphere );  


  /** Load from XML file. */
  void LoadFromXML( const char * xmlFile );
  void LoadFromXML( TiXmlNode * node );
  void LoadMetaDataFromXML( TiXmlNode * node );
  void LoadSubImagesFromXML( TiXmlNode * node );
  void LoadTreatmentsFromXML( TiXmlNode * node );

  /** Sort subimages. */
  void SortSubImages();

  /** Get the reference to the vector of sub images. */
  const SubImageArray & GetSubImages() const
    {
    return m_SubImages;
    }

  /** Insert sub image into sub image array */
  void InsertSubImage( SubImage * si );

  /** Get the reference to the string array of treatments. */
  const StringArray & GetTreatments() const
    {
    return m_Treatments;
    }

  /** Insert a treatment */
  void InsertTreatment( const char * treatment );

  /** Reset treatment */
  void ResetTreatments()
    {
    m_Treatments.resize(0);
    }

  /** Has specify treatment? */
  bool HasTreatment( const char * treatment )
    {
    StringArray::iterator it;
    it = std::find( m_Treatments.begin(), m_Treatments.end(), treatment );
    return ( it != m_Treatments.end() );
    }

  /** Copy information. */
  void CopyInformation( const Pointer & other );
    
protected:
  ImageSeries();
  ~ImageSeries();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  ImageSeries(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:
  unsigned long                   m_Id;
  std::string                     m_SpecimenName;
  std::string                     m_PlaneOfSection;
  std::string                     m_Age;
  SubImageArray                   m_SubImages;
  unsigned long                   m_InterSectionInterval;
  StringArray                     m_Treatments;
  std::string                     m_WorkflowState;
  std::string                     m_StorageDirectory;
  double                          m_ImageResolution;
  double                          m_SectionThickness;
  ReferenceSpace::Pointer         m_ReferenceSpace;
  std::string                     m_Hemisphere;


};


/** \class Specimen
 *  \brief Encapsulate lims specimen object
 *
 */
class Specimen : public Object
{
public:
  /** Standard typedefs. */
  typedef Specimen      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(Specimen, Object);

  /** Other types. */
  typedef std::vector<ImageSeries::Pointer> ImageSeriesArray;

  /** Set/Get the image id. */
  itkSetMacro( Id, unsigned long );
  itkGetConstMacro( Id, unsigned long );

  /** Set/Get the plane of section. */
  itkSetStringMacro( PlaneOfSection );
  itkGetStringMacro( PlaneOfSection );

  /** Set/Get the plane of section. */
  itkSetStringMacro( Age );
  itkGetStringMacro( Age );

  /** Set/Get the specimen name. */
  itkSetStringMacro( Name );
  itkGetStringMacro( Name );

  /** Set/Get the storage directory. */
  itkSetStringMacro( StorageDirectory );
  itkGetStringMacro( StorageDirectory );


  /** Load from XML file. */
  void LoadFromXML( TiXmlNode * node );
  void LoadMetaDataFromXML( TiXmlNode * node );
  void LoadImageSeriesFromXML( TiXmlNode * node );

  /** Get the reference to the vector of image series. */
  const ImageSeriesArray & GetImageSeries() const
    {
    return m_ImageSeries;
    }

    
protected:
  Specimen();
  ~Specimen();
  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  Specimen(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:
  unsigned long                   m_Id;
  std::string                     m_Name;
  std::string                     m_PlaneOfSection;
  std::string                     m_Age;
  std::string                     m_StorageDirectory;
  ImageSeriesArray                m_ImageSeries;

};


} // end namespace idp
} //end namespace itk

#endif
