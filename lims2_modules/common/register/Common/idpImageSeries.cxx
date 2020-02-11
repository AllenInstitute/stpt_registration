/*=========================================================================

  idpImageSeries.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpImageSeries_cxx
#define __idpImageSeries_cxx

#include "idpImageSeries.h"
#include "itkNumericTraits.h"
#include "tinyxml.h"
#include <algorithm>

#include <itksys/SystemTools.hxx>
#include "idpRegistrationUtilities.h"

namespace itk
{
namespace idp
{

/**
 * Constructor
 */
ReferenceSpace::ReferenceSpace()
{
  m_Id = 0;
  m_Name = "";
  m_Age = "";
}

/**
 * Destructor
 */
ReferenceSpace::~ReferenceSpace()
{

}

/**
 * PrintSelf
 */
void 
ReferenceSpace::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Id: " << this->GetId() << std::endl;
  os << indent << "Name: " << this->GetName() << std::endl;
  os << indent << "Age: " << this->GetAge() << std::endl;
}

/** 
 * Load data from XML node
 */
void
ReferenceSpace::LoadFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();

      if ( parameter.compare( "age" ) == 0 )
        {
        this->SetAge( e->GetText() );
        }
      else if ( parameter.compare("id") == 0 )
        {
        this->SetId( atoi( e->GetText() ) );
        }
      else if ( parameter.compare("name") == 0 )
        {
        this->SetName( e->GetText() );
        }


      } // end if (e)   

    } // end while

}



/**
 * Constructor
 */
SubImage::SubImage()
{
  m_Id = 0;
  m_SpecimenTissueIndex = NumericTraits<long>::min();
  m_Index.Fill( 0 );
  m_Size.Fill( 0 );
  m_ImageResolution  = 0.0;
  m_Jp2Filename = "";
  m_ExpressionJp2Filename = "";
  m_ProjectionMaskJp2Filename = "";
  m_DownsampleFilename = "";
  m_MaskFilename = "";
  m_CSVFilename = "";
  m_TifFilename = "";
  m_Failed = false;
  m_Treatments.resize(0);

  m_Transform = TransformType::New();
  TransformType::ParametersType p( m_Transform->GetNumberOfParameters() );
  p.Fill( 0.0 );
  m_Transform->SetParameters( p );
  
  m_Transform2D = Transform2DType::New();
  Transform2DType::ParametersType q( m_Transform2D->GetNumberOfParameters() );
  q.Fill( 0.0 );
  m_Transform2D->SetParameters( q );

  m_TissueSize = 0.0;
  m_Metric = 0.0;
  m_TierCount = 1;
  
  m_ClosestNissl = 0;
  
  m_GraphicLayers.clear();
  
}

/**
 * Destructor
 */
SubImage::~SubImage()
{

}

/**
 * PrintSelf
 */
void 
SubImage::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Id: " << this->GetId() << std::endl;
  os << indent << "SpecimenTissueIndex: " << this->GetSpecimenTissueIndex() << std::endl;
  os << indent << "Index: " << this->GetIndex() << std::endl;
  os << indent << "Size: " << this->GetSize() << std::endl;
  os << indent << "ImageResolution: " << this->GetImageResolution() << std::endl;
  os << indent << "Jp2Filename: " << this->GetJp2Filename() << std::endl;
  os << indent << "ExpressionJp2Filename: " << this->GetExpressionJp2Filename() << std::endl;
  os << indent << "ProjectionMaskJp2Filename: " << this->GetProjectionMaskJp2Filename() << std::endl;
  os << indent << "DownsampleFilename: " << this->GetDownsampleFilename() << std::endl;
  os << indent << "MaskFilename: " << this->GetMaskFilename() << std::endl;
  os << indent << "CSVFilename: " << this->GetCSVFilename() << std::endl;
  os << indent << "TifFilename: " << this->GetTifFilename() << std::endl;
  os << indent << "Failed: " << this->GetFailed() << std::endl;
  os << indent << "Transform: " << this->GetTransform() << std::endl;
  os << indent << "Transform2D: " << this->GetTransform2D() << std::endl;
  os << indent << "TissueSize: " << this->GetTissueSize() << std::endl;
  os << indent << "ClosestNissl: " << this->GetClosestNissl() << std::endl;
}

/** 
 * Load data from XML node
 */
void
SubImage::LoadFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;
  std::vector<unsigned int> zooms;
  std::vector<RegionType> exclusionRegions;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();
      
      //std::cout << parameter << std::endl;

      if ( parameter.compare( "specimen-tissue-index" ) == 0 )
        {
        this->SetSpecimenTissueIndex( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "id" ) == 0 )
        {
        this->SetId( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "x" ) == 0 )
        {
        m_Index[0] = atoi( e->GetText() );
        }
      else if ( parameter.compare( "y" ) == 0 )
        {
        m_Index[1] = atoi( e->GetText() );
        }
      else if ( parameter.compare( "width" ) == 0 )
        {
        m_Size[0] = atoi( e->GetText() );
        }
      else if ( parameter.compare( "height" ) == 0 )
        {
        m_Size[1] = atoi( e->GetText() );
        }
      else if ( parameter.compare( "resolution" ) == 0 )
        {
        m_ImageResolution = atof( e->GetText() );
        }
      else if ( parameter.compare( "primary-image-path" ) == 0 )
        {
        this->SetJp2Filename( e->GetText() );
        }
      else if ( parameter.compare( "primary-image" ) == 0 )
        {
        // new generic style primary-image path
        TiXmlNode *n = e->FirstChildElement("path");
        if ( n )
          { 
          TiXmlElement *p = n->ToElement();
          if ( p )
            {
            std::string ext = itksys::SystemTools::GetFilenameLastExtension( p->GetText() );
            if ( ext.compare(".jp2") == 0 )
              {
              this->SetJp2Filename( p->GetText() );
              }
            else
              {
              this->SetTifFilename( p->GetText() );
              }
            } // if p
          } // if n
        }
      else if ( parameter.compare( "mask-image-path" ) == 0 )
        {
        this->SetMaskFilename( e->GetText() );
        }
      else if ( parameter.compare( "expression-image-path" ) == 0 )
        {
        this->SetExpressionJp2Filename( e->GetText() );
        }
      else if ( parameter.compare( "projection-mask-image-path" ) == 0 )
        {
        this->SetProjectionMaskJp2Filename( e->GetText() );
        }
      else if ( parameter.compare( "csv-file-name" ) == 0 )
        {
        this->SetCSVFilename( e->GetText() );
        }
      else if ( parameter.compare( "tiff-image-path" ) == 0 )
        {
        this->SetTifFilename( e->GetText() );
        }
      else if ( parameter.compare( "failed" ) == 0 )
        {
        const char * txt = e->GetText();
        if ( txt != NULL )
          {
          std::string str( txt );
          this->SetFailed( str.compare("true") == 0 );
          }
        else
          {
          this->SetFailed( false );
          }
        }
      else if ( parameter.compare( "tier-count" ) == 0 )
        {
	      this->SetTierCount( ::atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "closest-nissl-sub-image" ) == 0 )
        {
        m_ClosestNissl = Self::New();
        m_ClosestNissl->LoadFromXML( n );
        }        
      else if ( parameter.compare( "exclusions" ) == 0 )
        {
	TiXmlNode *excsNode = NULL;
	    
	// annotations contains a list of tissue bubble/debris annotations
	while ( ( excsNode = n->IterateChildren( excsNode ) ) )
	  {
	  TiXmlNode *excNode = NULL;
	  TiXmlElement *excElement = NULL;
		
	  IndexType excOrigin; excOrigin.Fill( 0 );
	  SizeType excSize; excSize.Fill( 0 );
	  unsigned int zoom = 1;
	  
	  // each annotation has position and size parameters
	  while ( ( excNode = excsNode->IterateChildren( excNode ) ) )
	    {
	    excElement = excNode->ToElement();
		  
	    if (excElement)
	      {
	      std::string excParameter = excElement->Value();
		  
	      if ( excParameter.compare( "x" ) == 0 )
	        {
		excOrigin[0] = ::atoi( excElement->GetText() );
		}
	      else if ( excParameter.compare( "y" ) == 0 ) 
	        {
		excOrigin[1] = ::atoi( excElement->GetText() );
		}
	      else if ( excParameter.compare( "width" ) == 0 )
	        {
		excSize[0] = ::atoi( excElement->GetText() );
		}
	      else if ( excParameter.compare( "height" ) == 0 )
	        {
		excSize[1] = ::atoi( excElement->GetText() );
		}
	      else if ( excParameter.compare( "zoom" ) == 0 )
	        {
		zoom = ::atoi( excElement->GetText() );
		}
	      } // if (excElement)
	    } // while (excNode)
	      
	  RegionType excRegion;
	  excRegion.SetSize(excSize);
	  excRegion.SetIndex(excOrigin);
	  
	  exclusionRegions.push_back(excRegion);
	  zooms.push_back(zoom);
	  } // while (excsNode)
	} // end if (exclusions)
      else if ( parameter.compare( "missing-tile-polygons" ) == 0 )
        {
        PolyLineArray dummy;
        m_GraphicLayers[ "missing_tile" ] = dummy;
	ParsePolygonsXML(n, m_GraphicLayers[ "missing_tile" ]);
	}
      else if ( parameter.compare( "no-data-polygons" ) == 0 )
        {
        PolyLineArray dummy;
        m_GraphicLayers[ "no_data" ] = dummy;
	ParsePolygonsXML(n, m_GraphicLayers[ "no_data" ]);
	}
      else if ( parameter.compare( "aav-exclusion-polygons" ) == 0 )
        {
        PolyLineArray dummy;
        m_GraphicLayers[ "aav_exclusion" ] = dummy;
	ParsePolygonsXML(n, m_GraphicLayers[ "aav_exclusion" ]);
	}
  
      else if ( parameter.compare( "graphic-layers" ) == 0 )
        {
        this->ParseGraphicLayers( n );
        }
      } // end if (e)   
    } // end while

  // rescale the exclusion region sizes by the zoom level scale factor
  // this has to be done down here because m_TierCount may not be initialized before exclusions
  for (unsigned int i=0; i<zooms.size(); i++)
    {
    unsigned int zoomFactor = 1 << (this->GetTierCount() - zooms[i] - 1);
    RegionType::SizeType size = exclusionRegions[i].GetSize();
    size[0] *= zoomFactor;
    size[1] *= zoomFactor;
    exclusionRegions[i].SetSize(size);
    this->m_ExclusionRegions.push_back(exclusionRegions[i]);
    }

    // --- to do: move mask/downsample image logic at the end ------//
    std::string fname = this->GetJp2Filename();
    std::string mname = this->GetMaskFilename();

    if ( !fname.empty() && mname.empty() )
      {
      // no explicit mask file name
      // try appending "-mask.png" first else fall back to old style
      mname = fname;
      mname += "-mask.png";
      ChangePaths( mname );
      if ( !itksys::SystemTools::FileExists( mname.c_str(), true ) )
        {
        std::string otherName = this->GetJp2Filename();
        otherName = itksys::SystemTools::GetFilenamePath( otherName );
        otherName += "/mask.png";
        ChangePaths( otherName );
        if ( itksys::SystemTools::FileExists( otherName.c_str(), true ) )
          {
          mname = otherName;
          this->SetMaskFilename( mname.c_str() );
          }
        }
      else
        {
        this->SetMaskFilename( mname.c_str() );
        }
     }
}

void SubImage::ParsePolygonsXML(TiXmlNode *n, SubImage::PolyLineArray& polygons)
{
    polygons.clear();
    
    std::vector< std::vector<float> > polygonCoords;

    // format is:
    // <missing-tile-polygons>
    //   <missing-tile-polygon>
    //     <Polygon>
    //       <d>x,y,x,y,x,y,...</d>
    
    TiXmlNode *excsNode = NULL;
	    
    // iterate through the list of missing polygons (<exclusion-polygon>)
    while ( ( excsNode = n->IterateChildren( excsNode ) ) )
      {
      // <Polygon>
      TiXmlNode *polyNode = excsNode->FirstChild();

      // <d>
      if (polyNode) 
        {
	TiXmlNode *dNode = polyNode->FirstChildElement("d");
	    
	if (dNode) 
	  {
	  TiXmlElement *dElement = dNode->ToElement();

	  if (dElement)
	    {
	    std::vector<float> coords;

	    // format is "x,y,x,y,x,y,..."
	    std::string coordString = dElement->GetText();
	    std::istringstream iss(coordString);
	    std::string token;
	    
	    while (std::getline(iss, token, ','))
		coords.push_back(::atof(token.c_str()));
	    
	    polygonCoords.push_back(coords);
	    
	    } // if (dElement)
	  } // if (dNode)
	} // if (polyNode)
      } // while (excsNode)

    polygons.clear();
    for (unsigned int i = 0; i < polygonCoords.size(); i++)
      {
      std::vector<float>& polygon = polygonCoords[i];
      
      int ncoords = polygon.size();
    
      if (ncoords <= 2 || ncoords % 2 != 0)
        {
	std::cout << "skipping exclusion polygon " << i << ", incorrect number of vertex coordinates." << std::endl;
	continue;
	}

      PolyLineType::VertexListType::Pointer verts = PolyLineType::VertexListType::New();    
      for (unsigned int j = 0; j < polygon.size(); j+=2)
        {
	VertexType v;
	v[0] = polygon[j];
	v[1] = polygon[j+1];
//      p->AddVertex(v);   
	verts->push_back(v);
	}
    
      polygons.push_back(verts);
    }  
}

void SubImage::ParseGraphicLayers( TiXmlNode * node )
{
  //std::cout << "parsing graphic layers" << std::endl;
  
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode * n = NULL;

  while ( ( n = node->IterateChildren( "graphic-layer", n ) ) )
    {

    TiXmlNode * cn = NULL;
    TiXmlElement * ce = NULL;
    
    std::string group_label;    
    cn = n->FirstChildElement("group-label-name");

    if ( cn )
      {
      ce = cn->ToElement();
      if ( ce )
        {
        //std::cout << "load layer: " << ce->GetText() << std::endl;
        group_label = ce->GetText();
        }
      }
    
    if ( group_label == "" )
      {
      continue;
      }
    
    PolyLineArray dummy;
    m_GraphicLayers[ group_label ] = dummy;
    
    cn = NULL;
    cn = n->FirstChildElement("graphic-objects");
    if ( cn )
      {
      this->ParsePolygonsXML( cn, m_GraphicLayers[group_label] );
      }
    
    }
}

/** Operator to comparing the specimen index of the sub images
 *
 */
bool SubImageLessThan(
  SubImage::Pointer s1,
  SubImage::Pointer s2 )
{
  return s1->GetSpecimenTissueIndex() < s2->GetSpecimenTissueIndex();
}


/**
 * Constructor
 */
ImageSeries::ImageSeries()
{
  m_Id = 0;
  m_SpecimenName = "";
  m_PlaneOfSection = "";
  m_Age = "";
  m_InterSectionInterval = 0;
  m_ImageResolution = 0.0;
  m_SectionThickness = 0.0;
  m_SubImages.resize( 0 );
  m_Treatments.resize( 0 );
  m_WorkflowState = "";
  m_StorageDirectory = "";
  m_ReferenceSpace = 0;
  m_Hemisphere = "";
}

/**
 * Destructor
 */
ImageSeries::~ImageSeries()
{

}

/**
 * PrintSelf
 */
void 
ImageSeries::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Id: " << this->GetId() << std::endl;
  os << indent << "SpecimenName: " << this->GetSpecimenName() << std::endl;
  os << indent << "PlaneOfSection: " << this->GetPlaneOfSection() << std::endl;
  os << indent << "Age: " << this->GetAge() << std::endl;
  os << indent << "InterSectionInterval: " << this->GetInterSectionInterval() << std::endl;
  os << indent << "ImageResolution: " << this->GetImageResolution() << std::endl;
  os << indent << "SectionThickness: " << this->GetSectionThickness() << std::endl;
  os << indent << "WorkflowState: " << this->GetWorkflowState() << std::endl;
  os << indent << "StorageDirectory: " << this->GetStorageDirectory() << std::endl;
  os << indent << "Hemisphere: " << this->GetHemisphere() << std::endl;

  os << indent << "NumberOfTreatments: " << m_Treatments.size() << std::endl;
  for ( unsigned int k = 0; k < m_Treatments.size(); k++ )
    {
    os << indent.GetNextIndent() << "Name: " << m_Treatments[k] << std::endl;
    }

  os << indent << "NumberOfSubImages: " << m_SubImages.size() << std::endl;
  for ( unsigned int k = 0; k < m_SubImages.size(); k++ )
    {
  //  m_SubImages[k]->Print( os, indent.GetNextIndent() );
    }

}

/** 
 * Load data from XML file
 */
void
ImageSeries::LoadFromXML( const char * xmlFile )
{
  // Open xml file
  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    itkExceptionMacro( << "Could not load xml file " << xmlFile );
    }

  TiXmlNode * imageSeriesNode;

  // Locate the root node
  imageSeriesNode = doc.FirstChild( "image-series" );
  if ( !imageSeriesNode )
    {
    itkExceptionMacro( << "Can not find image-series node in file " << xmlFile );
    }

  this->LoadFromXML( imageSeriesNode );

}

/** 
 * Load data from XML node
 */
void
ImageSeries::LoadFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();
      if ( parameter.compare( "meta-data" ) == 0 )
        {
        this->LoadMetaDataFromXML( n );
        }
      else if ( parameter.compare("sub-images") == 0 )
        {
        this->LoadSubImagesFromXML( n );
        }

      } // end if (e)   

    } // end while

}

/** 
 * Load meta data from XML node
 */
void
ImageSeries::LoadMetaDataFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();
      if ( parameter.compare( "id" ) == 0 )
        {
        this->SetId( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "plane-of-section" ) == 0 )
        {
        this->SetPlaneOfSection( e->GetText() );
        }
      else if ( parameter.compare( "age" ) == 0 )
        {
        this->SetAge( e->GetText() );
        }
      else if ( parameter.compare( "treatments" ) == 0 )
        {
        this->LoadTreatmentsFromXML( n );
        }
      else if ( parameter.compare( "workflow-state" ) == 0 )
        {
        this->SetWorkflowState( e->GetText() );
        }
      else if ( parameter.compare( "storage-directory" ) == 0 )
        {
        this->SetStorageDirectory( e->GetText() );
        }
      else if ( parameter.compare( "resolution" ) == 0 )
        {
        this->SetImageResolution( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "section-thickness" ) == 0 )
        {
        this->SetSectionThickness( atof( e->GetText() ) );
        }
      else if ( parameter.compare("reference-space") == 0 )
        {
        ReferenceSpace::Pointer rs = ReferenceSpace::New();
        rs->LoadFromXML( n );
        this->SetReferenceSpace( rs );
        }
      else if ( parameter.compare( "hemisphere" ) == 0 )
        {
        this->SetHemisphere( e->GetText() );
        }        

      } // end if (e)   

    } // end while

}

/** 
 * Load sub images data from XML node
 */
void
ImageSeries::LoadSubImagesFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode * n = NULL;

  while ( ( n = node->IterateChildren( "sub-image", n ) ) )
    {
    SubImage::Pointer subimg = SubImage::New();
    subimg->LoadFromXML( n );
    m_SubImages.push_back( subimg );
    }

  this->SortSubImages();

}

/** 
 * Sort subimages
 */
void
ImageSeries::SortSubImages()
{
    
  // sanity check for number of images
  if ( m_SubImages.size() == 0 )
    {
    return;
    }
    
  // sort the sub images by the specimen tissue index
  sort( m_SubImages.begin(), m_SubImages.end(), SubImageLessThan );

  // compute the inter-section interval
  signed int interval = 1000;
  for ( unsigned int k = 1; k < m_SubImages.size(); k++ )
    {
    signed int diff = m_SubImages[k]->GetSpecimenTissueIndex() -
                      m_SubImages[k-1]->GetSpecimenTissueIndex();

    if ( diff == 0 )
      {
      // Need to handle duplicate 
      }
    else if ( diff < interval )
      {
      interval = diff;
      }
    }
    
  if ( m_SubImages.size() == 1 )
    {
    interval = 1;
    }    

  // sanity check to see if check if this interval make sense in edge cases
  int minIndex = m_SubImages[0]->GetSpecimenTissueIndex();
  bool intervalOk = true;
  for ( unsigned int k = 1; k < m_SubImages.size(); k++ )
    {
    int index = m_SubImages[k]->GetSpecimenTissueIndex() - minIndex;
    if ( index % interval )
      {
      intervalOk = false;
      //std::cout << m_SubImages[k] << std::endl;
      break;
      }
    }
  if (!intervalOk)
    {
    interval = 1;
    }

  this->SetInterSectionInterval( static_cast<unsigned long>( interval ) );

}


/**
 * Insert one sub image into sub image array
 */
void
ImageSeries::InsertSubImage( SubImage * si )
{
  m_SubImages.push_back( si );
}

/**
 * Copy information from another series
 */
void
ImageSeries::CopyInformation( const Pointer & other )
{
  this->SetId( other->GetId() );
  this->SetSpecimenName( other->GetSpecimenName() );
  this->SetPlaneOfSection( other->GetPlaneOfSection() );
  this->SetAge( other->GetAge() );
  this->SetStorageDirectory( other->GetStorageDirectory() );
  this->SetWorkflowState( other->GetWorkflowState() );
  this->SetImageResolution( other->GetImageResolution() );
  this->SetSectionThickness( other->GetSectionThickness() );
  this->SetReferenceSpace( other->GetReferenceSpace() );

  this->m_Treatments.resize(0);
  for( unsigned int j = 0; j < other->GetTreatments().size(); j++ )
    {
    this->InsertTreatment( other->GetTreatments()[j].c_str() );
    }
  
}


/**
 * Insert treatment into array
 */
void
ImageSeries::InsertTreatment( const char * treatment )
{
  m_Treatments.push_back( treatment );
}


/** 
 * Load treatment data from XML
 */
void
ImageSeries::LoadTreatmentsFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode * n = NULL;

  while ( ( n = node->IterateChildren( "treatment", n ) ) )
    {

    TiXmlNode * sn = NULL;
    TiXmlElement *e;

    while ( ( sn = n->IterateChildren( sn ) ) )
        {
        e = sn->ToElement();

        if ( e )
          {

          std::string parameter = e->Value();
          if ( parameter.compare( "name" ) == 0 )
            {
            this->InsertTreatment( e->GetText() );
            }

          } // if e

        } // while sn

    } // while n

}

/**
 * Constructor
 */
Specimen::Specimen()
{
  m_Id = 0;
  m_Name = "";
  m_PlaneOfSection = "";
  m_Age = "";
  m_StorageDirectory = "";
  m_ImageSeries.resize( 0 );
}

/**
 * Destructor
 */
Specimen::~Specimen()
{

}

/**
 * PrintSelf
 */
void 
Specimen::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Id: " << this->GetId() << std::endl;
  os << indent << "Name: " << this->GetName() << std::endl;
  os << indent << "PlaneOfSection: " << this->GetPlaneOfSection() << std::endl;
  os << indent << "Age: " << this->GetAge() << std::endl;
  os << indent << "StorageDirectory : " << this->GetStorageDirectory() << std::endl;
  os << indent << "NumberOfImageSeries: " << m_ImageSeries.size() << std::endl;
  for ( unsigned int k = 0; k < m_ImageSeries.size(); k++ )
    {
   // m_ImageSeries[k]->Print( os, indent.GetNextIndent() );
    }
}


/** 
 * Load data from XML node
 */
void
Specimen::LoadFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();

      if ( parameter.compare( "meta-data" ) == 0 )
        {
        this->LoadMetaDataFromXML( n );
        }
      else if ( parameter.compare("image-series") == 0 )
        {
        this->LoadImageSeriesFromXML( n );
        }

      } // end if (e)   

    } // end while

}


/** 
 * Load meta data from XML node
 */
void
Specimen::LoadMetaDataFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode *n;
  TiXmlElement *e;

  n = NULL;
  while ( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {
      std::string parameter = e->Value();

      if ( parameter.compare( "id" ) == 0 )
        {
        this->SetId( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "plane-of-section" ) == 0 )
        {
        this->SetPlaneOfSection( e->GetText() );
        }
      else if ( parameter.compare( "age" ) == 0 )
        {
        this->SetAge( e->GetText() );
        }
      else if ( parameter.compare("name") == 0 )
        {
        this->SetName( e->GetText() );
        }
      else if ( parameter.compare("storage-directory") == 0 )
        {
        this->SetStorageDirectory( e->GetText() );
        }


      } // end if (e)   

    } // end while

}

/** 
 * Load image series data from XML node
 */
void
Specimen::LoadImageSeriesFromXML( TiXmlNode * node )
{
  if ( !node )
    {
    itkExceptionMacro( << "TiXmlNode is NULL" );
    }

  TiXmlNode * n = NULL;

  while ( ( n = node->IterateChildren( "image-series", n ) ) )
    {
    ImageSeries::Pointer is = ImageSeries::New();
    is->LoadFromXML( n );
    m_ImageSeries.push_back( is );
    }

}


} // end namespace idp
} //end namespace itk

#endif
