/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $File: //depot/bioapps/infoapps/DevMouse/Code/main/Utilities/ValidatePolygonMapping.cxx $
  Language:  C++
  Date:      $DateTime: 2010/05/07 16:35:34 $
  Version:   $Revision: #1 $

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"
#include "itkVector.h"
#include <string>
#include <map>
#include "tinyxml.h"
#include "itkImageRegionIteratorWithIndex.h"


// --- define parameter based classes and read functions ---//
class Parameters
{
public:
  std::string BrainName;
  std::string SlabType;
  unsigned int SlabDirection;
  itk::Vector<double,3> Offset;
  unsigned long IdOffset;

Parameters()
  {
  BrainName = "";
  SlabType = "";
  SlabDirection = 0;
  Offset.Fill( 0.0 );
  IdOffset = 0;
  }
};

int FetchParametersFromXMLFile( const char *xmlFile, Parameters & p );


int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "parameterXMLFile chunkIdVolume ";
    std::cout << "inputCSVFile outputCSVFile outputVolume ";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<unsigned char, 3> VolumeType;
  typedef itk::Image<float, 3> FloatVolumeType;
  typedef itk::Image<unsigned short, 3> ShortVolumeType;

  typedef VolumeType::RegionType RegionType;
  typedef VolumeType::PointType PointType;
  typedef VolumeType::IndexType IndexType;

  typedef itk::ImageRegionIteratorWithIndex<VolumeType> Iterator;
  typedef itk::ImageRegionIteratorWithIndex<FloatVolumeType> FloatIterator;
  typedef itk::ImageRegionIteratorWithIndex<ShortVolumeType> ShortIterator;

  typedef itk::AffineTransform<double,3> TransformType;

  std::string parameterXMLFile = argv[1];
  std::string chunkIdVolume = argv[2];
  std::string inputCSVFile  = argv[3];
  std::string outputCSVFile  = argv[4];
  std::string outputVolume = argv[5];

  Parameters p;
  ShortVolumeType::Pointer cidVolume;
  VolumeType::Pointer mappedPositions;

  try
    {
    int err;

    err = FetchParametersFromXMLFile( parameterXMLFile.c_str(), p );
    if ( err ) { return EXIT_FAILURE; }

    std::cout << p.Offset << " " << p.IdOffset << std::endl;

    itk::idp::ReadImage<ShortVolumeType>( chunkIdVolume.c_str(), cidVolume );
    itk::idp::MakeImage<VolumeType,ShortVolumeType>( cidVolume, mappedPositions, 0 );

    // open input and output files
    std::ifstream input( inputCSVFile.c_str(), std::ios::in );
    if ( !input.fail() )
      {
      std::ofstream output( outputCSVFile.c_str(), std::ios::out );
      if ( !output.fail() )
        {

        // loop input file
        std::string buffer;
        getline( input, buffer );

        output << buffer;
        output << ",mapped_chunk_id,compare" << std::endl;

        while( getline( input, buffer ) )
          {
          std::istringstream instream;
          instream.str( buffer );
          std::string item;
          std::vector<std::string> components;
        
          // Parse comment separated text
          while( getline( instream, item, ',') )
            {
            components.push_back( item );
            }

          if( components.size() < 14 )
            {
            continue;
            }

          output << buffer;

          ShortVolumeType::IndexType index;
          index[0] = atoi(components[11].c_str());
          index[1] = atoi(components[12].c_str());
          index[2] = atoi(components[13].c_str());

          unsigned int cid = atoi(components[0].c_str());
          unsigned int value = cidVolume->GetPixel( index );
          value += p.IdOffset;

          mappedPositions->SetPixel( index, 255 );

          output << "," << value;
          output << "," << (cid == value) << std::endl;

          } // end while

        }// end if output

      output.close();

      }// end if input

    input.close();

    itk::idp::WriteImage<VolumeType>( outputVolume.c_str(), mappedPositions );

    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}


int FetchParametersFromXMLFile( const char *xmlFile, Parameters & p )
{

  // open xml file
  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
    std::cout << "Could not load xml file " << xmlFile;
    return -1;
    }

  TiXmlNode * node;

  // Locate the root node
  node = doc.FirstChild( "parameters" );
  if ( !node )
    {
    std::cout << "Can not find parameters node in file " << xmlFile;
    return -1;
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
      if ( parameter.compare( "brain-name" ) == 0 )
        {
        p.BrainName = e->GetText();
        }
      else if ( parameter.compare( "slab-type" ) == 0 )
        {
        p.SlabType = e->GetText();
        }
      else if ( parameter.compare( "slab-direction" ) == 0 )
        {
        p.SlabDirection = ( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "x-offset" ) == 0 )
        {
        p.Offset[0] = ( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "y-offset" ) == 0 )
        {
        p.Offset[1] = ( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "z-offset" ) == 0 )
        {
        p.Offset[2] = ( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "id-offset" ) == 0 )
        {
        p.IdOffset = ( atoi( e->GetText() ) );
        }

      } // end if (e)

    } // end while

  return 0;

}