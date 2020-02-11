/*=========================================================================

  Create2DMappingTable.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"
#include "itkVector.h"
#include <string>
#include <map>
#include "tinyxml.h"
#include "itkImageRegionIteratorWithIndex.h"

// ----- define specimen based classes and read functions ----//
class Specimen
{
public:
  unsigned long Id;
  std::string   Name;
  unsigned long ParentId;
  std::string   ParentName;
  itk::Vector<unsigned int> ParentCoord;

Specimen()
  {
  Id = 0;
  Name = "";
  ParentId = 0;
  ParentName = "";
  ParentCoord.Fill( 0 );
  }
};

typedef std::map< std::string, Specimen> SpecimenIdMapType;
int CreateSpecimenIdMap ( const char * csvFile, SpecimenIdMapType & map );

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
  if ( argc < 11 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "parameterXMLFile specimenIdCSVFile ";
    std::cout << "chunkTypeVolume chunkNumberVolume ";
    std::cout << "mriXMap mriYMap mriZMap bulkAffineTransform ";
    std::cout << "csvOutputDirectory volumeOutputDirectory [suffix]";
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
  std::string specimenCSVFile  = argv[2];
  std::string chunkTypeVolume = argv[3];
  std::string chunkNumberVolume = argv[4];
  std::string mriXMapVolume = argv[5];
  std::string mriYMapVolume = argv[6];
  std::string mriZMapVolume = argv[7];
  std::string bulkAffineTransform = argv[8];
  std::string csvOutputDirectory = argv[9];
  std::string volumeOutputDirectory = argv[10];

  std::string suffix = "";
  if ( argc > 11 )
    {
    suffix = argv[11];
    }

  SpecimenIdMapType specimenIdMap;
  Parameters p;
  std::string prefix;

  VolumeType::Pointer chunkType;
  VolumeType::Pointer chunkNumber;
  FloatVolumeType::Pointer mriXMap;
  FloatVolumeType::Pointer mriYMap;
  FloatVolumeType::Pointer mriZMap;
  TransformType::Pointer transform;

  ShortVolumeType::Pointer xVolume;
  ShortVolumeType::Pointer yVolume;
  ShortVolumeType::Pointer sidVolume;
  ShortVolumeType::Pointer cidVolume;
  ShortVolumeType::Pointer sNumberVolume;

  try
    {
    int err;

    err = CreateSpecimenIdMap( specimenCSVFile.c_str(),  specimenIdMap );
    if ( err ) { return EXIT_FAILURE; }

    std::cout << "No. specimens = " << specimenIdMap.size() << std::endl;

    err = FetchParametersFromXMLFile( parameterXMLFile.c_str(), p );
    if ( err ) { return EXIT_FAILURE; }

    prefix = p.BrainName;
    prefix += ".";
    prefix += p.SlabType;

    std::cout << prefix << std::endl;
    std::cout << specimenIdMap[prefix].Id << std::endl;
    std::cout << p.Offset << " " << p.IdOffset << std::endl;

    itk::idp::ReadImage<VolumeType>( chunkTypeVolume.c_str(), chunkType );
    itk::idp::ReadImage<VolumeType>( chunkNumberVolume.c_str(), chunkNumber );
    itk::idp::ReadImage<FloatVolumeType>( mriXMapVolume.c_str(), mriXMap );
    itk::idp::ReadImage<FloatVolumeType>( mriYMapVolume.c_str(), mriYMap );
    itk::idp::ReadImage<FloatVolumeType>( mriZMapVolume.c_str(), mriZMap );
    itk::idp::ReadTransform<TransformType>( bulkAffineTransform.c_str(), transform );

    RegionType region = chunkType->GetBufferedRegion();

    itk::idp::MakeImage<ShortVolumeType,VolumeType>( chunkType, xVolume, 0 );
    itk::idp::MakeImage<ShortVolumeType,VolumeType>( chunkType, yVolume, 0 );
    itk::idp::MakeImage<ShortVolumeType,VolumeType>( chunkType, sidVolume, 0 );
    itk::idp::MakeImage<ShortVolumeType,VolumeType>( chunkType, cidVolume, 0 );
    itk::idp::MakeImage<ShortVolumeType,VolumeType>( chunkType, sNumberVolume, 0 );

    Iterator titer( chunkType, region );
    Iterator niter( chunkNumber, region );
    FloatIterator xiter( mriXMap, region );
    FloatIterator yiter( mriYMap, region );
    FloatIterator ziter( mriZMap, region );
    ShortIterator x2iter( xVolume, region );
    ShortIterator y2iter( yVolume, region );
    ShortIterator siter( sidVolume, region );
    ShortIterator citer( cidVolume, region );
    ShortIterator sniter( sNumberVolume, region );

    unsigned long bid = specimenIdMap[p.BrainName].Id;

    int slice = -1;
    int slabNumber = 0;
    int cn = 0;
    int ct = 0;

    unsigned long sid = 0;
    unsigned long cid = 0;
    std::string cname = "";
    
    unsigned long count = 0;

    std::ofstream output;

    while( !titer.IsAtEnd() )
      {
      if ( titer.Get() )
        {

        IndexType index = titer.GetIndex();
        PointType point;
        chunkType->TransformIndexToPhysicalPoint( index, point );
        point = point + p.Offset;

        if ( slice != index[2] )
          {
          // update slab id and start new output file
          std::stringstream ss;
          ss << prefix << ".";
          ss.fill( '0' ); ss.width( 2 );
          ss << vnl_math_rnd( point[2] / chunkType->GetSpacing()[2] );
          ss << suffix;
          sid = specimenIdMap[ ss.str() ].Id;  
          slabNumber = vnl_math_rnd( point[2] / chunkType->GetSpacing()[2] );
          slice = index[2];
          std::cout << ss.str() << ": " << sid << std::endl;

          if ( output.is_open() )
            {
            output.close();
            } 

          std::string opath = csvOutputDirectory;
          opath += "/";
          opath += ss.str();
          opath += ".csv";

          output.open( opath.c_str(), std::ios::out );
          if ( output.fail() )
            {
            std::cout << "Error while opening output file " << opath << std::endl;
            return EXIT_FAILURE;
            }

          output << "brain_id,slab_id,chunk_id,mri_z,mri_y,mri_x,y,x" << std::endl;


          }
      
        int x = vnl_math_rnd( point[0] / chunkType->GetSpacing()[0] );
        int y = vnl_math_rnd( point[1] / chunkType->GetSpacing()[1] );

        if ( cn != niter.Get() || ct != titer.Get() )
          {
          // update chunk id
          std::stringstream ss;
          ss << prefix << ".";
          ss.fill( '0' ); ss.width( 2 );
          ss << vnl_math_rnd( point[2] / chunkType->GetSpacing()[2] ); 
          if ( p.SlabType != "bs" )
            {
            ss << ".";
            if ( titer.Get() == 40 )
              { ss << "cx" << "."; }
            else if ( titer.Get() == 80 )
              { ss << "s1" << "."; }
            else if ( titer.Get() == 120 )
              { ss << "s2" << "."; }
            ss.fill( '0' ); ss.width( 2 );
            ss << ( niter.Get() / 20 );
            }
          cid = specimenIdMap[ ss.str() ].Id;
          cn = niter.Get();
          ct = titer.Get();
          cname = ss.str();
          //std::cout << ss.str() << ": " << cid << std::endl;
          }

       PointType ip;

//       ip = point - p.Offset;
       ip[0] = xiter.Get();
       ip[1] = yiter.Get();
       ip[2] = ziter.Get();

       PointType op = transform->TransformPoint( ip );

       output << bid << "," << sid << "," << cid << ",";
       output << vnl_math_rnd( op[2] ) << ",";
       output << vnl_math_rnd( op[1] ) << ",";
       output << vnl_math_rnd( op[0] ) << ",";
       output << y << "," << x << std::endl;

       x2iter.Set( x );
       y2iter.Set( y );
       siter.Set( sid - p.IdOffset );
       citer.Set( cid - p.IdOffset );
       sniter.Set( slabNumber );

       count++;

        }

      ++titer;
      ++niter;
      ++xiter;
      ++yiter;
      ++ziter;
      ++x2iter;
      ++y2iter;
      ++siter;
      ++citer;
      ++sniter;
      }

    if ( output.is_open() )
      {
      output.close();
      }

    std::cout << count << std::endl;

    std::string opath;

    opath = volumeOutputDirectory;
    opath += "/x.mhd";
    itk::idp::WriteImage<ShortVolumeType>( opath.c_str(), xVolume );

    opath = volumeOutputDirectory;
    opath += "/y.mhd";
    itk::idp::WriteImage<ShortVolumeType>( opath.c_str(), yVolume );
    
    opath = volumeOutputDirectory;
    opath += "/slab_id.mhd";
    itk::idp::WriteImage<ShortVolumeType>( opath.c_str(), sidVolume );

    opath = volumeOutputDirectory;
    opath += "/chunk_id.mhd";
    itk::idp::WriteImage<ShortVolumeType>( opath.c_str(), cidVolume );

    opath = volumeOutputDirectory;
    opath += "/slab_number.mhd";
    itk::idp::WriteImage<ShortVolumeType>( opath.c_str(), sNumberVolume );


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


int CreateSpecimenIdMap ( const char * csvFile, SpecimenIdMapType & map )
{
  map.clear();

  std::ifstream gfile( csvFile, std::ios::in );
  if ( !gfile.fail() )
    {
    std::string buffer;

    // skip header line
    getline( gfile, buffer );

    while( getline( gfile, buffer ) )
      {

      if ( 0 )
        {  
        std::cout << buffer << std::endl;
        }

      std::istringstream instream;
      instream.str( buffer );
      std::string item;
      std::vector<std::string> components;
    
      // Parse comment separated text
      while( getline( instream, item, ',') )
        {
        components.push_back( item );
        }

      if( components.size() < 2 )
        {
        continue;
        }

      Specimen s;
      s.Id = atoi( components[0].c_str() );
      s.Name = components[1];

      if ( components.size() > 3 )
        {
        s.ParentId = atoi( components[2].c_str() );
        s.ParentName = components[3];
        }

      if ( components.size() > 6 )
        {
        s.ParentCoord[0] = atoi( components[4].c_str() );
        s.ParentCoord[1] = atoi( components[5].c_str() );
        s.ParentCoord[2] = atoi( components[6].c_str() );
        }
   
      map[s.Name] = s;

      } // while

    gfile.close();
    }
  else
    {
    return -1;
    }

   return 0;

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