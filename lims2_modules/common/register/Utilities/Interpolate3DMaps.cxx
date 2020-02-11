/*=========================================================================

  Interpolate3DMaps.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "idpRegistrationUtilities.h"
#include <string>
#include <vector>
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkReinitializeLevelSetImageFilter.h"
#include <map>

typedef itk::Image<unsigned short, 3> VolumeType;
typedef itk::Image<float, 3> FloatVolumeType;
typedef itk::ImageRegionIterator<VolumeType> Iterator;
typedef itk::ImageRegionIterator<FloatVolumeType> FloatIterator;
typedef std::vector<VolumeType::PixelType> ArrayType;

// --- helper function to generate signed distance map --- //
int SignedDistance(
FloatVolumeType::Pointer & input,
FloatVolumeType::Pointer & output,
float maximumDistance = 30 );

// -- helper function to resolve interpolation using distance maps within mask --- //
int ResolveSpecimen(
VolumeType::Pointer & idVolume,
VolumeType::Pointer & mask,
VolumeType::Pointer & outputVolume,
ArrayType & uniqueIds );

// -- helper function to compute propagation field -- //
int PropagationField(
VolumeType::Pointer & idVolume,
VolumeType::Pointer & mask,
FloatVolumeType::Pointer & fieldx,
FloatVolumeType::Pointer & fieldy,
FloatVolumeType::Pointer & fieldz );


int main( int argc, char *argv[] )
{

  if ( argc < 9 )
    {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "slabIDVolume chunkIDVolume MaskVolume ";
    std::cout << "outputSlabID outputChunkId fieldx fieldy fieldz [maxdist]";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }

  std::string slabIdFile = argv[1];
  std::string chunkIdFile = argv[2];
  std::string maskFile = argv[3];
  std::string outputSlabIdFile = argv[4];
  std::string outputChunkIdFile = argv[5];

  std::vector< std::string > fieldFiles;
  for ( unsigned int i = 0; i < 3; i++ )
    {
    fieldFiles.push_back( argv[i+6] );
    }

  float maxdist = 30;
  if ( argc > 9 )
    {
    maxdist = atof( argv[9] );
    }


  try
    {
      int err;

      VolumeType::Pointer inputSlabId;
      VolumeType::Pointer inputChunkId;
      VolumeType::Pointer inputMask;

      VolumeType::Pointer outputSlabId;
      VolumeType::Pointer outputChunkId;
      std::vector< FloatVolumeType::Pointer > field( fieldFiles.size() );

      // read in the files
      itk::idp::ReadImage<VolumeType>( slabIdFile.c_str(), inputSlabId );
      itk::idp::ReadImage<VolumeType>( chunkIdFile.c_str(), inputChunkId );
      itk::idp::ReadImage<VolumeType>( maskFile.c_str(), inputMask );


      if ( inputSlabId.IsNull() )
        {
        std::cout << "Error: Can not open slab id volume " << slabIdFile << std::endl;
        return EXIT_FAILURE;
        }

      if( inputChunkId.IsNull() )
        {
        std::cout << "Error: Can not open chunk id volume " << chunkIdFile << std::endl;
        return EXIT_FAILURE;
        }

      if( inputMask.IsNull() )
        {
        std::cout << "Error: Can not open mask volume " << maskFile << std::endl;
        return EXIT_FAILURE;
        }

      // -- resolve slabs first using structure mask ---- //
      ArrayType uniqueSlabIds;
      std::cout << "Resolving slabs with global mask" << std::endl;
      itk::idp::MakeImage<VolumeType,VolumeType>( inputSlabId, outputSlabId, 0 );
      err = ResolveSpecimen( inputSlabId, inputMask, outputSlabId, uniqueSlabIds );
      if ( err )
        {
        std::cout << "Error: while resolving slabs " << std::endl;
        return EXIT_FAILURE;
        }
      itk::idp::WriteImage<VolumeType>( outputSlabIdFile.c_str(), outputSlabId );

      // -- for each slab resolve chunk using slab mask --- //
      itk::idp::MakeImage<VolumeType,VolumeType>( inputSlabId, outputChunkId, 0 );
      ArrayType uniqueChunkIds;

      std::cout << "Resolving chunks within created slab mask" << std::endl;
      for ( unsigned int s = 0; s < uniqueSlabIds.size(); s++ )
        {
        std::cout << "slab: " << uniqueSlabIds[s] << std::endl;
        VolumeType::Pointer mask;
        ArrayType ids;
        itk::idp::BinaryThreshold<VolumeType>( outputSlabId, mask, 
                                               uniqueSlabIds[s], uniqueSlabIds[s], 
                                               uniqueSlabIds[s], 0 );
        err = ResolveSpecimen( inputChunkId, mask, outputChunkId, ids );
        if ( err )
          {
          std::cout << "Error: while resolving chunks in slab" << std::endl;
          return EXIT_FAILURE;
          }
        for( unsigned int c = 0; c < ids.size(); c++ )
          {
          uniqueChunkIds.push_back( ids[c] );
          }
        }

      itk::idp::WriteImage<VolumeType>( outputChunkIdFile.c_str(), outputChunkId );

      // -- for each resolve chunk propogate information --- //
      for( unsigned int j = 0; j < 3; j++ )
        {
        itk::idp::MakeImage<FloatVolumeType,VolumeType>( inputSlabId, field[j], 0.0 );
        }

      std::cout << "Create propgation field within each create chunk mask" << std::endl;
      for( unsigned int s = 0; s < uniqueChunkIds.size(); s++ )
        {
        std::cout << "chunk: " << uniqueChunkIds[s] << std::endl;
        VolumeType::Pointer mask;
        itk::idp::BinaryThreshold<VolumeType>( outputChunkId, mask, 
                                               uniqueChunkIds[s], uniqueChunkIds[s], 
                                               uniqueChunkIds[s], 0 );
        err = PropagationField( inputChunkId, mask, field[0], field[1], field[2] );
        if ( err )
          {
          std::cout << "Error: while propagating in chunk" << std::endl;
          return EXIT_FAILURE;
          }
        }

      itk::idp::WriteImage<FloatVolumeType>( fieldFiles[0].c_str(), field[0] );
      itk::idp::WriteImage<FloatVolumeType>( fieldFiles[1].c_str(), field[1] );
      itk::idp::WriteImage<FloatVolumeType>( fieldFiles[2].c_str(), field[2] );

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



// -- helper function to resolve specimens with a mask --
int ResolveSpecimen(
VolumeType::Pointer & idVolume,
VolumeType::Pointer & mask,
VolumeType::Pointer & outputVolume,
ArrayType & values )
{

  typedef std::map<VolumeType::PixelType, VolumeType::PixelType> MapType;
  MapType uniqueIds;

  // iterate through the id Volume to find all unique ids
  Iterator iter( idVolume, idVolume->GetBufferedRegion() );
  Iterator kiter( mask, mask->GetBufferedRegion() );
  while ( !iter.IsAtEnd() )
    {

    MapType::iterator mend = uniqueIds.end();

    if ( iter.Get() && kiter.Get()  )
      {
      MapType::iterator miter = uniqueIds.find( iter.Get() );
      if (miter == mend)
        {
        uniqueIds[iter.Get()] = iter.Get();
        mend = uniqueIds.end();
        //std::cout << "new id: " << iter.Get() << std::endl;
        }
      }
    ++iter;
    ++kiter;
    }

  std::cout << "no. ids: " << uniqueIds.size() << std::endl;

  std::vector<FloatVolumeType::Pointer> distances;
  values.clear();

  // for each unique id create mask file and compute signed distance
  for ( MapType::iterator iter = uniqueIds.begin(); iter != uniqueIds.end(); iter++ )
    {
      std::cout << "distance: " << (*iter).first << std::endl;
      FloatVolumeType::Pointer dvol;
      {
      typedef itk::BinaryThresholdImageFilter<VolumeType,FloatVolumeType> FilterType;
      FilterType::Pointer filter = FilterType::New();
      filter->SetInput( idVolume );
      filter->SetUpperThreshold( (*iter).first );
      filter->SetLowerThreshold( (*iter).first );
      filter->SetInsideValue( -0.5 );
      filter->SetOutsideValue( 0.5 );
      filter->Update();
      dvol = filter->GetOutput();
      dvol->DisconnectPipeline();
      }

      int err = SignedDistance( dvol, dvol, 30 );
      if ( err ) { return err; }
      distances.push_back( dvol );
      values.push_back( (*iter).first );

    //  char buffer[200]= "";
    //  sprintf( buffer, "distance%d.mhd",  (*iter).first  );
    //  itk::idp::WriteImage< FloatVolumeType >( buffer, dvol );
    }

  // populate the output volume
  VolumeType::RegionType region = outputVolume->GetBufferedRegion();
  Iterator oiter( outputVolume, region );
  Iterator miter( mask, region );
  std::vector< FloatIterator > diters( values.size() );

  for ( unsigned int i = 0; i < values.size(); i++ )
    {
    diters[i] = FloatIterator( distances[i], region );
    }

  while ( !oiter.IsAtEnd() )
    {
    if ( miter.Get() )
      {
     // std::cout << miter.GetIndex() << std::endl;
      float minValue = itk::NumericTraits<float>::max();
      unsigned int minIndex = values.size();
      for ( unsigned int i = 0; i < values.size(); i++ )
        {
        if ( diters[i].Get() < minValue )
          {
          minValue = diters[i].Get();
          minIndex = i;
          }
        }
      oiter.Set( values[minIndex] );
      }

    ++oiter;
    ++miter;
    for( unsigned int i = 0; i < values.size(); i++ )
       {
       ++(diters[i]);
       }
    }

  return 0;

}


// -- helper function to compute propagation field -- //
int PropagationField(
VolumeType::Pointer & idVolume,
VolumeType::Pointer & mask,
FloatVolumeType::Pointer & fieldx,
FloatVolumeType::Pointer & fieldy,
FloatVolumeType::Pointer & fieldz )
{
  VolumeType::RegionType region = idVolume->GetBufferedRegion();
  Iterator iiter( idVolume, region );
  Iterator miter( mask, region );

  std::vector< VolumeType::IndexType > ilist;
  std::vector< VolumeType::PointType > plist;

  while ( !iiter.IsAtEnd() )
    {
    if ( iiter.Get() && miter.Get() )
      {
      VolumeType::IndexType index = iiter.GetIndex();
      VolumeType::PointType point;
      idVolume->TransformIndexToPhysicalPoint( index, point );
      ilist.push_back( index );
      plist.push_back( point );
      }
    ++iiter;
    ++miter;
    }

  if ( ilist.size() == 0 ) 
    {
    return 0;
     }
  std::cout << "no. points: " << ilist.size() << std::endl;

  iiter.GoToBegin();
  miter.GoToBegin();
  FloatIterator xiter( fieldx, region );
  FloatIterator yiter( fieldy, region );
  FloatIterator ziter( fieldz, region );

  while ( !iiter.IsAtEnd() )
    {
    if ( miter.Get() )
      {
      VolumeType::IndexType index = iiter.GetIndex();
      VolumeType::PointType point;
      idVolume->TransformIndexToPhysicalPoint( index, point );

      float minDist = point.SquaredEuclideanDistanceTo( plist[0] );
      VolumeType::PointType minPoint = plist[0];

      for( unsigned int i = 0; i < ilist.size(); i++ )
        {
        float d = point.SquaredEuclideanDistanceTo( plist[i] );
        if ( d < minDist )
          {
          minDist = d;
          minPoint = plist[i];
          }
        }

      xiter.Set( minPoint[0] );
      yiter.Set( minPoint[1] );
      ziter.Set( minPoint[2] );
      
      }
    ++iiter;
    ++miter;
    ++xiter;
    ++yiter;
    ++ziter;
    }

  return 0;

}


// -- helper function to compute signed distance
int SignedDistance(
FloatVolumeType::Pointer & input,
FloatVolumeType::Pointer & output,
float maximumDistance )
{

  typedef itk::ReinitializeLevelSetImageFilter<FloatVolumeType> DistanceType;

  DistanceType::Pointer distance = DistanceType::New();

  distance->SetInput( input );
  distance->NarrowBandingOn();
  distance->SetNarrowBandwidth( 2 * maximumDistance );

  try
    {
    distance->Update();
    output = distance->GetOutput();
    output->DisconnectPipeline();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return -1;
    }
  catch( ... )
    {
    std::cerr << "Caught unknown exception" << std::endl;
    return -1;
    }

  return 0;

}