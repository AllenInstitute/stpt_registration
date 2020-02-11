/*=========================================================================

  idpRegistrationUtilities.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpRegistrationUtilities_h
#define __idpRegistrationUtilities_h

#include <stddef.h>

#include "itkExceptionObject.h"
#include "itkInterpolateImageFunction.h"
#include "itkTransform.h"
#include "itkCommand.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include <itksys/SystemTools.hxx>

namespace itk
{
namespace idp
{

// helper function to change between linux and windows paths
void ChangePaths( std::string & str );



// helper function to read image from file
template <typename ImageType>
void
ReadImage( 
const char * filename,
typename ImageType::Pointer & image ) 
throw (ExceptionObject);

// helper function to write image to file
template <typename ImageType>
void
WriteImage( 
const char * filename,
const ImageType * image,
bool compression = false ) 
throw( ExceptionObject);

// helper function make an image using a reference
template <typename ImageType, typename RefType>
void
MakeImage(
typename RefType::Pointer & ref,
typename ImageType::Pointer & output,
typename ImageType::PixelType defaultPixel
)
throw( ExceptionObject );

// helper function to resample image
template <typename ImageType>
void
ResampleImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & ref,
const Transform< double, ImageType::ImageDimension, ImageType::ImageDimension> * trans,
typename ImageType::Pointer & output,
typename ImageType::PixelType pad = 0,
const char * interpolatorType = "Linear"
 )
throw( ExceptionObject );

// helper function to resample image
template <typename ImageType>
void
VectorResampleImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & ref,
const Transform< double, ImageType::ImageDimension, ImageType::ImageDimension> * trans,
typename ImageType::Pointer & output,
typename ImageType::PixelType pad = 0,
const char * interpolatorType = "VectorLinear"
 )
throw( ExceptionObject );

// helper function to compute the maximum projection image
template <typename ImageType>
void
MaximumProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int projectionDimension = 2
 )
throw( ExceptionObject );

// helper function to compute the minimum projection image
template <typename ImageType>
void
MinimumProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int projectionDimension = 2
)
throw( ExceptionObject );

// helper function to compute the average projection image
template <typename ImageType>
void
AverageProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int projectionDimension = 2
)
throw( ExceptionObject );

// helper function to compute the average projection image
template <typename ImageType, typename MaskImageType>
void
MaskedAverageProjection(
typename ImageType::Pointer & input,
typename MaskImageType::Pointer & mask,
typename ImageType::Pointer & output,
unsigned int projectionDimension = 2
)
throw( ExceptionObject );

// helper function to compute the median projection image
template <typename ImageType>
void
MedianProjection(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
unsigned int projectionDimension = 2
)
throw( ExceptionObject );



// copy region from one image to another
template <typename ImageType>
void
CopyRegion(
typename ImageType::Pointer & fromImage,
const typename ImageType::RegionType & fromRegion,
typename ImageType::Pointer & toImage,
const typename ImageType::RegionType & toRegion
)
throw( ExceptionObject );



// helper function to compute the correlation between two images
template <typename ImageType, typename MaskImageType>
void
ComputeCorrelation( const ImageType * image1,
                    const ImageType * image2,
                    const MaskImageType * mask1,
                    const MaskImageType * mask2,
                    double & value ) 
throw (ExceptionObject);

// helper function to invert the intensity
template <typename ImageType>
void
InvertIntensity(
const ImageType * input,
typename ImageType::Pointer & output
 )
throw( ExceptionObject );



// helper function to extract one slice
template <typename InputImageType, typename OutputImageType>
void
ExtractSlice(
typename InputImageType::Pointer & input,
typename OutputImageType::Pointer & output,
unsigned int slice,
unsigned int dimension
)
throw( ExceptionObject );

// helper function to extract a region of interest
template <typename ImageType>
void
ExtractRegion(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
typename ImageType::RegionType  region
)
throw( ExceptionObject );

// helper function to pad an image
template <typename ImageType>
void
PadImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
typename ImageType::SizeType  upperExtendSize,
typename ImageType::SizeType  lowerExtendSize,
typename ImageType::PixelType  padValue        
)
throw( ExceptionObject );

// helper function to fill a specifiy region
template <typename ImageType>
void
FillRegion(
typename ImageType::Pointer & ref,
typename ImageType::Pointer & output,
typename ImageType::RegionType region,
typename ImageType::PixelType foreground,
typename ImageType::PixelType background
)
throw( ExceptionObject );

// helper function to change the spacing of a volume
template <typename ImageType>
void
ChangePixelSpacing(
typename ImageType::Pointer & input,
const FixedArray<double,ImageType::ImageDimension> & spacing,
const FixedArray<double,ImageType::ImageDimension> & variance,
typename ImageType::Pointer & output,
const char * interpolatorType = "Linear",
typename ImageType::PixelType pad = 0
)
throw( ExceptionObject );

// helper function to threshold the image
template <typename ImageType>
void
BinaryThreshold(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
typename ImageType::PixelType lower,
typename ImageType::PixelType upper,
typename ImageType::PixelType inside,
typename ImageType::PixelType outside
 )
throw( ExceptionObject );

// helper function to cast image type to another
template <typename InputImageType, typename OutputImageType>
void
CastImage(
typename InputImageType::Pointer & input,
typename OutputImageType::Pointer & output
)
throw( ExceptionObject );


// helper function to write transform to file
template< typename TransformType >
void WriteTransform(
const char * filename,
typename TransformType::Pointer & transform
)
throw( ExceptionObject );



// helper function to read transform from file
template< typename TransformType >
void ReadTransform(
const char * filename,
typename TransformType::Pointer & transform
)
throw( ExceptionObject );



// helper function to mask an image
template <typename ImageType, typename MaskType>
void
MaskImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & mask,
typename ImageType::Pointer & output,
typename ImageType::PixelType outsideValue
)
throw( ExceptionObject );


// helper function to gaussian smooth an image
template <typename ImageType>
void
GaussianSmoothImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
const double * variance 
)
throw( ExceptionObject );

// helper function to divide one image by another
template <typename ImageType>
void
DivideImage(
typename ImageType::Pointer & input1,
typename ImageType::Pointer & input2,
typename ImageType::Pointer & output
)
throw( ExceptionObject );

// helper function to add one image by another
template <typename ImageType>
void
AddImage(
typename ImageType::Pointer & input1,
typename ImageType::Pointer & input2,
typename ImageType::Pointer & output,
bool inPlace = false
)
throw( ExceptionObject );

// help function to compute statistics
template <typename ImageType>
void
ImageStatistics(
typename ImageType::Pointer & input,
double & minimum,
double & maximum,
double & sum,
double & mean,
double & variance,
double & sigma
)
throw( ExceptionObject );

// helper function to not an image
template <typename ImageType>
void
NotImage(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output
)
throw( ExceptionObject );

// helper function to shift and scale an image
template <typename InputImageType, typename OutputImageType>
void
ShiftScale(
typename InputImageType::Pointer & input,
typename OutputImageType::Pointer & output,
double shift,
double scale
)
throw( ExceptionObject );

// select on component out of a vector
template <typename InputImageType, typename OutputImageType>
void
VectorIndexSelection( 
typename InputImageType::Pointer & input,
unsigned int index,
typename OutputImageType::Pointer & output
) throw( ExceptionObject );

// helper function transform pixel intensity by log10
template <typename ImageType>
void
Log10Image(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output
)
throw( ExceptionObject );

// helper function to threshold below
template <typename ImageType>
void
ThresholdBelow(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
double thresh,
double outsideValue
)
throw( ExceptionObject );

// helper function to threshold above
template <typename ImageType>
void
ThresholdAbove(
typename ImageType::Pointer & input,
typename ImageType::Pointer & output,
double thresh,
double outsideValue
)
throw( ExceptionObject );

// Resample image using deformation field
template <typename ImageType, typename RefImageType, typename DeformationFieldType>
void ResampleImageByDfmfld(
typename ImageType::Pointer & input,
typename RefImageType::Pointer & ref,
typename DeformationFieldType::Pointer & dfmfld,
typename ImageType::Pointer & output,
const char * interpolatorType
 )
throw( ExceptionObject );



class Relation
{
public:
	double weight;
	int vert1;
	int vert2;
	Relation(double w=0,int v1=0,int v2=0);
	bool operator<(const Relation & b) const;
	void print();
};


int newRelation(Relation edge1,Relation& edgeout,bool* vertices);

} // end namespace idp
} //end namespace itk
           
#ifndef ITK_MANUAL_INSTANTIATION
#include "idpRegistrationUtilities.txx"
#endif

#endif

  
