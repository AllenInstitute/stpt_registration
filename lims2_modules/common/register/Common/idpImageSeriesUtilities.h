/*=========================================================================

  idpImageSeriesUtilities.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpImageSeriesUtilities_h
#define __idpImageSeriesUtilities_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "idpImageSeries.h"
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkCenteredAffineTransform.h"
#include <tinyxml.h>

namespace itk
{
namespace idp
{

/** \class ImageSeriesUtilities
 *
 */
template <class TPixel=unsigned char>
class ImageSeriesUtilities : public Object
{
public:

  /** Standard typedefs. */
  typedef ImageSeriesUtilities      Self;
  typedef Object        Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageSeriesUtilities, Object);
  
  /** Other types. */
  typedef TPixel PixelType;
  static const unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension > VolumeType;
  typedef typename VolumeType::Pointer VolumePointer;
  typedef itk::Image< unsigned char, Dimension > MaskVolumeType;
  typedef typename MaskVolumeType::Pointer MaskVolumePointer;
  typedef itk::Image< itk::RGBPixel<PixelType>, Dimension > RGBVolumeType;
  typedef typename VolumeType::SizeType SizeType;
  typedef typename VolumeType::SpacingType SpacingType;
  typedef typename VolumeType::PointType   PointType;

  typedef RGBPixel<PixelType>         RGBPixelType;
  typedef Image<RGBPixelType,2>       RGBImageType;
  typedef typename RGBImageType::Pointer       RGBImagePointer;
  typedef Image<PixelType,2>          ImageType;
  typedef typename ImageType::Pointer          ImagePointer;

  typedef ImageSeries::SubImageArray SubImageArray;

  typedef itk::Image< double, 2 > SizeVolumeType;
  typedef SizeVolumeType::Pointer SizeVolumePointer;
  typedef itk::Image< unsigned char, 2 > SizeMaskType;
  typedef SizeMaskType::Pointer SizeMaskPointer;

  typedef SubImage::TransformType TransformType;
  typedef SubImage::TransformPointer TransformPointer;
  typedef SubImage::Transform2DType Transform2DType;
  typedef SubImage::Transform2DPointer Transform2DPointer;
  typedef TransformType::ParametersType ParametersType;
  typedef std::vector<ParametersType>   ParametersVector;

  typedef itk::CenteredAffineTransform<double,3>  Affine3DTransformType;
  typedef Affine3DTransformType::Pointer Affine3DTransformPointer;

  /** Image series related types. */
  typedef ImageSeries::Pointer            ImageSeriesPointer;
  typedef std::vector<ImageSeriesPointer> ImageSeriesVector;

  /** Get volume region. */
  itkGetConstMacro( Size, SizeType );

  /** Get volume spacing. */
  itkGetConstMacro( Spacing, SpacingType );

  /** Get volume origin. */
  itkGetConstMacro( Origin, PointType );

  /** Get volume center. */
  itkGetConstMacro( Center, PointType );

  /** Set/Get the output width. */
  itkSetMacro( OutputWidth, double );
  itkGetConstMacro( OutputWidth, double );

  /** Get the output height. */
  itkSetMacro( OutputHeight, double );
  itkGetConstMacro( OutputHeight, double );

  /** Set/Get the model directory. */
  itkSetStringMacro( ModelDirectory );
  itkGetStringMacro( ModelDirectory );

  /** Set/Get image series. */
  itkSetObjectMacro( ImageSeries, ImageSeries );
  itkGetObjectMacro( ImageSeries, ImageSeries );

  /** Set/Get verbose flag. */
  itkSetMacro( Verbose, bool );
  itkGetConstMacro( Verbose, bool );

  /** Set/Get downsample factor. */
  itkSetMacro( DownsampleFactor, unsigned long )
  itkGetConstMacro( DownsampleFactor, unsigned long )

  /** Set/Get generate mask boolean. If true, mask file from each sub-image is used to generate a mask.
      If false, an all non-zero mask volume is generated.  */
  itkSetMacro( GenerateMask, bool );
  itkGetConstMacro( GenerateMask, bool );

  /** Set/Get use standard size boolean. If true, the reconstructed volume is set to the output
   * size in physical coordinates as define in "model.xm" in the model directory. If false,
   * the output size is set to be the maximum width/height of the image-series. */
  itkSetMacro( UseStandardSize, bool );
  itkGetConstMacro( UseStandardSize, bool );

  /** Set/Get flag to detect and invert dark field images. Only uses the blue channel assuming that it is DAPI. */
  itkSetMacro( InvertDarkFieldImages, bool );
  itkGetConstMacro( InvertDarkFieldImages, bool );

  /** Compute the meta information for the 3D volume.*/
  void ComputeVolumeMetaInformation();

  /** Compute centroid alignment. */
  void CentroidAlignment();

  /** Compute centroid alignment based on existing mask volume. */
  void CentroidAlignment( MaskVolumePointer & mask, PointType & target );

  /** Initialize volume to the right size and defaults. */
  void InitializeVolume( VolumePointer & volume, MaskVolumePointer & mask, MaskVolumePointer & sliceMask, PixelType imgBkgnd = 255 );  

  /** Populate volume by resampling sections. */
  void PopulateVolume( VolumePointer & volume, MaskVolumePointer & mask, MaskVolumePointer & sliceMask  );
  void PopulateVolume( VolumePointer & volume, MaskVolumePointer & mask, MaskVolumePointer & sliceMask, 
                       double rScale, double gScale, double bScale, PixelType imgBkgnd = 255 );
  /** Populate size volume. Generates tissue size data structure from tissue size
   * in the SubImage object. Use CentroidAlignment to computes tissue size. */
  void PopulateSizeVolume( SizeVolumePointer & sizeVolume, SizeMaskPointer & sizeMask );

  /** Populate size volume. Generates tissue size data structures from input volume. */
  static void PopulateSizeVolume( MaskVolumePointer & mask, 
                                  SizeVolumePointer & sizeVolume, SizeMaskPointer & sizeMask );

  /** Compose this transform to each subimage transform. */
  void ComposeTransform( const TransformType * transform );
  void ComposeTransform( const Transform2DType * transform );

  /** Compose an array of transform to each subimage transform. */
  void ComposeTransforms( const std::vector<TransformPointer> & transforms );
  void ComposeTransforms( const std::vector<Transform2DPointer> & transforms );

  /** Back up the current transform parameters. */
  void BackupTransformParameters();

  /** Revert to backup transform parameters. */
  void RevertTransformParameters();

  /** Read/Write transform parameters to XML. */
  static void WriteTransformParametersToXML( ImageSeriesVector & series, TiXmlNode * node );
  static void WriteTransformParametersToXML( ImageSeriesVector & series,
                                             const Affine3DTransformType * transform,
                                             const Affine3DTransformType * atog, 
                                             double metric, TiXmlNode * node  );
  void ReadTransformParametersFromXML( TiXmlNode * node, Affine3DTransformPointer & transform );
  void ReadTransformParametersFromXML( TiXmlNode * node );
  void ReadSubImageTransformParametersFromXML( TiXmlNode * node );


  /** Read in high resolution JP2 from file. */
  template <class TImage>
    void ReadImage( const SubImage * subimage, typename TImage::Pointer & image, int reduce = 0, int channel = -1 );

  template <class TImage>
    void ReadExpressionImage( const SubImage * subimage, typename TImage::Pointer & image, int reduce = 0, int channel = -1 );

  template <class TImage>
    void ReadProjectionMaskImage( const SubImage * subimage, typename TImage::Pointer & image, int reduce = 0, int channel = -1 );

  template <class TImage>
    void ReadImageFromFile( const std::string& filename, typename TImage::Pointer & image );



  /** Resample an image into a volume. */
  template <class TOutputPixel>
    void ResampleIntoVolume( const SubImage * subimage, const Image<TOutputPixel,2> * image, 
                             Image<TOutputPixel,3> * volume, const char * interpolationType = "Linear" );

  /** Set sub_image metric. */
  void SetSubImageMetric( const std::vector<double> & metric );


protected:
  ImageSeriesUtilities();
  ~ImageSeriesUtilities();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Get standard size information from XML. */
  void GetStandardSizeFromXML();


private:
  ImageSeriesUtilities(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ImageSeries::Pointer    m_ImageSeries;
  SpacingType             m_Spacing;
  PointType               m_Origin;
  SizeType                m_Size;
  PointType               m_Center;
  double                  m_OutputWidth;
  double                  m_OutputHeight;
  std::string             m_ModelDirectory;
  bool                    m_Verbose;
  ParametersVector        m_BackupParameters;
  unsigned long           m_DownsampleFactor;
  bool                    m_GenerateMask;
  bool                    m_UseStandardSize;
  bool                    m_InvertDarkFieldImages;


};

} // end namespace idp
} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "idpImageSeriesUtilities.txx"
#endif

#endif
