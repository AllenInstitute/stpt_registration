/*=========================================================================

  itkSliceImageConstIterator.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkSliceImageConstIterator_h
#define _itkSliceImageConstIterator_h

#include "itkImage.h"

namespace itk
{

/** \class SliceImageIterator
 * \brief Helper class to iterate through the export of one (N-1)D slice image
 * at a time.
 */
template < typename TImage >
class ITK_EXPORT SliceImageConstIterator
{
public:

  /** Standard class typedefs. */
  typedef SliceImageConstIterator Self;

  /** Dimension of the input image the iterator walks.   */
  itkStaticConstMacro(ImageIteratorDimension, unsigned int,
                      TImage::ImageDimension );

  /** Dimension of the output wrapped (N-1)D image. */
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TImage::ImageDimension - 1);

  virtual const char *GetNameOfClass() const 
    {return "SliceImageConstIterator";}

  /** Image typedef support. */
  typedef TImage                              ImageType;
  typedef typename TImage::PixelType          PixelType;
  typedef itk::Image<PixelType, OutputImageDimension> OutputImageType;

  /** Default Constructor. Need to provide a default constructor since we
   * provide a copy constructor. */
  SliceImageConstIterator() : m_InputRegion(), m_OutputRegion()
    {
    m_Image = 0;
    m_Buffer = 0;
    m_Offset = 0;
    m_BeginOffset = 0;
    m_EndOffset = 0;

    m_PixelPerSlice = 0;
    m_OutputSpacing.Fill( 1.0 );
    m_OutputOrigin.Fill( 0.0 );
    m_StartIndex = 0;

    m_OutputImage = 0;
    }

  /** Default Destructor. */
  virtual ~SliceImageConstIterator() {};

  /** Constructore establishes an iterator to walk a particular image. */
  SliceImageConstIterator( const ImageType * ptr )
    {
    m_Image = ptr;
    m_Buffer = m_Image->GetBufferPointer();
    
    m_InputRegion = m_Image->GetBufferedRegion(); 
    m_StartIndex = m_InputRegion.GetIndex( ImageIteratorDimension - 1 );

    for( unsigned int j = 0; j < OutputImageDimension; j++ )
      {
      m_OutputRegion.SetIndex( j, m_InputRegion.GetIndex(j) );
      m_OutputRegion.SetSize( j, m_InputRegion.GetSize(j) );
      m_OutputSpacing[j] = m_Image->GetSpacing()[j];
      m_OutputOrigin[j] = m_Image->GetOrigin()[j];
      }   
   
    m_PixelPerSlice = m_OutputRegion.GetNumberOfPixels();

    m_Offset = 0;
    m_BeginOffset = m_Offset;

    if( m_InputRegion.GetNumberOfPixels() == 0 )
      {
      m_EndOffset = m_Offset;
      }
    else
      {
      m_EndOffset = m_InputRegion.GetSize( ImageIteratorDimension - 1 ) *
                    m_PixelPerSlice;
      }

    m_OutputImage = 0;

    }

  /** Copy Constructor. */
  SliceImageConstIterator( const Self& it )
    {
    m_Image = it.m_Image;
    m_Buffer = it.m_Buffer;
    m_Offset = it.m_Offset;
    m_BeginOffset = it.m_BeginOffset;
    m_EndOffset = it.m_EndOffset;

    m_PixelPerSlice = it.m_PixelPerSlice;
    m_OutputSpacing = it.m_OutputSpacing;
    m_OutputOrigin = it.m_OutputOrigin;

    m_InputRegion = it.m_InputRegion;
    m_OutputRegion = it.m_OutputRegion;

    m_StartIndex = it.m_StartIndex;

    m_OutputImage = 0;
    }


  /** Comparison operator. Two iterators are the same if they "point to" the
   * same memory location */
  bool
  operator!=(const Self &it) const
    {
    // two iterators are the same if they "point to" the same memory location
    return (m_Buffer + m_Offset) != (it.m_Buffer + it.m_Offset);
    };

  /** Comparison operator. Two iterators are the same if they "point to" the
   * same memory location */
  bool
  operator==(const Self &it) const
    {
    // two iterators are the same if they "point to" the same memory location
    return (m_Buffer + m_Offset) == (it.m_Buffer + it.m_Offset);
    };

  /** Comparison operator. An iterator is "less than" another if it "points to"
   * a lower memory location. */
  bool
  operator<=(const Self &it) const
    {
    // an iterator is "less than" another if it "points to" a lower
    // memory location
    return (m_Buffer + m_Offset) <= (it.m_Buffer + it.m_Offset);
    };

  /** Comparison operator. An iterator is "less than" another if it "points to"
   * a lower memory location. */
  bool
  operator<(const Self &it) const
    {
    // an iterator is "less than" another if it "points to" a lower
    // memory location
    return (m_Buffer + m_Offset) < (it.m_Buffer + it.m_Offset);
    };

  /** Comparison operator. An iterator is "greater than" another if it
   * "points to" a higher location. */
  bool
  operator>=(const Self &it) const
    {
    // an iterator is "greater than" another if it "points to" a higher
    // memory location
    return (m_Buffer + m_Offset) >= (it.m_Buffer + it.m_Offset);
    };

  /** Comparison operator. An iterator is "greater than" another if it
   * "points to" a higher location. */
  bool
  operator>(const Self &it) const
    {
    // an iterator is "greater than" another if it "points to" a higher
    // memory location
    return (m_Buffer + m_Offset) > (it.m_Buffer + it.m_Offset);
    };

  /** Get the slice index. This provides a read only reference to the index.
   * This causes the index to be calculated from pointer arithmetic and is
   * therefore an expensive operation.
   * \sa SetSliceIndex */
  long GetSliceIndex() const
    { 
    return ( ( m_Offset / m_PixelPerSlice ) + m_StartIndex );
    }

  /** Set the slice index. No bounds checking is performed.
   * \sa GetSliceIndex */
  virtual void SetSliceIndex( long ind )
    { 
    m_Offset = ( ind - m_StartIndex ) * m_PixelPerSlice;
    }
   
  /** Get the output image */
  virtual const OutputImageType * GetOutput();

  /** Move an iterator to the beginning of the region. "Begin" is
   * defined as the first pixel in the region. */
  void GoToBegin()
    {
    m_Offset = m_BeginOffset;
    }

  /** Move an iterator to the end of the region. "End" is defined as
   * one pixel past the last pixel of the region. */
  void GoToEnd()
    {
    m_Offset = m_EndOffset;
    }

  /** Is the iterator at the beginning of the region? "Begin" is defined
   * as the first pixel in the region. */
  bool IsAtBegin(void) const
    {
    return (m_Offset == m_BeginOffset);
    }

  /** Is the iterator at the end of the region? "End" is defined as one
   * pixel past the last pixel of the region. */
  bool IsAtEnd(void) const
    {
    return (m_Offset == m_EndOffset);
    }

  /** Increment (prefix) */
  Self &
  operator++()
    {
    this->m_Offset += m_PixelPerSlice;
    return *this;
    }

  /** Decrement (prefix) */
  Self & operator--()
    {
    this->m_Offset -= m_PixelPerSlice;
    return *this;
    }

protected: 
  typename TImage::ConstWeakPointer   m_Image;
  unsigned long                       m_Offset;
  unsigned long                       m_BeginOffset;
  unsigned long                       m_EndOffset;
  const PixelType                  *  m_Buffer;

  typename TImage::RegionType           m_InputRegion;
  long                                  m_StartIndex;
  typename OutputImageType::RegionType  m_OutputRegion;
  typename OutputImageType::SpacingType m_OutputSpacing;
  typename OutputImageType::PointType   m_OutputOrigin;

  unsigned long                         m_PixelPerSlice;
  typename OutputImageType::Pointer     m_OutputImage;

};

} // namespace itk

#if ITK_TEMPLATE_TXX
#include "itkSliceImageConstIterator.txx"
#endif

#endif
