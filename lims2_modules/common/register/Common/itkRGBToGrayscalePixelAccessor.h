/*=========================================================================

  itkRGBToGrayscalePixelAccessor.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkRGBToGrayscalePixelAccessor_h
#define __itkRGBToGrayscalePixelAccessor_h


#include "itkRGBPixel.h"


namespace itk
{

/**
 * \class RGBToGrayscalePixelAccessor
 * \brief Give access to the grayscale image form from a linear
 * combination of the red, green and blue channels. 
 *
 * This class is intended to be used as parameter of 
 * an ImageAdaptor to make an RGBPixel image appear as being
 * of scalar type T.
 *
 * \sa ImageAdaptor
 * \ingroup ImageAdaptors
 *
 */

template <class T>
class ITK_EXPORT RGBToGrayscalePixelAccessor
{
public:
  /** Standard class typedefs. */
  typedef   RGBToGrayscalePixelAccessor        Self;

  /** External typedef. It defines the external aspect
   * that this class will exhibit */
  typedef T ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data */
  typedef   RGBPixel<T>    InternalType;

  /** Write access to the Blue component */
  inline void Set( InternalType & output, const ExternalType & input ) const
    {
    output.SetRed( input );
    output.SetGreen( input );
    output.SetBlue( input ); 
    }

  /** Read access to the Blue component */
  inline ExternalType  Get( const InternalType & input ) const
    { return static_cast<ExternalType>( m_RedScale * input.GetRed() +
                                        m_GreenScale * input.GetGreen() +
                                        m_BlueScale * input.GetBlue() ); }

  inline void BlueOnly()
    {
    m_RedScale = 0.00; m_GreenScale = 0.00; m_BlueScale = 1.00;
    }

  inline void RedOnly()
    {
    m_RedScale = 1.00; m_GreenScale = 0.00; m_BlueScale = 0.00;
    }

  inline void GreenOnly()
    {
    m_RedScale = 0.00; m_GreenScale = 1.00; m_BlueScale = 0.00;
    }

  inline void PixelLuminance()
    {
    m_RedScale = 0.30;
    m_GreenScale = 0.59;
    m_BlueScale = 0.11;
    }

  inline void SetScales( double redScale, double greenScale, double blueScale )
    {
    m_RedScale = redScale;
    m_GreenScale = greenScale;
    m_BlueScale = blueScale;
    }

  bool operator!=( const Self & other ) const
    {
    return false;
    }
  bool operator==( const Self & other ) const
    {
    return !(*this != other);
    }

  RGBToGrayscalePixelAccessor()
    {
    m_RedScale = 0.30;
    m_GreenScale = 0.59;
    m_BlueScale = 0.11;
    };

  RGBToGrayscalePixelAccessor( const Self & other )
    {
    m_RedScale = other.RedScale;
    m_GreenScale = other.GreenScale;
    m_BlueScale = other.BlueScale;
    };

private:

  double  m_RedScale;
  double  m_BlueScale;
  double  m_GreenScale;

};
  
}  // end namespace itk

#endif
