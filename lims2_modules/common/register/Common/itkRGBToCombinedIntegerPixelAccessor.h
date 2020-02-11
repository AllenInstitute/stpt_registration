/*=========================================================================

  itkRGBToCombinedIntegerPixelAccessor.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __itkRGBToCombinedIntegerPixelAccessor_h
#define __itkRGBToCombinedIntegerPixelAccessor_h


#include "itkRGBPixel.h"


namespace itk
{

/**
 * \class RGBToCombinedIntegerPixelAccessor
 * \brief Give access to the grayscale image form from a concatentating 
 * the red, green and blue channels. 
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
class ITK_EXPORT RGBToCombinedIntegerPixelAccessor
{
public:
  /** Standard class typedefs. */
  typedef   RGBToCombinedIntegerPixelAccessor        Self;

  /** External typedef. It defines the external aspect
   * that this class will exhibit */
  typedef T ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data */
  typedef   RGBPixel<unsigned char>    InternalType;

  /** Write access to the combined component */
  inline void Set( InternalType & output, const ExternalType & input ) const
    {    
    T rem = input % 65536;
    output.SetRed( input / 65536 );
    output.SetGreen( rem / 256 );
    rem = rem % 256;
    output.SetBlue( rem ); 
    }

  /** Read access to the combined component */
  inline ExternalType  Get( const InternalType & input ) const
    { return static_cast<ExternalType>( 65536 * input.GetRed() +
                                        256 * input.GetGreen() +
                                        input.GetBlue() ); }


  bool operator!=( const Self & other ) const
    {
    return false;
    }
  bool operator==( const Self & other ) const
    {
    return !(*this != other);
    }


};
  
}  // end namespace itk

#endif
