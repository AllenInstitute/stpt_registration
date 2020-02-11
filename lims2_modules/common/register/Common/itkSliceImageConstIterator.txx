/*=========================================================================

  itkSliceImageConstIterator.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _itkSliceImageConstIterator_txx
#define _itkSliceImageConstIterator_txx

#include "itkSliceImageConstIterator.h"
#include "itkImportImageFilter.h"

namespace itk
{

/* 
 * Return slice wrapped up as an itk image
 */
template< typename TImage >
const
typename SliceImageConstIterator<TImage>
::OutputImageType *
SliceImageConstIterator<TImage>
::GetOutput()
{

  typedef typename TImage::PixelType PixelType;
  typedef ImportImageFilter< PixelType, OutputImageDimension > ImporterType;
  typename ImporterType::Pointer importer = ImporterType::New();
 
  importer->SetSpacing( m_OutputSpacing );
  importer->SetOrigin( m_OutputOrigin );
  importer->SetRegion( m_OutputRegion );

  importer->SetImportPointer( const_cast<PixelType *>(m_Buffer + m_Offset), m_PixelPerSlice, false );

  try
    {
    importer->Update();
    m_OutputImage = importer->GetOutput();
    m_OutputImage->DisconnectPipeline();
    }
  catch( itk::ExceptionObject & err )
    {
    throw err;
    }
  catch( ... )
    {
    itkExceptionMacro( << "Caught unknown error" );
    }

  return m_OutputImage;

}


} // namespace itk

#endif
