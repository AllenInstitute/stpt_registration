/*=========================================================================

  itkJP2ImageIO.h

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#ifndef __itkJP2ImageIO_h
#define __itkJP2ImageIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <stddef.h>
#include <fstream>
#include "itkImageIOBase.h"
#include <stdio.h>
#include "itkImageRegion.h"
#include "itkArray.h"

class kdu_thread_env;


namespace itk
{
	
	/** 
	 * \class JP2ImageIO class 
	 * \brief Read JPEG2000 Image file format
	 *
	 * This is simple wrapper class for the Kakadu JPEG2000 library 
	 *
	 *  \ingroup IOFilters
	 *
	 */
	class ITK_EXPORT JP2ImageIO : public ImageIOBase
	{
	public:
		/** Standard class typedefs. */
		typedef JP2ImageIO            Self;
		typedef ImageIOBase  Superclass;
		typedef SmartPointer<Self>  Pointer;
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** Run-time type information (and related methods). */
		itkTypeMacro(JP2ImageIO, Superclass);
		
		/*-------- This part of the interfaces deals with reading data. ----- */
		
		/** Determine the file type. Returns true if this ImageIO can read the
		 * file specified. */
		virtual bool CanReadFile(const char*) ;
		
		/** Set the spacing and dimension information for the set filename. */
		virtual void ReadImageInformation();
		
		/** Reads the data from disk into the memory buffer provided. */
		virtual void Read(void* buffer);
		
		/** Set/Get the transpose flag. If Transpose is true, vertical coordinates
		 * is transposed with the horizontal coordinates. */
		itkSetMacro( Transpose, bool );
		itkGetMacro( Transpose, bool );
		itkBooleanMacro( Transpose );
		
		/** Set/Get the vertical flip flag. If true, the image is flipped vertically. */
		itkSetMacro( VerticalFlip, bool );
		itkGetMacro( VerticalFlip, bool );
		itkBooleanMacro( VerticalFlip );
		
		/** Set/Get the horizonal flip flag. If true, the image is flipped horizontally. */
		itkSetMacro( HorizontalFlip, bool );
		itkGetMacro( HorizontalFlip, bool );
		itkBooleanMacro( HorizontalFlip );
		
		/** Set/Get the color channel compression flag. If true, the image's color channels are compressed independently.
		 I.e., the luminance/chrominance color space conversion is skipped. */
		itkSetMacro( KeepColorChannelsSeparate, bool );
		itkGetMacro( KeepColorChannelsSeparate, bool );
		itkBooleanMacro( KeepColorChannelsSeparate );
		
		/** Set/Get the image reduction. Image dimensions are essentially divided by
		 * 2 to the power of this number. This number argument affects the apparent
		 * dimensions and number of DWT levels. */
		itkSetMacro( Reduce, int );
		itkGetMacro( Reduce, int ); 
		
		/** There are two way of setting the region of interest for decomposition.
		 * The first method is set the region directly using SetRegionOfInterest().
		 * The second method is to specifiy the ROI as fractions of the original
		 * image. The latter method is used if UseFractionROISpecification is true. */
		itkSetMacro( UseFractionROISpecification, bool );
		itkGetMacro( UseFractionROISpecification, bool );
		itkBooleanMacro( UseFractionROISpecification );
		
		/** Image region typedef. */
		typedef ImageRegion<2> RegionType;
		
		/** Set/Get the region of interest with respect to original geometry. 
		 If UseFractionROISpecification is true, a region of interest is generated
		 using the top,left,height,width fraction and any user specified 
		 RegionOfInterest will be overridden */
		itkSetMacro( RegionOfInterest, RegionType );
		itkGetMacro( RegionOfInterest, RegionType );
		
		/** Set/Get the region of interest for decompression. 
		 * The region of interest is defined by 4 floating point numbers:
		 * RegionOfInterestTop, RegionOfInterestLeft, RegionOfInterestHeight
		 * and RegionOfInterestWidth. These number represent the fraction of
		 * the full resolution image. For example, if top=0.3, left=0.2,
		 * height=0.6 and width=0.4 then the extracted region starts
		 * 30% down and 20% in from the left, and extends for 60%
		 * of the original height and 40% of the original width. */
		itkSetMacro( RegionOfInterestTop, double );
		itkSetMacro( RegionOfInterestLeft, double );
		itkSetMacro( RegionOfInterestHeight, double );
		itkSetMacro( RegionOfInterestWidth, double );
		itkGetMacro( RegionOfInterestTop, double );
		itkGetMacro( RegionOfInterestLeft, double );
		itkGetMacro( RegionOfInterestHeight, double );
		itkGetMacro( RegionOfInterestWidth, double );
		
		/** Get the full resolution image region. The output is only valid after a
		 * call to ReadImageInformation. */
		itkGetMacro( FullResolutionImageRegion, RegionType );
		
		/*-------- This part of the interfaces deals with writing data. ----- */
		
		/** Determine the file type. Returns true if this ImageIO can write the
		 * file specified. */
		virtual bool CanWriteFile(const char*);
		
		/** Set the spacing and dimension information for the set filename. */
		virtual void WriteImageInformation();
		
		/** Writes the data to disk from the memory buffer provided. Make sure
		 * that the IORegions has been set properly. */
		virtual void Write(const void* buffer);
		
		/** Set compression parameter string. This string will be split into blank delimited 
		 * substrings. Each substrings will be passed directly to the compression algorithm
		 * via the Kakadu siz_params. See Kakadu documentation for a description
		 * and syntax of allowable parameters. */
		itkSetStringMacro( CompressionParameterString );
		itkGetStringMacro( CompressionParameterString );
		
		/** Set the rate specificiation. This string will be parsed to
		 * determine the compression rate. Use the value UseAllBits to
		 * specify the final layers should include all bits.
		 *
		 * Example 1: Rates = [1.0,0.5,0.25] specifies irreversible compression to a 3 layer
		 * code-stream with three embedded rate.
		 *
		 * Example 2: rates = [UseAllBits,1.0,0.5,0.25] and
		 * CompressionParameterString = "Creversible=yes" 
		 * specifies reversible (lossless) compression with a progressive lossy to
		 * lossless code-stream with 4 layers. The use of "UseAllBits" specify that
		 * the final layer should include all remaining uncompressed bits, 
		 * not included in previous layers.
		 *
		 * For more detail see Kakadu documentation in kdu_compress.cpp.
		 */
		typedef Array<double> RatesType;
		itkSetMacro( Rates, RatesType );
		itkGetMacro( Rates, RatesType );
		itkGetMacro( UseAllBits, double );
		
		/** Set/Get the number of threads to use. */
		itkGetMacro( NumThreads, int );
		itkSetMacro( NumThreads, int );

		/** Read in a single component only by setting this index */
		itkGetMacro( ReadSingleComponentIndex, int );
		itkSetMacro( ReadSingleComponentIndex, int );

		
		JP2ImageIO();
		~JP2ImageIO();
		void PrintSelf(std::ostream& os, Indent indent) const;
		
	private:
		JP2ImageIO(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
		bool          m_Transpose;
		bool          m_VerticalFlip;
		bool          m_HorizontalFlip;
		bool          m_KeepColorChannelsSeparate;
		int           m_Reduce;
		bool          m_UseFractionROISpecification;
		double        m_RegionOfInterestTop;
		double        m_RegionOfInterestLeft;
		double        m_RegionOfInterestHeight;
		double        m_RegionOfInterestWidth;
		RegionType    m_RegionOfInterest;
		int           m_ReadSingleComponentIndex;
		RegionType    m_FullResolutionImageRegion;
		
		std::string   m_CompressionParameterString;
		RatesType     m_Rates;
		
		int           m_NumThreads;
		
		virtual void CreateThreads(kdu_thread_env& thread_env);
		virtual void CleanupThreads(kdu_thread_env& thread_env);
		
		static const double  m_UseAllBits;  
	};
	
} // end namespace itk

#endif // __itkJP2ImageIO_h
