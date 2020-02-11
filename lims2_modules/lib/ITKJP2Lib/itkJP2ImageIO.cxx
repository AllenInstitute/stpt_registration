/*=========================================================================

  itkJP2ImageIO.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/

#include "itkJP2ImageIO.h"
#include "itkExceptionObject.h"
#include "itkByteSwapper.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include <iostream>
#include <list>
#include <string>
#include <float.h>
#include <math.h>

#include "jp2.h"
#include "kdu_stripe_decompressor.h"
#include "kdu_stripe_compressor.h"

#include <algorithm>

static const int ZOOMIFY_TILE_SIZE = 256;


namespace itk
{
	
	const double JP2ImageIO::m_UseAllBits = NumericTraits<double>::max();
	
	/** Constructor */
	JP2ImageIO
	::JP2ImageIO()
	{
		m_ByteOrder = BigEndian;
		this->SetNumberOfDimensions(2);
		m_PixelType = SCALAR;
		m_ComponentType = UCHAR;
		m_Spacing[0] = 1.0;
		m_Spacing[1] = 1.0;
		
		m_Origin[0] = 0.0;
		m_Origin[1] = 0.0;
		
		m_Transpose      = false;
		m_VerticalFlip   = false;
		m_HorizontalFlip = false;
		m_KeepColorChannelsSeparate = false;
		m_Reduce         = 0;
		
		m_UseFractionROISpecification = false;
		
		int large = NumericTraits<int>::max();
		
		m_RegionOfInterest.SetIndex( 0, 0 );
		m_RegionOfInterest.SetIndex( 1, 0 );
		m_RegionOfInterest.SetSize( 0, large );
		m_RegionOfInterest.SetSize( 1, large );
		
		m_RegionOfInterestTop      = 0.0;
		m_RegionOfInterestLeft     = 0.0;
		m_RegionOfInterestHeight   = 1.0;
		m_RegionOfInterestWidth    = 1.0;
		
		/*
		 // Lossless defaults
		 m_CompressionParameterString = "Clayers=20 Creversible=yes Clevels=8 Cprecincts={256,256},{256,256},{128,128} Corder=RPCL ORGgen_plt=yes ORGtparts=R Cblk={32,32}";
		 m_Rates.SetSize( 2 );
		 m_Rates[0] = GetUseAllBits();
		 m_Rates[1] = 1.5;
		 */
		
		/*
		 // Kakadu defaults
		 m_CompressionParameterString = "Clayers=20 Creversible=yes Clevels=8 Cprecincts={256,256},{256,256},{128,128} Corder=RPCL ORGgen_plt=yes ORGtparts=R Cblk={32,32}";
		 m_Rates.SetSize( 0 );
		 */
		
		m_CompressionParameterString = "Clayers=20 Creversible=yes Clevels=8 Cprecincts={256,256},{256,256},{128,128} Corder=RPCL ORGgen_plt=yes ORGtparts=R Cblk={32,32}";
		//m_Rates.SetSize(2);
		//m_Rates[0] = DBL_MAX;
		//m_Rates[1] = 1.5;
		m_Rates.SetSize(1);
		m_Rates[0] = 1.5;
		
		m_NumThreads = 6;

		m_ReadSingleComponentIndex = -1;
	}
	
	
	/** Destructor */
	JP2ImageIO
	::~JP2ImageIO()
	{
	}
	
	
	bool JP2ImageIO
	::CanReadFile( const char* filename ) 
	{ 
		// First check the filename extension
		std::string fname = filename;
		if ( fname == "" )
		{
			itkDebugMacro(<< "No filename specified.");
		}
		
		bool extensionFound = false;
		std::string::size_type JP2Pos = fname.rfind(".jp2");
		if ((JP2Pos != std::string::npos)
			&& (JP2Pos == fname.length() - 4))
		{
			extensionFound = true;
		}
		
		JP2Pos = fname.rfind(".JP2");
		if ((JP2Pos != std::string::npos)
			&& (JP2Pos == fname.length() - 4))
		{
			extensionFound = true;
		}
		
		if( !extensionFound )
		{
			itkDebugMacro(<<"The filename extension is not recognized");
			return false;
		}
		
		// Now check the content
		jp2_source jp2_in;
		jp2_family_src jp2_ultimate_src;
		
		jp2_ultimate_src.open( filename );
		if ( !jp2_in.open( &jp2_ultimate_src ) )
		{
			itkDebugMacro(<<"File: " << filename <<
						  "does not appear to contain any boxes.");
			jp2_in.close();
			return false;
		}
		
		try
		{
			jp2_in.read_header();
		}
		catch( ExceptionObject & err )
		{
			itkDebugMacro( << err.GetDescription() );
			jp2_in.close();
			return false;
		}
		
		kdu_codestream codestream;
		codestream.create( &jp2_in );
		
		int depth = codestream.get_bit_depth( 0 );
		switch ( depth )
		{
			case 8:
				m_ComponentType = UCHAR;
				break;
			case 16:
				m_ComponentType = USHORT;
				break;
			default:
				itkDebugMacro( << "Can only support 8-bit or 16-bit components" );
				codestream.destroy();
				jp2_in.close();
				return false;
		}
		
		codestream.destroy();
		jp2_in.close();
		return true;		
	}
	
	
	
	bool JP2ImageIO
	::CanWriteFile( const char * name )
	{		
		std::string filename = name;
		if ( filename == "" )
		{
			itkDebugMacro(<< "No filename specified.");
		}
		
		bool extensionFound = false;
		std::string::size_type JP2Pos = filename.rfind(".jp2");
		if ((JP2Pos != std::string::npos)
			&& (JP2Pos == filename.length() - 4))
		{
			extensionFound = true;
		}
		
		JP2Pos = filename.rfind(".JP2");
		if ((JP2Pos != std::string::npos)
			&& (JP2Pos == filename.length() - 4))
		{
			extensionFound = true;
		}
		
		if( !extensionFound )
		{
			itkDebugMacro(<<"The filename extension is not recognized");
			return false;
		}
		
		if( extensionFound )
		{
			return true;
		}
		
		return false;		
	}
	
	
	
	void JP2ImageIO
	::Read(void* buffer)
	{
		// Create a code stream
		jp2_source jp2_in;
		jp2_family_src jp2_ultimate_src;
		
		jp2_ultimate_src.open( this->GetFileName() );
		jp2_in.open( &jp2_ultimate_src );
		jp2_in.read_header();
		
		kdu_codestream codestream;
		codestream.create( &jp2_in );
		codestream.set_fussy();
		
		// Reconstruct region of interest
		kdu_dims region;
		region.pos.x  = m_RegionOfInterest.GetIndex( 0 );
		region.pos.y  = m_RegionOfInterest.GetIndex( 1 );
		region.size.x = m_RegionOfInterest.GetSize( 0 );
		region.size.y = m_RegionOfInterest.GetSize( 1 );

		int requestedComponent = this->GetReadSingleComponentIndex();
		int startComponent = (requestedComponent >= 0) ? requestedComponent : 0;
		
		// Do apply_input_restructions
		codestream.apply_input_restrictions( startComponent, this->GetNumberOfComponents(), m_Reduce, 0, &region );
		
		// Do change_appearance
		codestream.change_appearance( m_Transpose, m_VerticalFlip, m_HorizontalFlip );
		
		// Get image meta info
		kdu_dims dims;
		codestream.get_dims( 0, dims );
		
		// Create threads
		// kdu_thread_env thread_env;
		// CreateThreads( thread_env );
		
		// Decompress image in one hit
		kdu_stripe_decompressor decompressor;
		// decompressor.start( codestream, false, false, &thread_env );
		decompressor.start( codestream );
		int stripe_heights[3] = { dims.size.y, dims.size.y, dims.size.y };
		
		if ( this->GetComponentType() == UCHAR )
		{
			kdu_byte *p = static_cast<kdu_byte *>(buffer);
			
			decompressor.pull_stripe( p, stripe_heights );
		}
		else if ( this->GetComponentType() == USHORT )
		{
			kdu_int16 *p = static_cast<kdu_int16 *>(buffer);
			bool is_signed[3] = { false, false, false };
			
			decompressor.pull_stripe( p, stripe_heights, NULL, NULL, NULL, NULL, is_signed );
		}
		
		decompressor.finish();
		
		// Clean up
		codestream.destroy();
		jp2_in.close();
		
		// CleanupThreads(thread_env);
	}
	
	
	/** 
	 *  Read Information about the JPEG file
	 *  and put the cursor of the stream just before the first data pixel
	 */
	void JP2ImageIO
	::ReadImageInformation()
	{
		
		// Create a code stream
		jp2_source jp2_in;
		jp2_family_src jp2_ultimate_src;
		
		jp2_ultimate_src.open( this->GetFileName() );
		jp2_in.open( &jp2_ultimate_src );
		jp2_in.read_header();
		
		kdu_codestream codestream;
		codestream.create( &jp2_in );
		
		// Get number of components
		int num_components = codestream.get_num_components();
		int requestedComponent = this->GetReadSingleComponentIndex();

		if (requestedComponent >= num_components)
		{
		    itkExceptionMacro( << "Requested component (" << requestedComponent << ") is out of bounds." );
		}

		if ( num_components == 1 )
		{
			this->SetNumberOfComponents( 1 );
			this->SetPixelType( SCALAR );
			this->SetReadSingleComponentIndex( 0 );
		}
		else if ( num_components == 3 )
		{
		    if ( requestedComponent >= 0 )
		    {
			// There are multiple channels, but a single channel was requested.
			num_components = 1;
			this->SetNumberOfComponents( 1 );			
			this->SetPixelType( SCALAR );
		    }
		    else
		    {
			// Multiple channels, but no channel was requested.
			this->SetNumberOfComponents( 3 );
			this->SetPixelType( RGB );
		    }
		}
		
		// Get bit depth
		int depth = codestream.get_bit_depth( 0 );
		switch ( depth )
		{
			case 8:
				m_ComponentType = UCHAR;
				break;
			case 16:
				m_ComponentType = USHORT;
				break;
			default:
				break;
		}
		
		// Compute region of interest
		kdu_dims region;
		siz_params * siz = codestream.access_siz();
		
		if (!(siz->get(Sorigin,0,0,region.pos.y) &&
			  siz->get(Sorigin,0,1,region.pos.x) &&
			  siz->get(Ssize,0,0,region.size.y) &&
			  siz->get(Ssize,0,1,region.size.x)))
		{
			itkExceptionMacro( << "Error with getting image size" );
		}
		
		region.size.y -= region.pos.y;
		region.size.x -= region.pos.x;
		
		m_FullResolutionImageRegion.SetIndex( 0, region.pos.x );
		m_FullResolutionImageRegion.SetIndex( 1, region.pos.y );
		m_FullResolutionImageRegion.SetSize( 0, region.size.x );
		m_FullResolutionImageRegion.SetSize( 1, region.size.y );
		
		if ( m_UseFractionROISpecification )
		{
			if ( ( m_RegionOfInterestTop < 0.0 ) ||
				( m_RegionOfInterestLeft < 0.0 ) ||
				( m_RegionOfInterestHeight < 0.0 ) ||
				( m_RegionOfInterestWidth < 0.0 ) )
			{
				itkExceptionMacro( << "ROI specification must non-zero" );
			} 
			
			
			region.pos.y += (int) floor( region.size.y * m_RegionOfInterestTop );
			region.pos.x += (int) floor( region.size.x * m_RegionOfInterestLeft );
			region.size.y = (int) ceil( region.size.y  * m_RegionOfInterestHeight );
			region.size.x = (int) ceil( region.size.x  * m_RegionOfInterestWidth );
			
			m_RegionOfInterest.SetIndex( 0, region.pos.x );
			m_RegionOfInterest.SetIndex( 1, region.pos.y );
			m_RegionOfInterest.SetSize( 0, region.size.x );
			m_RegionOfInterest.SetSize( 1, region.size.y );
			
		}
		else
		{
			region.pos.x  = m_RegionOfInterest.GetIndex( 0 );
			region.pos.y  = m_RegionOfInterest.GetIndex( 1 );
			region.size.x = m_RegionOfInterest.GetSize( 0 );
			region.size.y = m_RegionOfInterest.GetSize( 1 );
		}
		
		
		// Do change_appearance
		codestream.change_appearance( m_Transpose, m_VerticalFlip, m_HorizontalFlip );
		
		// Do apply_input_restructions
		codestream.apply_input_restrictions( 0, 0, m_Reduce, 0, &region );
		
		// Get image meta info
		kdu_dims dims;
		codestream.get_dims( 0, dims );
		
		this->SetNumberOfDimensions(2);
		m_Dimensions[0] = dims.size.x;
		m_Dimensions[1] = dims.size.y;
		
		// Check that components have consistent dimension
		if ( num_components == 3 )
		{
			kdu_dims dims1, dims2;
			codestream.get_dims( 1, dims1 );
			codestream.get_dims( 2, dims2 );
			if ( (dims1 != dims) || (dims2 != dims) )
			{
				num_components = 1;
				this->SetReadSingleComponentIndex( 0 );
				this->SetNumberOfComponents( 1 );
				this->SetPixelType( SCALAR );
			}
		}
		
		codestream.destroy();
		jp2_in.close();
		
	}
	
	
	void 
	JP2ImageIO
	::WriteImageInformation(void)
	{		
	}
	
	
	void 
	JP2ImageIO
	::Write( const void* buffer) 
	{
		int bitDepth;
		if ( this->GetComponentType() == UCHAR )
			bitDepth = 8;
		else if ( this->GetComponentType() == USHORT )
			bitDepth = 16;
		else
			itkExceptionMacro(<<"JP2ImageIO supports unsigned char or unsigned short components only");

		if ( this->GetNumberOfComponents() != 1 &&
			this->GetNumberOfComponents() != 3 )
		{
			itkExceptionMacro(<<"JP2ImageIO supports 1 or 3 components only" );
		}
		
		if ( this->GetNumberOfDimensions() != 2 )
			itkExceptionMacro(<<"JP2ImageIO supports 2D images only" );
		
		if ( m_Rates.size() == 0 )
			itkExceptionMacro(<<"JP2ImageIO error -- no target bit rate specified" );
		
		// Set up code stream siz parameter
		siz_params siz;
		siz.set( Scomponents, 0, 0, (int) this->GetNumberOfComponents() );
		siz.set( Sdims, 0, 0, (int) m_Dimensions[1] );  // height
		siz.set( Sdims, 0, 1, (int) m_Dimensions[0] );  // width
		siz.set( Sprecision, 0, 0, bitDepth );  // component bit depth
		siz.set( Ssigned, 0, 0, false );  // unsigned data
				
		siz_params *siz_ref = &siz;
		reinterpret_cast<kdu_params *>(siz_ref)->finalize();
		
		siz_params siz_scratch;
		if (m_Transpose)
		{
			siz_scratch.copy_from(siz_ref, -1, -1, -1, 0, 0, m_Transpose, false, false);
			siz_ref = &siz_scratch;
		}
		
		// Construct a target object and codestream
		jp2_target output;
		jp2_family_tgt ultimate_target;
		ultimate_target.open( this->GetFileName() );
		output.open( &ultimate_target );
		
		kdu_codestream codestream;
		codestream.create(siz_ref, &output);
				
		// Make sure the number of tiers/levels is what is expected by ZoomifyEncoder
		// Zoomify expects the smallest tier's dimensions to be <= the tile size
		double tiers = std::max(::log(static_cast<double>(m_Dimensions[0]) / static_cast<double>(ZOOMIFY_TILE_SIZE)) / ::log(2.0),
								::log(static_cast<double>(m_Dimensions[1]) / static_cast<double>(ZOOMIFY_TILE_SIZE)) / ::log(2.0));
		int numTiers = static_cast<int>(tiers);
		numTiers++;             // Add the full resolution tier
		
		double intPart;
		double fracPart = ::modf(tiers, &intPart);
		if (fracPart != 0.0)
			numTiers++;
		
		std::stringstream compressionParams;
		compressionParams << "Clayers=20 Creversible=yes Clevels=" << numTiers << " Cycc=" << (m_KeepColorChannelsSeparate ? "no" : "yes")
						  << " Cprecincts={256,256},{256,256},{128,128} Corder=RPCL ORGgen_plt=yes ORGtparts=R Cblk={32,32}";
		compressionParams.seekg(0);		
		m_CompressionParameterString = compressionParams.str();
		
		// Break up parameter string into substrings and pass them to the encoder
		std::istringstream instream;
		instream.str( m_CompressionParameterString );
		while( 1 )
		{
			std::string substring;
			instream >> substring;
			codestream.access_siz()->parse_string( substring.c_str() );
			if ( substring.length() == 0 ) 
			{
				break;
			}
		}
				
		// Set the compression rate as specify by user
		std::vector<kdu_long> layerBytes;
		if ( m_Rates.GetSize() > 0 )
		{			
			int argSpecs = m_Rates.GetSize();
			
			int codSpecs; 
			kdu_params *cod = codestream.access_siz()->access_cluster( COD_params );
			cod->get( Clayers, 0, 0, codSpecs );
			
			// Deterimine the number of required layers
			int numSpecs = codSpecs;
			if ( numSpecs == 0 )
			{
				numSpecs = 1;
			}
			
			if ( ( argSpecs != numSpecs ) &&
				( ( argSpecs > 2 ) ||  
				 ( ( argSpecs == 2 ) && ( numSpecs == 1 ) ) ) )
			{
				itkExceptionMacro( << "The number of bit rates and the number of quality layers "
								  << "does not conform to the rules. " );
			}
			
			cod->set( Clayers, 0, 0, numSpecs );
			
			layerBytes.resize( numSpecs, 0 );
			
			kdu_long totalPixels = static_cast<kdu_long>(m_Dimensions[0]) * static_cast<kdu_long>(m_Dimensions[1]);
			bool useAllBits = false;
			
			for ( int k = 0; k < argSpecs; k++ )
			{
				if( m_Rates[k] == m_UseAllBits )
				{
					useAllBits = true;
					layerBytes[k] = NumericTraits<signed long>::max();
				}
				else
				{
					layerBytes[k] = static_cast<kdu_long>( floor( m_Rates[k] * static_cast<double>(totalPixels) / 8.0 ) );
				}
			}
			
			// sort in ascending order
			std::sort( layerBytes.begin(), layerBytes.begin() + argSpecs );
			
			if ( argSpecs && ( argSpecs != numSpecs ) )
			{
				layerBytes[numSpecs - 1] = layerBytes[argSpecs - 1];
				layerBytes[argSpecs - 1] = 0;
			}
			
			if ( useAllBits )
			{
				// force rate to allocate all remaining compressed bits to last layer
				layerBytes[numSpecs - 1] = 0;
			}
			
			if ( ( numSpecs > 0 ) && ( layerBytes[numSpecs - 1] > 0 ) )
			{
				codestream.set_max_bytes( layerBytes[numSpecs - 1] );
			}
		}
				
		// Finalize 
		codestream.access_siz()->finalize_all(); // Set up coding defaults
		
		if ( ultimate_target.exists() )
		{
			// Do JP2 file initialization
			jp2_dimensions dimensions = output.access_dimensions();
			dimensions.init( codestream.access_siz() );
			jp2_colour colour = output.access_colour();
			colour.init((this->GetNumberOfComponents() >= 3 ) ? JP2_sRGB_SPACE : JP2_sLUM_SPACE);
			output.write_header();
			output.open_codestream(true);
		}
		else
		{
			itkExceptionMacro(<< "Error while creating codestream object" );
		}
		
		// Do change_appearance
		codestream.change_appearance( m_Transpose, m_VerticalFlip, m_HorizontalFlip );
		
		// *** Activating threading currently (Kakadu 6.4.1) causes the compressor to ignore the target bit rate
		// kdu_thread_env thread_env;
		// CreateThreads(thread_env);
		
		// Create compressor and run
		kdu_stripe_compressor compressor;
		/*
		// if (maxRate == m_UseAllBits)
		if (!layerBytes.size())
			compressor.start(codestream, 0, NULL, NULL, 0,
							 false, false, true, 0.0, 0,
							 false, &thread_env, NULL, 0);

		else
			compressor.start(codestream, layerBytes.size(), &layerBytes[0], NULL, 0,
							 false, false, true, 0.0, 0,
							 false, &thread_env, NULL, 0);
		*/
		
		compressor.start(codestream, 0, NULL, NULL, 0,
						 false, false, true, 0.0, 0,
						 false, NULL, NULL, 0);
		
		void *tempPtr = const_cast<void *>( buffer );
		
		int height = (int) m_Dimensions[1];
		int stripe_heights[3] = { height, height, height };
		
		if (bitDepth == 8)
		{
			kdu_byte *p = static_cast<kdu_byte *>( tempPtr );
			compressor.push_stripe(p, stripe_heights);
		}
		else if (bitDepth == 16)
		{
			kdu_int16 *p = static_cast<kdu_int16 *>( tempPtr );		
			int precisions[3] = { bitDepth, bitDepth, bitDepth };
			bool is_signed[3] = { false, false, false };
			compressor.push_stripe(p, stripe_heights, NULL, NULL, NULL, precisions, is_signed, 0);
		}
		
		compressor.finish();
		
		// Clean up
		codestream.destroy();
		output.close();
		
		// CleanupThreads(thread_env);
	}
	
	void JP2ImageIO
	::CreateThreads(kdu_thread_env& thread_env)
	{
		thread_env.create();
		for (int nt = 0; nt < m_NumThreads; nt++)
			if (!thread_env.add_thread())
				break;
	}
	
	void JP2ImageIO
	::CleanupThreads(kdu_thread_env& thread_env)
	{
		thread_env.terminate(NULL, true); // Makes sure any outstanding flush job finishes
		thread_env.destroy();
	}
	
	/** Print Self Method */
	void JP2ImageIO
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Transpose: " << m_Transpose << std::endl;
		os << indent << "VerticalFlip: " << m_VerticalFlip << std::endl;
		os << indent << "HorizontalFlip: " << m_HorizontalFlip << std::endl;
		os << indent << "Reduce: " << m_Reduce << std::endl;
		os << indent << "UseFractionROISpecification: " << m_UseFractionROISpecification << std::endl;
		os << indent << "RegionOfInterest: " << m_RegionOfInterest << std::endl;
		os << indent << "RegionOfInterestTop: " << m_RegionOfInterestTop << std::endl;
		os << indent << "RegionOfInterestLeft: " << m_RegionOfInterestLeft << std::endl;
		os << indent << "RegionOfInterestHeight: " << m_RegionOfInterestHeight << std::endl;
		os << indent << "RegionOfInterestWidth: " << m_RegionOfInterestWidth << std::endl;
		os << indent << "FullResolutionImageRegion: " << m_FullResolutionImageRegion << std::endl;
		os << indent << "Rates: " << m_Rates << std::endl;
		
	}
	
} // end namespace itk
