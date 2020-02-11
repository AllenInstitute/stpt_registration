PROJECT(Register)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


# Find ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Find VTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Include TinyXML
SET(TINY_XML_PATH CACHE PATH "Top directory of your TinyXML installation")
INCLUDE(${TINY_XML_PATH}/UseTinyXml.cmake)

# Include Boost
SET (BOOST_DIR CACHE PATH "Boost library path")
SET (BOOST_LIBRARIES_USED boost_filesystem boost_program_options boost_system)

# Set include and link directories for boost libraries
INCLUDE_DIRECTORIES(${BOOST_DIR}/include)
LINK_DIRECTORIES(${BOOST_DIR}/lib)

# ITK Kakadu wrapper library
# Check if it's been pulled down into the source directory for build
# or if it's been prebuilt somewhere else.
SET(ITK_JP2_SRC_DIR ${Register_SOURCE_DIR}/ITKJP2Lib)
IF(IS_DIRECTORY ${ITK_JP2_SRC_DIR})
    # Add directory for build
    ADD_SUBDIRECTORY(${ITK_JP2_SRC_DIR})
    INCLUDE(${ITKJP2Lib_BINARY_DIR}/UseITKJP2.cmake)
ELSE(IS_DIRECTORY ${ITK_JP2_SRC_DIR})
    # Get location from the user
    SET(ITK_JP2_DIR CACHE PATH "ITK Kakadu wrapper library")

    INCLUDE(${ITK_JP2_DIR}/UseITKJP2.cmake)
ENDIF(IS_DIRECTORY ${ITK_JP2_SRC_DIR})


# Set output directories
SET (EXECUTABLE_OUTPUT_PATH ${Register_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")


# Build Options
OPTION( BUILD_UTILITIES "Build command line utilities directory." OFF )
OPTION( BUILD_REGISTRATION "Build registration directory." ON )


# Command line applications
IF ( BUILD_UTILITIES )
	SUBDIRS( Utilities )
ENDIF ( BUILD_UTILITIES )

IF ( BUILD_REGISTRATION )
	SUBDIRS( Registration )
ENDIF ( BUILD_REGISTRATION )

