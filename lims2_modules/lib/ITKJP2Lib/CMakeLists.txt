CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ITKJP2Lib)

FIND_PACKAGE(ITK)

IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)


# Kakadu sources
SET(KAKADU_PATH CACHE PATH "Version directory of your Kakadu library installation")

SET (KAKADU_SRC
  ${KAKADU_PATH}/apps/support/kdu_stripe_compressor.cpp
  ${KAKADU_PATH}/apps/support/kdu_region_decompressor.cpp
  ${KAKADU_PATH}/apps/jp2/jp2.cpp
  ${KAKADU_PATH}/apps/support/kdu_stripe_decompressor.cpp 
  itkJP2ImageIO.cxx
  itkJP2ImageIOFactory.cxx
  itkVaa3dRawImageIO.cxx
  itkVaa3dRawImageIOFactory.cxx
  itkExtendedImageIOFactory.cxx)

INCLUDE_DIRECTORIES(
  ${KAKADU_PATH}/coresys/common
  ${KAKADU_PATH}/apps/compressed_io
  ${KAKADU_PATH}/apps/support)

ADD_LIBRARY(ITKJP2 ${KAKADU_SRC})

IF (WIN32)
  LINK_DIRECTORIES(${KAKADU_PATH}/../lib)
  TARGET_LINK_LIBRARIES(ITKJP2 ITKCommon ITKIO kdu_v64R)
ELSE (WIN32)
  LINK_DIRECTORIES(${KAKADU_PATH}/lib)
  TARGET_LINK_LIBRARIES(ITKJP2 ITKCommon ITKIO kdu)
ENDIF (WIN32)


CONFIGURE_FILE(UseITKJP2.cmake.in
               ${ITKJP2Lib_BINARY_DIR}/UseITKJP2.cmake)


INCLUDE(${ITKJP2Lib_BINARY_DIR}/UseITKJP2.cmake)
ADD_EXECUTABLE(ImageReadWrite ImageReadWrite.cxx)
TARGET_LINK_LIBRARIES(ImageReadWrite ITKCommon ITKIO ${ITK_JP2_LIBRARIES})

ADD_EXECUTABLE(ReduceRGBJP2Image ReduceRGBJP2Image.cxx)
TARGET_LINK_LIBRARIES(ReduceRGBJP2Image ITKCommon ITKIO ${ITK_JP2_LIBRARIES})

ADD_EXECUTABLE(ConvertVaa3dRaw ConvertVaa3dRaw.cxx)
TARGET_LINK_LIBRARIES(ConvertVaa3dRaw ITKCommon ITKIO ${ITK_JP2_LIBRARIES})
