INCLUDE_DIRECTORIES(
${Register_SOURCE_DIR}/Common
${Register_SOURCE_DIR}/Registration
${Register_SOURCE_DIR}/Grid
)

SET (ITK_LIBRARIES_USED ITKCommon ITKIO ITKNumerics)

ADD_EXECUTABLE( DownsampleImage DownsampleImage.cxx )
TARGET_LINK_LIBRARIES( DownsampleImage ${ITK_LIBRARIES_USED} ${ITK_JP2_LIBRARIES})

