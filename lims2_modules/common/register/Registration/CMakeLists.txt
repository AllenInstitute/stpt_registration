INCLUDE_DIRECTORIES(${Register_SOURCE_DIR}/Common)

SET (ITK_LIBRARIES_USED ITKCommon ITKIO ITKNumerics)

ADD_EXECUTABLE( idpProjectionAlignmentModule
idpProjectionAlignmentModule.cxx
${Register_SOURCE_DIR}/Common/idpDecomposeAffineMatrix2D.cxx
${Register_SOURCE_DIR}/Common/idpDecomposeAffineMatrix3D.cxx
${Register_SOURCE_DIR}/Common/idpImageSeries.cxx 
${Register_SOURCE_DIR}/Common/idpRegistrationUtilities.cxx
${Register_SOURCE_DIR}/Common/idpTransformUtilities.cxx
${Register_SOURCE_DIR}/Common/idpReferenceSpaceUtilities.cxx
)
TARGET_LINK_LIBRARIES( idpProjectionAlignmentModule ${ITK_LIBRARIES_USED} ${ITK_JP2_LIBRARIES} ${TINY_XML_LIBRARIES})


ADD_EXECUTABLE( idpProjectionLocalAlignmentModule
idpProjectionLocalAlignmentModule.cxx
${Register_SOURCE_DIR}/Common/idpDecomposeAffineMatrix2D.cxx
${Register_SOURCE_DIR}/Common/idpDecomposeAffineMatrix3D.cxx
${Register_SOURCE_DIR}/Common/idpImageSeries.cxx 
${Register_SOURCE_DIR}/Common/idpRegistrationUtilities.cxx
${Register_SOURCE_DIR}/Common/idpTransformUtilities.cxx
${Register_SOURCE_DIR}/Common/idpReferenceSpaceUtilities.cxx
)
TARGET_LINK_LIBRARIES( idpProjectionLocalAlignmentModule ${ITK_LIBRARIES} ${ITK_JP2_LIBRARIES} ${TINY_XML_LIBRARIES})

ADD_EXECUTABLE( idpProjectionResampleVolumeRGBModule 
idpProjectionResampleVolumeRGBModule.cxx
${Register_SOURCE_DIR}/Common/idpDecomposeAffineMatrix3D.cxx
${Register_SOURCE_DIR}/Common/idpImageSeries.cxx 
${Register_SOURCE_DIR}/Common/idpRegistrationUtilities.cxx
${Register_SOURCE_DIR}/Common/idpTransformUtilities.cxx
${Register_SOURCE_DIR}/Common/idpReferenceSpaceUtilities.cxx
)
TARGET_LINK_LIBRARIES( idpProjectionResampleVolumeRGBModule ${ITK_LIBRARIES_USED} ${ITK_JP2_LIBRARIES} ${TINY_XML_LIBRARIES})
