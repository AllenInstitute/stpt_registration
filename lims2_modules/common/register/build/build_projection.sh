#!/bin/bash
#
# Runs the configure script if the bin directory does not exist.
# Runs make in the bin directory.
#

set -o errexit

if [[ -z "$1" ]]; then
	BIN_DIR="bin"
else
	BIN_DIR=$1
fi

# Build ITK JP2 library
ITK_JP2_SRC=../../lib/ITKJP2Lib
ITK_JP2_BIN=${ITK_JP2_SRC}/bin
bash build/build_itkjp2.sh ${ITK_JP2_SRC} ${ITK_JP2_BIN}
export ITK_JP2_DIR=${ITK_JP2_BIN}

if [[ ! -d ${BIN_DIR} || ! -f ${BIN_DIR}/Makefile ]]; then
	bash build/configure.sh ${BIN_DIR}
fi

pushd ${BIN_DIR} > /dev/null

make idpProjectionAlignmentModule idpProjectionLocalAlignmentModule idpProjectionResampleVolumeRGBModule 

popd > /dev/null

