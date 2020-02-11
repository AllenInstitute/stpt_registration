#!/bin/bash

set -o errexit

if [[ ${MACHTYPE} =~ 'apple' ]]; then
	CMAKE=cmake
	CMAKE_ARCH="-DCMAKE_OSX_ARCHITECTURES:STRING=x86_64"
else
	# Set CMake and library locations
	CMAKE=${CMAKE:-/shared/bioapps/cmake-latest/bin/cmake}
	TINY_XML_DIR=${TINY_XML_DIR:-/shared/bioapps/tinyxml/tinyxml.x86_64}
	BOOST_DIR=${BOOST_DIR:-/shared/bioapps/boost/boost.x86_64}

	ITK_DIR=${ITK_DIR:-/shared/bioapps/itk/latest}
	export ITK_DIR

	VTK_DIR=${VTK_DIR:-/shared/bioapps/vtk/latest}
	export VTK_DIR
fi

# Current directory contains the source code
SOURCE_DIR=`pwd`

# Set the build directory to the first argument if given
if [[ -z "$1" ]]
then
	BIN_DIR="bin"
else
	BIN_DIR=$1
fi

# Show directory paths for debugging
echo "Source directory: " ${SOURCE_DIR}
echo "Build directory: " ${BIN_DIR}
echo "Using ITK: " ${ITK_DIR}
echo "Using VTK: " ${VTK_DIR}
echo "Using boost: " ${BOOST_ROOT}
echo "Using tinyxml: " ${TINY_XML_DIR}

# Create build directory and run CMake to generate the Makefile
if [[ ! -d ${BIN_DIR} ]]; then
	mkdir ${BIN_DIR}
fi

pushd ${BIN_DIR} > /dev/null

# The second argument to the script can set to the build system
${CMAKE} ${CMAKE_ARCH} -DCMAKE_BUILD_TYPE:STRING=Release \
  -DITK_JP2_DIR:PATH=${ITK_JP2_DIR} \
  -DTINY_XML_PATH=${TINY_XML_DIR} \
  -DBOOST_DIR:PATH=${BOOST_DIR} \
  $2 ${SOURCE_DIR}

popd > /dev/null

