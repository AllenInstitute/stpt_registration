#!/bin/bash

if [[ ${MACHTYPE} =~ 'apple' ]]; then
	CMAKE=cmake
	CMAKE_ARCH="-DCMAKE_OSX_ARCHITECTURES:STRING=x86_64"
else
	# Set CMake and library locations
	CMAKE=${CMAKE:-/shared/bioapps/cmake-latest/bin/cmake}
	KAKADU_DIR=${KAKADU_DIR:-/shared/bioapps/kakadu/latest}

	ITK_DIR=${ITK_DIR:-/shared/bioapps/itk/latest}
	export ITK_DIR
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
echo "Using kakadu: " ${KAKADU_DIR}
echo "Using ITK: " ${ITK_DIR}

# Create build directory and run CMake to generate the Makefile
if [[ ! -d ${BIN_DIR} ]]; then
	mkdir ${BIN_DIR}
fi

pushd ${BIN_DIR} > /dev/null

# The second argument to the script can set to the build system
${CMAKE} ${CMAKE_ARCH} -DCMAKE_BUILD_TYPE:STRING=Release \
  -DKAKADU_PATH:PATH=${KAKADU_DIR} \
  $2 ${SOURCE_DIR}

# Return the CMake exit status
status=$?

popd > /dev/null

exit ${status}
