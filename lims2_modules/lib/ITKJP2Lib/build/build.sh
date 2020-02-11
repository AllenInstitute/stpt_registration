#!/bin/bash
#
# Runs the configure script if the bin directory does not exist.
# Runs make in the bin directory.
#

if [[ -z "$1" ]]; then
	BIN_DIR="bin"
else
	BIN_DIR=$1
fi

if [[ ! -d ${BIN_DIR} || ! -f ${BIN_DIR}/Makefile ]]; then
	bash build/configure.sh ${BIN_DIR}
fi

pushd ${BIN_DIR} > /dev/null

make ITKJP2

# Save the make exit status
status=$?

popd > /dev/null

exit ${status}
