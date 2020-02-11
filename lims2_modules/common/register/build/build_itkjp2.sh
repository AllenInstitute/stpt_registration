#!/bin/bash
#
# Build the ITKJP2 library
#

set -o errexit

pushd $1 > /dev/null

bash build/configure.sh $2
cd $2
make ITKJP2

popd > /dev/null

