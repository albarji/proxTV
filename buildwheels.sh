#!/usr/bin/env bash
# Source: https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
# Must be run inside container quay.io/pypa/manylinux1_x86_64
set -v

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Parameters: <python version to use in build> <python packaging (m or mu)>"
    exit
fi
pythonversion=$1
packaging=$2

# Install system packages required by the library
yum install -y blas-devel lapack-devel

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/lib64"

# Compile wheel
"/opt/python/cp${pythonversion}-cp${pythonversion}${packaging}/bin/pip" wheel /io/ -w wheelhouse/

# Bundle external shared libraries into the wheel
auditwheel repair wheelhouse/prox_tv-*-cp${pythonversion}-cp${pythonversion}${packaging}-linux_x86_64.whl -w /io/wheelhouse/
