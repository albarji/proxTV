#!/usr/bin/env bash
# Source: https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
# Must be run inside container quay.io/pypa/manylinux1_x86_64

# Install system packages required by the library
yum install -y blas-devel lapack-devel

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/lib64"

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/ -vv
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/prox_tv*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done
