#!/usr/bin/env bash
set -v

# Linux build and install
if [ "${TRAVIS_OS_NAME}" == "linux" ]
then
    docker pull quay.io/pypa/manylinux1_x86_64
    PYTHONVERSION=$(echo ${TRAVIS_PYTHON_VERSION} | tr -d .)
    PACKAGING=$([[ ${PYTHONVERSION} = "27" ]] && echo "mu" || echo "m")
    # Build precompiled python wheels
    docker run --rm -it -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 bash -c "cd /io; bash buildwheels.sh"
    # Install wheel for this with version
    pip install wheelhouse/prox_tv-*-cp${PYTHONVERSION}-cp${PYTHONVERSION}${PACKAGING}-manylinux1_x86_64.whl

# Mac build and install
elif [ "${TRAVIS_OS_NAME}" == "osx" ]
then
    #brew install gcc  # TODO: needed?
    pip wheel . -w wheelhouse/
    pip install delocate
    delocate-listdeps wheelhouse/*.whl
    pip install wheelhouse/prox_tv-*-cp*-cp*-macosx_*.whl
fi
