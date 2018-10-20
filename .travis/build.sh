#!/usr/bin/env bash
set -v

# Get just the first two numbers of python version, without dots
python --version
PYTHONVERSION=$(echo ${TRAVIS_PYTHON_VERSION} | tr -d .)
echo ${PYTHONVERSION}
# Linux build and install
if [ "${TRAVIS_OS_NAME}" == "linux" ]
then
    docker pull quay.io/pypa/manylinux1_x86_64
    PACKAGING=$([[ ${PYTHONVERSION} = "27" ]] && echo "mu" || echo "m")
    # Build precompiled python wheels
    docker run --rm -it -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 bash -c "cd /io; bash buildwheels.sh ${PYTHONVERSION} ${PACKAGING}"
    # Copy to dist folder
    mkdir dist
    cp wheelhouse/prox_tv-*-cp${PYTHONVERSION}-cp${PYTHONVERSION}${PACKAGING}-manylinux1_x86_64.whl dist
    # Install wheel for this with version
    pip install wheelhouse/prox_tv-*-cp${PYTHONVERSION}-cp${PYTHONVERSION}${PACKAGING}-manylinux1_x86_64.whl

# Mac build and install
elif [ "${TRAVIS_OS_NAME}" == "osx" ]
then
    # Build wheel
    pip wheel . -w wheelhouse/
    # Bundle dependencies
    pip install delocate
    delocate-listdeps wheelhouse/*.whl
    # Copy to dist folder
    mkdir dist
    cp wheelhouse/prox_tv-*-cp*-cp*-macosx_*.whl dist
    # Install wheel for this with version
    pip install wheelhouse/prox_tv-*-cp${PYTHONVERSION}-cp${PYTHONVERSION}m-macosx_*.whl
fi
