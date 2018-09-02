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
    #apt-get install -y libblas-devel liblapack-devel
    #sudo launchctl load -w /System/Library/LaunchDaemons/com.apple.locate.plist
    #locate libgfortran.a
    brew install gcc  # TODO: needed?
    #locate libgfortran.a
    pip wheel . -w wheelhouse/
    pip install delocate
    delocate-listdeps wheelhouse/*.whl
    ls wheelhouse
    pip install wheelhouse/prox_tv-*-cp*-cp*-macosx_*.whl
    which python
    export PATH=$PATH:/usr/local/bin  # TODO Requirement already satisfied. Updating path needed?
fi
