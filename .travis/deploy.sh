#!/usr/bin/env bash
set -v

# Select deploy repo
if [[ ! -z "$TRAVIS_TAG" ]]; then
    repo=""
else
    repo="--repository-url https://test.pypi.org/legacy/"
fi

# Deploys built wheels to PyPI
if [[ ! -z "$TRAVIS_TAG" ]]
then
    pip install twine
    ls -l dist/*  # List files to upload
    twine upload ${repo} --username albarji dist/*
fi
