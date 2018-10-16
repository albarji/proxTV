#!/usr/bin/env bash
set -v

# Deploys built wheels to PyPI
if [[ ! -z "$TRAVIS_TAG" ]]
then
    pip install twine
    ls -l dist/*  # List files to upload
    twine upload --username albarji dist/*
fi
