#!/usr/bin/env bash
set -v

# Deploys built wheels to PyPI
pip install twine
ls -l dist/*  # List files to upload
if [[ ! -z "$TRAVIS_TAG" ]]
then
    twine upload --username albarji dist/*
else  # Allow failure in test PyPI
    twine upload --repository-url https://test.pypi.org/legacy/ --username albarji dist/* || :
fi
