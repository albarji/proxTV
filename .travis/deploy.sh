#!/usr/bin/env bash
set -v

# TODO: run only for tags
# TODO: use real PyPI server

# Deploys built wheels to PyPI
pip install twine
ls -l dist/*  # List files to upload
twine upload --repository-url https://test.pypi.org/legacy/ --username albarji dist/*
