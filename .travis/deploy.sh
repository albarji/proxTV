#!/usr/bin/env bash
set -v

# Deploys built wheels to PyPI
pip install twine
twine upload --repository-url https://test.pypi.org/legacy/ --username albarji dist/*
