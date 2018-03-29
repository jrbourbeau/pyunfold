#!/usr/bin/env bash

set -e

echo "Building documentation..."
pip install -r requirements/docs.txt
cd docs
make clean
make html
cd ../
echo "Successfully built the documentation!"
