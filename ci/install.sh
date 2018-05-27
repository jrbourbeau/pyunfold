#!/bin/bash

# Install conda
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    MINICONDA_FILENAME="Miniconda3-latest-MacOSX-x86_64.sh"
else
    MINICONDA_FILENAME="Miniconda3-latest-Linux-x86_64.sh"
fi

wget https://repo.continuum.io/miniconda/$MINICONDA_FILENAME -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda config --set always_yes yes --set changeps1 no
# Useful for debugging any issues with conda
conda info -a

# Create conda environment
conda create --quiet --name test-environment python=$PYTHON
conda env update --quiet --name=test-environment --file=environment.yml
source activate test-environment

# Install pyunfold + dev dependencies
pip install .[dev]
echo conda list
conda list
echo pip list
pip list
