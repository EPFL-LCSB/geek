#!/bin/bash

source $HOME/miniconda/bin/activate
conda config --add channels conda-forge

# optional: create environment for readdy, switch to that environment
conda create -n readdy python=3
source activate readdy

# install readdy
conda install -c readdy readdy

# install dependencies on the virtual environment
conda install pandas
conda install scipy