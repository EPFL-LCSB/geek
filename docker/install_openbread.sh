#!/bin/bash

source /home/user/openfpm_vars

git clone -b revision https://github.com/EPFL-LCSB/openbread
pip3 install -e openbread
