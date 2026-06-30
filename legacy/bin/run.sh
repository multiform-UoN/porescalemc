#!/bin/bash

. $WM_PROJECT_DIR/etc/bashrc
. $WM_PROJECT_DIR/bin/tools/RunFunctions

##--- ANACONDA PYTHON
#export PATH=~/anaconda/bin:$PATH
pyt=python3 # PYTHON EXECUTABLES

export PYTHONPATH=$PYTHONPATH:$MLMC_PORESCALE

sh $MLMC_PORESCALE/script/setup.sh

# run the standard scripts in the code directory
$pyt -u run.py
