#!/bin/bash

pyt=python3 # PYTHON EXECUTABLES

cmd="$pyt -u $MLMC_PORESCALE/code_dev/mlmc_hybrid.py $1"
cond=$(ps -Af | grep "$cmd" | wc -l)

# run the rpyc server
if [ $cond -le $2  ]; then
	$cmd
else
    echo "$2 rpyc servers already running"
fi

