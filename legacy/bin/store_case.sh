#!/bin/bash
arg=$1
mkdir $1
find . -maxdepth 1 -iname "${arg}_*" -exec mv {} $arg \;
mv -f $arg*Q $arg*W nohup.out $1
cp mlmcDict.py $1
