#!/bin/bash
arg=$1
find . -iname "${arg}_*" -exec rm -rf {} \;
