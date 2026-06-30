#!/bin/bash
cat $1 | grep areaAverage | cut -d' ' -f9 | tr -d ',(' > fluxinlet
#head -n 100 $1 | grep 'with volume' | cut -d' ' -f11 | tr -d ',(' > volume
cat $1 | grep 'getMeshStats volume' | cut -d' ' -f 3 > volume
cat $1 | grep 'getMeshStats area pores' | cut -d' ' -f 4 > surface
