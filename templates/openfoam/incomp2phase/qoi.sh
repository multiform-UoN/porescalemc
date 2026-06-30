#!/bin/bash
cat $1 | grep MassFlows | cut -d' ' -f7 | tr -d ',(' > massfluxinlet
cat $1 | grep MassFlows | cut -d' ' -f11 | tr -d ',(' > massfluxoutlet
cat $1 | grep 'Averages of alpha.water' | cut -d' ' -f9 | tr -d ',(' > wateroutlet
cat $1 | grep 'Averages of alpha.air' | cut -d' ' -f9 | tr -d ',(' > airoutlet
cat $1 | grep 'Integral of U' | cut -d' ' -f9 | tr -d ',(' > fluxoutlet
cat $1 | grep 'Integral of U' | cut -d' ' -f15 | tr -d ',(' > fluxinlet
cat $1 | grep 'Integral of alpha.water' | cut -d' ' -f6 | tr -d ',(' > watersaturation
cat $1 | grep 'Integral of alpha.air' | cut -d' ' -f6 | tr -d ',(' > airsaturation
cat $1 | grep 'Time =' | cut -d' ' -f3 | tr -d ',(' | paste - - -d' ' > t
#head -n 100 $1 | grep 'with volume' | cut -d' ' -f11 | tr -d ',(' > volume
cat $1 | grep 'getMeshStats volume' | cut -d' ' -f 3 > volume
cat $1 | grep 'getMeshStats area pores' | cut -d' ' -f 4 > surface
