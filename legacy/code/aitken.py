# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# January 2015
#--------------------------------------------------------------------------------------------
# module for aitken extrapolation

from math import *
from numpy import *


def aitken(x):
#    return [(x[i]*x[i-2]-x[i-1]**2)/(x[i]-2*x[i-1]+x[i-2]) for i in range(2,len(x))]
    return array([x[i-2] - (x[i-2]-x[i-1])**2/(x[i]-2*x[i-1]+x[i-2]) for i in range(2,len(x)) if abs(x[i]-2*x[i-1]+x[i-2])>1e-12])


def aitken_process(f):
    lines=[]
    with open(f, 'rb') as csvfile:
        for row in csvfile:
            #print(row.split())
            lines.append([float(i) for i in row.split()])
    x=transpose(lines)
    if len(x)<2:
        print("Warning: Not enough iterations")
        return x[-1]
    
    y=aitken(x[0])

    with open(f+'_aitken', 'w') as f:
        for row in y:
            f.write(str(row)+"\n")
    
    y=y[int(len(y)*.95):]
    if std(y)>abs(y[-1])*.1:
        print("Warning: aitken not converging: returning original value")
        return x[-1]
    return mean(y)


