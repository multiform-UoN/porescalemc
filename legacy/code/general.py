# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module for global variables and modules

# ----- EXTERNAL MODULES
import glob
import string
import subprocess
import os
env = os.environ

import shutil
import math
from numpy import *
import random
import string
import time
import pickle
import copy
import queue
# -----------------

# ---- OPTIMIZATION
try:
    from numba import jit
except:
    print("Module numba not available")
    def jit(func):
        return func
# -----------------

# ---- DEBUG AND LOGGING
#import pdb
import logging
LOGGING_LEVELS = {'critical': logging.CRITICAL,
                  'error': logging.ERROR,
                  'warning': logging.WARNING,
                  'info': logging.INFO,
                  'debug': logging.DEBUG}
# -----------------

# ---- PLOTTING MODULES
#import prettyplotlib as ppl
#from prettyplotlib import brewer2mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#font = {'family' : 'normal','weight' : 'normal','size'   : 18}
font = {'size'   : 18}
plt.rc('font', **font)
plt.rc('figure', dpi=120)#, figsize=(8,8))
plt.rc('lines', linestyle='dashed', linewidth=3, marker='o', markersize=10, mew=3)
#plt.Artist().set_clip_on(False)
ffloat=plt.ScalarFormatter()
ffloat.set_scientific(True)
ffloat.set_powerlimits((-2,2))
#plt.ticklabel_format(axis='both', style='sci', scilimits=(-2,2), useOffset=False)
#plt.get_major_formatter().set_powerlimits((0, 1))
# -----------------


# ---- MPLTOOLS
try:
    from mpltools import special
except:
    print("Module mpltools not available")
    special=plt
    special.errorfill=plt.errorbar
# -----------------

# ---- SEABORN
try:
    import seaborn as sns
    sns.set(font_scale=2)
    sns.set_style("white")
    def ncolors(n):
        return sns.husl_palette(n)
    colors = [(0, 0, 0)]+sns.color_palette()
except:
    print("Module seaborn not available")
    colors = ['k', 'c', 'm', 'y', 'k', 'r', 'g', 'b']
    def ncolors(n):
        return cm.rainbow(linspace(0, 1, n))
# ---------------------

class dummyStruct():
    def __init__(self,name=""):
        self.name = name
        self.level = 0
        self.xlen = 1
        self.ylen = 1
        self.zlen = 1
        self.mu   = 0

dummyObj=dummyStruct()

# ---- GLOBAL VARIABLES
small = 1e-16
third = 1.0/3.0
half  = 0.5
one   = 1.0
two   = 2.0
three = 3.0
pdf   = 'pdf'
png   = 'png'
eps   = 'eps'
unit_vect = ones(3)/sqrt(3.)
pi=math.pi
# -----------------

# ------- ADDITIONAL MATHEMATICAL FUNCTIONS
@jit
def samplevar(x,ax=None):
    return var(x,ax)*len(x)/(len(x)-1)

@jit
def samplestd(x,ax=None):
    return sqrt(samplevar(x,ax))

@jit
def modsign(x,y): # mod operation with sign
    #return (x%y)*math.copysign(1,x)
    return ((x+y)%(2*y))-y
# -----------------

# --------- OTHER FUNCTIONS OR VARIABLES
flat = lambda *n: (e for a in n
    for e in (flat(*a) if isinstance(a, (tuple, list, ndarray)) else (a,)))
def flatten(l):
    return list(flat(l))
# ----------------

# -------------- GENERIC SYSTEM FUNCTIONS
def get_cpu():
    return multiprocessing.cpu_count()

# get used RAM memory
def memory():
    try:
        mem=open('/proc/meminfo', 'r')
        tmp = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) == 'MemTotal:':
                total = int(sline[1])
            elif str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                tmp += int(sline[1])
        return 1-tmp/total
    except:
        return 0.
# -----------------

# ------ BASE CLASSES
class base_solver:
    # ---------- SOLVER - CLASS INITIALIZATION
    def __init__(self, geom=dummyObj,hierarchy=""):
        self.level     = geom.level
        self.name      = geom.name
        self.basename  = geom.name
        self.hierarchy = hierarchy
        self.set_input(geom)

# TODO operator overloading and other basic functions for base classes

# TODO other base classes
