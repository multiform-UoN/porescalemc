# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# random.py dictionary for random packing inputs

#  --- grain packing settings
psd        = "uniform"    # grain size distribution (uniform, lognormal, constant)
coeffvar   = 0.4          # coefficient of variation of the grain size distribution
detached   = False        # overlapping or detached grains (required for GMSH)
jodreytory = True         # jodrey-tory algorithm
detachedbc = True         # overlapping or detached to the boundaries (required for GMSH)
ellips     = False        # using ellipsoids of spheres (only openfoam now)
minpor     = 0.5          # minimum porosity allowed
mu         = 0.1          # mean grain size (overwritten if included in the MLMC hierarchy)
ngrains    = 10           # total grain size (overwritten if included in the MLMC hierarchy)
max_tries  = 10000        # number of tentatives to place randomly non-overlapping spheres
eps_jd     = 0.2          # jodrey-tory displacement (in overlapping distance units)
min_dist   = 1            # jodrey-tory distance threshold % (1=touching, <1 allow some overlap, >1 force distancing)
max_dist   = -1           # jodrey-tory distance threshold % (1=touching, <1 allow some overlap, >1 force distancing, <0 disabled)
cluster    = 3          # number of grains to consider for maximum distance
nmoves_jd  = 1          # number of grain to displace each JD iteration
regions    = 1          # number of separate regions (max 2 up to now, works only if detached=false)
region_ratio = 2        # ratio between regions
periodic = False        # works only for Jodrey-Tory so far #TODO

#  --- domain settings
xlen      = 1.0          # length of the domain
ylen      = 1.0          # width of the domain
zlen      = 1.0          # depth of the domain
voidspace = 1+coeffvar   # void space at the beginning and end of domain in x direction (mu units), set it to very small number if detachedbc is enable
dimension   = 3         # 2 or 3 dimensional simulations (overwritten by runDict)

#  --- other settings  to create a smooth continouos field
compute_field       = False
smoothing           = "gaussian"  # kernel for smoothing porosity field (gaussian, cell, none)
smoothing_length    = 2           # kernel size length if constant_smoothing=True otherwise number of cells
constant_smoothing  = True        # constant or mesh-dependent smoothing
resolution          = 20          # points per unit length (overwritten by the solver)
upscaling_type      = "kozeny"    # permeability form law: kozeny, poiseuille, hill, ergun, wenyu, mixed, custom
upscaling_constant  = 5           # a-priori knowledge of the constant scaling factor
upscaling_limiter   = "min"       # correct permeability to be consistent with a piecewise porosity by limiting porosity range
invkglob            = 1           # a-priori knowledge of the inverse permeability (or other effective parameter under study)

# gmsh, openfoam, blender or none
geom_out="gmsh"
workdir  = "./"  # overwritten by runDict


#  --- MEAN GRAIN SIZE
def get_mean(level):
    return mu*pow(2.0,-level*0.666)

#  --- NUMBER OF GRAINS
def get_n(level):
    return ngrains*(level+1)#pow(2,level)

#  --- DOMAIN SIZE SCALING
def get_len(level):
    return pow(2,level)
# ---------------------------------

#  --- GEOMETRY ITERATIONS
def get_max_tries(level):
    return max_tries*pow(2,level)
# ---------------------------------

# --- Nathan's project on fractures
EPS = 0.15*xlen
Prob = 0.5
nfractures = 0
randomize_frac=0.3

#Matrix=tmp[0]
#index=tmp[1]
#PointX = tmp[2]
#PointY = tmp[3]
#PointZ = tmp[4]
# --------------------------------

# --- Approximate upscaling function
# gives back the inverse permeability
def custom_porosity_transform(eps,mu,re):
    return 1e8*(1-eps)+1e3


# ------------------------------------------------- END FILE


#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
