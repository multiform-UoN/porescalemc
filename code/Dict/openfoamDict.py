# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# dictionary for OpenFOAM solver


#  --- mesh settings
dimension   = 3         # 2 or 3 dimensional simulations (overwritten by geom module)
regularmesh = True      # False to adapt on grains, True refine everywhere
gridres     = 10        # initial grid resolution (per unit length)
scalegrid   = 1.0
voxelized   = True      # skip snapping process 
refratio    = 2         # refinement ratio between level (overwritten by runDict)
inside      = True      # mesh also inside the grains

#  --- simulation settings
nprocs    = 1               # max number of processes for openfoam solver, or keywork "max"
pressure_source = False     # pressure imposed as a source term instead of a BCs
# BCs for mean variable (pressure in NS is imposed automatically)
# minX,maxX,lateral, pores
# symmetry,symmetryPlane,fixedValue,zeroGradient, fixedGradient, cyclic, cyclicAMI, etc.
# bcs      = ["pressureInletVelocity","pressureInletOutletVelocity", "symmetryPlane", "fixedValue"] STANDARD NS
bcs      = ["pressureInletVelocity","pressureInletOutletVelocity", "symmetryPlane", "fixedValue"]
bcsvalue = [0,0,0,0]
# same for transport (only for multiphase flows)
bcs3      = ["fixedValue","zeroGradient", "symmetryPlane", "zeroGradient"] # for scalar transport
# bcs3      = ["uniformFixedValue","zeroGradient", "symmetryPlane", "constantAlphaContactAngle"] # for interFoam
bcs3value = [1,0,0,0]

#  ---- thermophysical conditions
deltap            = 1        # deltap per unit meter
p0                = 0         # static pressure
capillarypressure = 0         # capillary pressure (unscaled)
contactangle      = 90        # for multiphase simulations
temperature       = 330       # constant temperature
viscosity         = 1         # fluid dynamic viscosity
perm1             = 1e-3      # permeabilities for discountinous coefficient simulations - inside
perm2             = 1         # outside
# NB: when using darcyOF (heterogeneous interpolated coefficient computed with an upscaling formula defined in randomgeoDict, it should be perm1=0 and perm2=1)

# --- time stepping and output (per unit length) (add *10**-4 for unsteady)
deltat    = 1        # initial time step
timeinj   = 200      # injection time (time needed to stabilize the initial conditions)
timetot   = 20000    # total simulation time # THIS CONTROLS ALSO MAX ITER OF STEADY SIMULATIONS
timeoutf  = 1000     # output time step (full data)
timeout   = 10       # output time step (averaged data) # THIS CONTROLS ALSO OUTPUT OF STEADY SIMULATIONS

# ---  other simulation settings for continuous field simulation (other settings are in randomgeoDict)
min_cells_per_grain = 0 # decide when to activate the continuous formulation, 0 to disable upscaling
field_resolution    = 1 # (per unit cell), overwrite resolution in randomgeoDict

# --- OTHER PARAMETERS, they may be overloaded by the study module (see examples in mlmcDict.py or singleDict.py)
aitken_extrapolation = False  # use aitken acceleration for steady state problem (to use in conjunction with a much smaller tolerance)
# TODO implement adaptive stopping based on aitken. when convergence is linear, aitken is not working

pdeproblem = "NavierStokes"
tolerance  = 1e-5
refinement = 0
min_ncells = 2  # minimum refinement layer (cells between 2 layers)
stretch_factor = 1.5  # refinement region stretch factor (each level, the refinement will grow of this factor)

workdir  = "./"  # overwritten by runDict

#  --- SOLVER TOLERANCE
def get_tol(level):
    return tolerance*pow(3,-level)

#  --- refinement speed, number of each between layers
def get_ncells(level):
     #return int(min_ncells+stretch_factor**level)
     return max(min_ncells,level)

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
