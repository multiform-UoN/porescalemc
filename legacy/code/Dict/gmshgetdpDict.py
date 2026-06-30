# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module for running GMSH and GETDP (Elliptic solvers)

#  --- mesh settings
dimension   = 3            # 2 or 3 dimensional simulations (overwritten by runDict)
regularmesh = True         # False to adapt on grains, True refine everywhere
gridres     = 10           # initial grid size x (y and z not needed for gmsh)
scalegrid   = 1.0
refratio    = 2            # refinement ratio between level (overwritten by runDict)
geom_type   = "perforated" # "disc_coeff" or "perforated"
#  --- simulation settings
nprocs   = 1             # max number of processes for openfoam solver, or keywork "max"
deltap    = 1            # Dirichlet BCs
forcing   = 0            # source terms
flux      = 1            # Neumann BCs
perm1 = 1         # coefficient in the elliptic operator
perm2 = 10         # coefficient in the elliptic operator


# ---  other simulation settings for single studies
tolerance = 1e-4
refinement= 1

# ----- random PDE
## uniform (0,1) or .....#TODO
random_force = "uniform"
random_bc    = "uniform"
random_force_cov = 0.1
random_bc_cov    = 0.1

# --- OTHER PARAMETERS overloaded by the study module (see examples in mlmcDict.py or singleDict.py)
pdeproblem = "NavierStokes"

workdir  = "./"  # overwritten by runDict

#  --- SOLVER TOLERANCE
def get_tol(level):
    return tolerance*pow(3,-level)

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
