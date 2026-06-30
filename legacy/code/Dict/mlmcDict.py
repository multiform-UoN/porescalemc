# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# mlmc.py dictionary


#  --- simulation settings
cleardat       = False        # for big simulations, this clear the final CFD data to save space
# nprocs         = 1            # number of processes for openfoam solver, or keywork "max"
# REMOVED ---- TODO move all in the solver dictionary
nthreads       = 1            # number of threads for parallel sampling, or keyword "adaptive"
writelog       = True         # write log files #TODO
distributed    = False        # TODO
nthreads_scale = 8            # scale down nthreads for each level for memory limitations #TODO
intermediate_output   = False # for long run write qoi and work objects during the run
intermediate_output_n = 100   # number of samples after which results are stored
intermediate_output_t = 3600  # time (in seconds) after which results are stored
ramlimit       = 0.8          # max percentage of RAM usage (only works in linux)
workdir  = "./"               # overwritten by runDict

#  --- MLMC settings
mlmc_algorithm = 'none'    # giles (adaptive), none (non-adaptive), mesh_convergence (single stochastic realization on multiple levels)
mlmc_estimator = 'pair'    # pair (default), triplet (experimental) or test (for debugging)
mlmc_mratio    = 4.0       # ratio between number of samples between two level
tot_level      = 3         # total number of levels (initial guess)
max_level      = 5         # maximum level (adaptive algorithm cannot exceed this level)
m0             = 2*mlmc_mratio**(tot_level)  # number of samples at level 0
confidence     = 0.99      # confidence
final_tol      = 1e-2      # overall tolerance TOL (can be a list as long as the final number of QoIs)
error_split    = 0.5       # constant theta, weight of bias wrt stat err
reuse_sample   = True      # reuse the samples from previous iterations or recompute each time
active_var     = -1        # variable to use for the error estimation (among the ones defined in final_qoi in runDict), -1 for all
scale_qoi      = "none"    # scale or shift the qoi at each level to match the expected value ("sum" or "prod")
refratio       = 2         # refinement ratio between level (everything not set separately in ad-hoc functions) (overwrite solver settings, overwritten by runDict setting)
update_ml_relaxation = 0.9 # relaxation rate (between 0 and 1) to update number of samples. 1=use formula 0=never change
min_samples    = 5         # minimum number of samples per level for stability issues

# ---  convergence properties
# extrapolation
extrapolation  = "richardson" # extrapolate to reduce bias richardson or aitken #TODO
# assumed convergence rate (#TODO)
beta        = 1.2
alpha       = 2.5
gamma       = 3.0

# ------------- type of MLMC hierarchy
# g=grid
# d=domain
# f=fourier space discretization (only for darcyOF)
# s=stone size
# n=stone numbers
# r=stone roughness
# p=placement tentatives
hierarchy   = "gt"         # classical MLMC
#hierarchy  = "dgsnrp"     # multi-scale multi-level MC
mimc        = False        # Multi Index Monte Carlo - #TODO


#------- TWO-PASS AVERAGING
# Modify samples with available average from first-pass
def qoi_twopass(avg):
    # It must return a function to modify the final_qoi vector
    def myfun(q):
        #-----------
        ## this is an example to compute variance
        #import numpy as np
        #q2=(q-avg)**2
        #return np.concatenate([q,q2])
        #-----------
        # this is an example that does nothing
        return q
    #-----------
    return myfun
# ---------------------------------

#------- POSTPROCESSING QUANTITIES
# After the averaging process
def final_postprocess(avg,err,bias):
    # It must return an array of result with consistent error estimation
    # result contains the average of the final_qoi vector and the absolute error associated
    #-----------
    ## this is an example that simply remove the second half of the vector that is added automatically by the code to compute covariances between levels (default)
    import numpy as np
    nvar=len(avg)/2
    result=avg[:nvar]
    error=np.abs(bias[:nvar])+np.abs(err[:nvar])
    #-----------
    ## add this part to compute variance or STD when final_qoi contains [qoi1, qoi2, qoi1**2, qoi2**2]
    #l=len(result)/2
    #def var(r):
        #return np.sqrt(np.maximum(1e-10,r[l:]-r[0:l]**2))
    #result=np.concatenate([result[0:l],var(result)])
    #error=np.concatenate([error[0:l], [0.5*(error[l:]+2.*np.abs(result[0:l])*error[0:l])*var(result)**-.5])
    #-----------
    return result, error
# ---------------------------------

#------- SAMPLE ACCEPTANCE/REJECTANCE criteria (should return a boolean)
def accept(result,geometry_object=None,solver_object=None):
    # result is a list of output (see function above)
    # it should return either True or False
    return True
# ---------------------------------


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
