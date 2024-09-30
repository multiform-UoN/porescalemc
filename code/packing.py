# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module to run packing statistical analysis (no PDE solver)

try:
    random
except:
    from general import *

from packingDefault import *

try:
    exec(open("packingDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file packingDict.py not found or not readable. Loaded default values")


try:
    exec(open("runDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file runDict.py not found or not readable. Using module specific options")

# this works only if packing.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("packingDict.py").read(),globals_)


# ------------ CLASS DEFINITION FOR A SOLVER
class solver(base_solver):
    # ---------- SOLVER - CLASS INITIALIZATION
    def __init__(self, geom=dummyObj,hierarchy=""):
        self.nprocs    = nprocs
        base_solver.__init__(self,geom,hierarchy)

    # ----------- SET SOLVER INPUT
    # this function is called to initialize the input (with update =0)
    # and to keep the same sample and solve it at different levels (update=1)
    def set_input(self,update=None):
        # ------ parameters always set, also in the update
        if (update is not None):
            self.level=update
    # --------------------------

    # ---------- RANDOM SAMPLING
    def sample(self):
        # now empty, it can be used to introduce further randomness in the solver parameters or physical parameters
        return
    #---------------------------------

    # -------- SETUP and RUN SOLVER
    def setup(self,geom):
        self.level=geom.level
        self.name=geom.name+"_l"+str(self.level)
        geom.name=self.name
        geom.sample()
        self.sample()
        self.p=geom.porosity_out
        self.n=geom.n_grains_out
        #self.n_intersections=geom.n_intersections
        #self.ntries_out=geom.ntries_out
        if ("gmsh" in geom_out):
            import gmshgetdp as gg
            gg.pdeproblem=""
            outg=gg.solver(geom)
            geom.name="."
            outg.snappy(geom)
        if ("openfoam" in geom_out):
            import openfoam as gg
            gg.pdeproblem="laplacian"
            outg=gg.solver(geom)
            geom.name="."
            outg.setup(geom,sample_geo=False)
            outg.snappy(geom)
        if ("openscad" in geom_out):
            geom.name="."
            geom.write_openscad()
        if ("stl" in geom_out):
            geom.name="."
            geom.write_stl()
        return ''.encode('ascii')

    def solve(self):
        log=''.encode('ascii')
        res=self.qoi()
        return res, log

    def close(self,cleardat):
        return


    #----------------- COMPUTE THE QUANTITY OF INTEREST
    def qoi(self):
        val=[self.p, self.n]#, self.n_intersections, self.ntries_out]
        return array(val)

    #---------------------

    def max_procs(self):
        return 1

# ------------------------- END SOLVER CLASS

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
