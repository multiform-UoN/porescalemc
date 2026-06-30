# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module for running single simulations without MLMC

try:
    random
except:
    from general import *

from singleDefault import *

try:
    exec(open("singleDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file singleDict.py not found or not readable. Loaded default values")

try:
    exec(open("runDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file runDict.py not found or not readable. Using module specific options")

# this works only if single.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("singleDict.py").read(),globals_)


class estimator:
    def __init__(self, name,sol,geo):
        self.name=name
        #self.pdeproblem=pdeproblem
        self.geo=geo
        self.sol=sol

    # ---------------------------------
    def solve(self):
        t=time.time()
        name=self.name
        g=geo_interface(0,name)
        s=sol_interface(g)
        log=s.setup(g)
        [res,log1]=s.solve()
        res=array(flatten(final_qoi(res,g,s)))
        log+=log1
        s.close(cleardat)
        self.q=res
        self.w=time.time()-t
        self.g=g
        self.s=s
        return res,time.time()-t

    # ---------------------------------
    def generate(self):
        name=self.name
        g=geo_interface(0,name)
        s=sol_interface(g)
        log=s.setup(g)
        g.write()
        return log





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
