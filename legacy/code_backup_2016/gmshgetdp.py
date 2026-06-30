# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module to run GMSH and GETDP simulations

try:
    random
except:
    from general import *

from geomodule import rotmat_to_euler
from gmshgetdpDefault import *

try:
    exec(open("gmshgetdpDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file gmshgetdpDict.py not found or not readable. Loaded default values")

try:
    exec(open("runDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file runDict.py not found or not readable. Using module specific options")

# this works only if gmshgetdp.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("gmshgetdpDict.py").read(),globals_)


# ------------ CLASS DEFINITION FOR A SOLVER
class solver(base_solver):
    # ---------- SOLVER - CLASS INITIALIZATION
    def __init__(self, geom=dummyObj,hierarchy=""):
        self.nprocs    = nprocs
        self.force     = forcing
        self.deltap    = deltap
        base_solver.__init__(self,geom,hierarchy)

    # ----------- SET SOLVER INPUT
    # this function is called to initialize the input (with update =0)
    # and to keep the same sample and solve it at different levels (update=1)
    def set_input(self,geom):
        if (geom.level==self.level):  # initialization (copy global variables)
            self.update=False
            # initial grid is not changing
            self.gridres    = gridres
            self.scalegrid  = scalegrid     # scale the domain
        else:
            self.totcells = self.totcells*(refratio**dimension)**(geom.level-self.level)
            self.update=True
            self.level=geom.level
        # ------ parameters always set, also in the update
        if (regularmesh and 'g' in self.hierarchy):
            self.gridres     = gridres * refratio**self.level
            self.refinement  = refinement
        elif ('g' in self.hierarchy):
            self.refinement  = refinement + self.level
        else:
            self.refinement  = refinement
        if ('t' in self.hierarchy):
            self.tol   =  get_tol(self.level) # adapt the tolerance to the level
        else:
            self.tol  =  tolerance
        self.totcells   = self.gridres**dimension*geom.xlen*geom.ylen*geom.zlen
        self.name       = geom.name
        #self.sample(update)
    # --------------------------

    # ---------- RANDOM SAMPLING
    def sample(self):
        if (not self.update):
            if (random_bc=="uniform"): #TODO other distributions
                self.deltap+=random.random()*random_force_cov*self.deltap
            if (random_force=="uniform"):
                self.force+=random.random()*random_bc_cov*self.force
        return
    #---------------------------------


    # -------- SETUP and RUN SOLVER
    def setup(self,geom,sample_geo=True):
       # TODO gestire self.nprocs automaticamente
        if geom.level!=self.level:
            self.set_input(geom)
        self.name=self.basename+"_l"+str(self.level)
        geom.name=self.name
        #print("------------- setup "+self.name)
        shutil.rmtree(workdir+self.name, ignore_errors=True)
        shutil.copytree("gmsh_template",workdir+self.name)
        if sample_geo:
            geom.sample()
        self.sample()
        self.nprocs=min(self.nprocs,self.max_procs())
        log=self.snappy(geom)
        return log

    def solve(self):
        basename=self.name
        log=runPos(self.name,pdeproblem,self.tol,self.nprocs)
        res=self.qoi()
        return res, log

    def close(self,clear):
        if (clear):
            for l in range(self.level+1):
                name0=self.name[:-1]+str(l)
                shutil.rmtree(workdir+name0, ignore_errors=True)


    # -------- GMSH/GETDP SPECIFIC ROUTINES
    # CREATE GEOMETRY FILES
    def snappy(self,g):
        #print("Create gmsh input files")
        name=g.name
        tolmerge=pow(refratio,-8-g.level) # geometrical tolerance

        # ------- SCALED QUANTITIES
        mu=g.mu*scalegrid
        self.xlen=g.xlen*scalegrid
        self.ylen=g.ylen*scalegrid
        self.zlen=g.zlen*scalegrid

        # ------- WRITE FILES
        filename=workdir+name+"/grains.geo"
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        f=open(filename,'w')

        strout1=""
        strout2=""
        strout3=""
        strout4=""
        strout5=""
        strout6=""
        strout7=""
        strout8=""
        strout9=""
        strout10=""
        rr=g.get_radii()
        f.write("ngrains = %05d"%len(g.g)+";\n")
        for i in range(len(g.g)):
            rx, ax, tx = rotmat_to_euler(g.get_transformation_matrix(i))
            if sum(array(ax)**2)==0:
                ax=[1,0,0]
                tx=0
            strout1+=str(scalegrid*g.g[i,0])+","
            strout2+=str(scalegrid*g.g[i,1])+","
            strout3+=str(scalegrid*g.g[i,2])+","
            strout4+=str(rx[0])+","
            strout5+=str(rx[1])+","
            strout6+=str(rx[2])+","
            strout7+=str(ax[0])+","
            strout8+=str(ax[1])+","
            strout9+=str(ax[2])+","
            strout10+=str(tx)+","
        strout = "xx = {"+strout1[:-1]+"};\nyy = {"+\
                strout2[:-1]+"};\nzz = {"+strout3[:-1]+"};\nrr1 = {"+\
                strout4[:-1]+"};\nrr2 = {"+strout5[:-1]+"};\nrr3 = {"+\
                strout6[:-1]+"};\nax1 = {"+strout7[:-1]+"};\nax2 = {"+\
                strout8[:-1]+"};\nax3 = {"+strout9[:-1]+"};\ntx = {"+\
                strout10[:-1]+"};\n"
        f.write(strout)
        f.close()

        filename=workdir+name+"/mesh_size.geo"
        f=open(filename,'w')
        refmax = 1./self.gridres
        refmin = refmax * (refratio**-self.refinement)
        f.write("ref1 = "+str(refmax)+";\n")
        f.write("ref3 = "+str(refmin)+";\n")
        f.close()
        filename=workdir+name+"/box.geo"
        f=open(filename,'w')
        f.write("x1 ="+str(-0.5*self.xlen)+";\n")
        f.write("y1 ="+str(-0.5*self.ylen)+";\n")
        f.write("z1 ="+str(-0.5*self.zlen)+";\n")
        f.write("x2 ="+str(0.5*self.xlen)+";\n")
        f.write("y2 ="+str(0.5*self.ylen)+";\n")
        f.write("z2 ="+str(0.5*self.zlen)+";\n")
        f.write("scalegrid = "+str(self.scalegrid)+";\n") #TODO not used now
        f.close()
        filename=workdir+name+"/DirBCs"
        f=open(filename,'w')
        f.write("Function {\n")
        f.write("inletp ="+str(self.deltap*self.xlen)+";\n")
        f.write("outletp ="+str(0)+";\n")
        f.write("force ="+str(self.force)+";\n")
        f.write("flux ="+str(flux)+";\n")
        f.write("permeability ="+str(perm1)+";\n")
        f.write("permeability2 ="+str(perm2)+";\n")
        f.write("}\n")
        f.close()
        # ------------------------------
        ## stats of output grain size distribution (not scaled)
        #log="Mean="+str(mean(g.g[:,3]))+"microns\n"
        #log=log+"Dev.Std.="+str(samplestd(g.g[:,3]))+"\n"
        #log=log+"Number of grains:"+str(g.n_grains)+"\n"
        #log+log+"Expected porosity:"+str(g.porosity_out)+"\n\n\n"

        return "".encode()
    #-----------------------

    #----------------- COMPUTE THE QUANTITY OF INTEREST
    # TO DO: now it's only flux (and related quantity), extend to porosity and others
    def qoi(self):
        name=self.name
        def read_dat_gmsh(ff):
            filename=workdir+name+"/"+ff+".dat"
            try:
                if os.stat(filename).st_size > 0:
                    f=open(filename,'r')
                    lines=f.readline()
                    try:
                        return abs(float(lines[2:]))
                    except ValueError:
                        return nan
                else:
                    return nan
            except OSError:
                return nan

        intout   = read_dat_gmsh("IntOut")
        intvol   = read_dat_gmsh("Int")
        flux     = read_dat_gmsh("IntFlux")
        intvol2   = read_dat_gmsh("Int2")
        flux2     = read_dat_gmsh("IntFlux2")
        porosity = read_dat_gmsh("Por")
        permeability=nan
        if (not isnan(flux).any()):
            permeability=flux/(self.deltap+small)
        # TODO other QoI
        val=[porosity, permeability, intout, flux, intvol, flux2, intvol2]
        return array(val)

    #---------------------



    def max_procs(self):
        return max(1,int(self.totcells/4000)) # 4000 is related to the setup of GAMG solver

    # -------------------------
# ------------------------- END SOLVER CLASS

def runPos(name,pde,tol,parallel=1):
    #print("---------running "+name)
    #print("gmsh")
    log=subprocess.check_output(["gmsh", "-3", workdir+name+"/"+geom_type+".geo"])
    #print("getdp")
    if (parallel>1):
        log+=subprocess.check_output(["mpirun", "-np", str(parallel), "getdp", workdir+name+"/"+pde+".pro", "-msh", workdir+name+"/"+geom_type+".msh", "-solve", "Pde" , "-pos", "Pde", "-log", "-ksp_type", "gmres", "-pc_type ilu", "-ksp_rtol", str(tol)])
    else:
        log+=subprocess.check_output(["getdp", workdir+name+"/"+pde+".pro", "-msh", workdir+name+"/"+geom_type+".msh", "-solve", "Pde", "-pos", "Pde", "-log", "-ksp_type", "gmres", "-pc_type ilu", "-ksp_rtol", str(tol)])
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
