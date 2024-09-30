# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module for running OpenFOAM solver

try:
    random
except:
    from general import *

from openfoamDefault import *

try:
    exec(open("openfoamDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file openfoamDict.py not found or not readable. Loaded default values")

try:
    exec(open("runDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file runDict.py not found or not readable. Using module specific options")

# this works only if openfoamDict.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("openfoamDict.py").read(),globals_)

# ------------ CLASS DEFINITION FOR A SOLVER
class solver(base_solver):
    # ---------- SOLVER - CLASS INITIALIZATION
    def __init__(self, geom=dummyObj,hierarchy=""):
        self.nprocs    = nprocs
        #self.force     = forcing
        self.deltap    = deltap
        base_solver.__init__(self,geom,hierarchy)

    # ----------- SET SOLVER INPUT
    # this function is called to initialize the input (with update =0)
    # and to keep the same sample and solve it at different levels (update=1)
    def set_input(self,geom):
        # hierarchy of meshes (classical MLMC)
        if (geom.level == self.level):  # initialization (copy global variables)
            self.update=False
            # initial grid is not changing
            self.gridres    = gridres
            self.scalegrid  = scalegrid     # scale the domain
        else:
            self.update=True
            self.level=geom.level
        # ------ parameters always set, also in the update
        if (regularmesh and 'g' in self.hierarchy):
            self.gridres      = gridres * refratio**self.level
            self.refinement = refinement
        elif ('g' in self.hierarchy):
            self.refinement=refinement + self.level
        else:
            self.refinement=refinement
        if ('t' in self.hierarchy):
            self.tol=get_tol(self.level) # adapt the tolerance to the level
        else:
            self.tol=tolerance
        if pressure_source:
            self.deltap=0
            self.sourcep=deltap
        else:
            self.deltap=deltap
            self.sourcep=0
        if (pdeproblem=="NavierStokes" and geom.mu*self.gridres*refratio**(self.refinement-refinement) <= min_cells_per_grain):
            self.pdeproblem   ="darcy"
            if ('g' in self.hierarchy):
                self.gridres      = self.gridres*refratio**(self.refinement-refinement)
        else:
            self.pdeproblem   = pdeproblem
        if 'f' in self.hierarchy:
            geom.resolution=field_resolution*self.gridres
        else:
            geom.resolution=field_resolution
        self.viscosity=viscosity
        self.totcells   = self.gridres**dimension*geom.xlen*geom.ylen*geom.zlen
        self.mu=geom.mu*scalegrid
        self.xlen=geom.xlen*scalegrid
        self.ylen=geom.ylen*scalegrid
        self.zlen=geom.zlen*scalegrid

        #self.sample(update)
    # --------------------------

    # ---------- RANDOM SAMPLING
    def sample(self):
        # now empty, it can be used to introduce further randomness in the solver parameters or physical parameters
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
        shutil.copytree(self.pdeproblem+"_template",workdir+self.name)
        print(workdir+self.name)
        if sample_geo:
            geom.sample()
        self.sample()
        self.nprocs=min(self.nprocs,self.max_procs())
        log=self.snappy(geom)
        return log

    def solve(self):
        name0=self.name[:-1]+str(self.level-1)
        if (os.path.isdir(name0)):
            log=runfoam(self.name,self.nprocs,name0)
        else:
            log=runfoam(self.name,self.nprocs)
        res=self.qoi()
        return res, log

    def close(self,clear):
        if (clear and '000000' not in self.name):
            for l in range(self.level+1):
                name0=self.name[:-1]+str(l)
                shutil.rmtree(workdir+name0, ignore_errors=True)
        #for fold in glob.glob(workdir+self.name+"/processor*"):
            #shutil.rmtree(fold, ignore_errors=True)



    # -------- OPENFOAM SPECIFIC ROUTINES
    # CREATE SNAPPYHEXMESH FILES and all the other files to run the simulation
    def snappy(self,g):
        log="Start snappy function\n"
        t=time.time();
        name=self.name
        tolmerge=pow(2,-8-g.level) # geometrical tolerance

        # ------- WRITE FILES
        filename=workdir+name+"/system/shm_objects"
        f=open(filename,'w')
        filename=workdir+name+"/system/shm_stl"
        f2=open(filename,'w')

        if False: #g.stl: #TODO
            log+=g.write_stl(workdir+name+"/constant/triSurface/")
            strout="\tgrain.stl\n\t{\n\t\ttype\ttriSurfaceMesh;\n\t\tname\tgrain;\n\t}\n"  # OPENFOAM 2.3
            f2.write(strout)
            strout="\tgrain\n\t{\n\t\tsurface grain.stl;\n\t\tscale\t(1 1 1);\n\t\ttransform\n\t\t{\n\t\t\tcoordinateSystem\n\t\t{\n\t\t\ttype\tcartesian;\n\t\t\torigin\t(0 0 0);\n\t\tcoordinateRotation\n\t\t{\n\t\ttype\taxesRotation;\n\t\te1\t(0 -1 0);\n\t\te2\t(1 0 0);\n\t\t}\n\t\t}\n\t\t}\n\t}\n"  # OPENFOAM 2.3
            f.write(strout)
        else:
            rr=g.get_radii()
            for i in range(len(g.g)):
                m=g.get_transformation_matrix(i)
                #strout="\tgrain"+"%05d"%i+"\n\t{\n\t\tsurface circle;\n\t\tscale ("+r1+" "+r2+" "+r3+");\n\t\ttransform\n\t\t{\n\t\t\ttype cartesian;\n\t\t\torigin ("+str(scalegrid*g.x[i])+" "+str(scalegrid*g.y[i])+" "+str(scalegrid*g.z[i])+");\n\t\t}\n\t}\n"  # OPENFOAM 2.2
                strout="\tgrain"+"%05d"%i+"\n\t{\n\t\tsurface\tcircle;\n\t\tscale\t("+str(scalegrid*array(rr[i])).strip("[]")+");\n\t\ttransform\n\t\t{\n\t\t\tcoordinateSystem\n\t\t{\n\t\t\ttype\tcartesian;\n\t\t\torigin\t("+str(scalegrid*g.g[i,:3]).strip('[]')+");\n\t\tcoordinateRotation\n\t\t{\n\t\ttype\taxesRotation;\n\t\te1\t("+str(m[:,0]).strip("[]")+");\n\t\te2\t("+str(m[:,1]).strip("[]")+");\n\t\t}\n\t\t}\n\t\t}\n\t}\n"  # OPENFOAM 2.3
                f.write(strout)
        f.close()
        f2.close()

        filename=workdir+name+"/system/shm_inside"
        f=open(filename,'w')
        if (inside):
            f.write("            cellZone    porous;\n")
            f.write("            faceZone    porous;\n")
            f.write("            cellZoneInside  outside;\n")
        else:
            f.write("\n")
        f.close()
        filename=workdir+name+"/system/shm_location"
        f=open(filename,'w')
        f.write("\tlocationInMesh ("+str(-0.99*g.xlen)+" "+str(0)+" "+str(0)+");")
        f.close()
        filename=workdir+name+"/system/refinement"
        f=open(filename,'w')
        f.write("refmin "+str(self.refinement)+";\n")
        f.write("refmax "+str(self.refinement)+";\n")
        f.write("ncells "+str(get_ncells(self.refinement))+";\n")
    #       f.write("d "+str(0.0001)+";\n")
    #       f.write("ref "+str(g.refinement)+";\n")
    #       f.write("dd "+str(0.01)+";\n")
    #       f.write("refref "+str(g.refinement)+";\n")
        f.write("tolmerge "+str(tolmerge)+";\n")
        if (voxelized):
            f.write("snap false;\n")
        else:
            f.write("snap true;\n")
        f.write("#inputMode merge\n")
        f.close()
        filename=workdir+name+"/constant/perm"
        f=open(filename,'w')
        f.write("perm1 "+str(perm1)+";\n")
        f.write("perm2 "+str(perm2)+";\n")
        f.write("#inputMode merge\n")
        f.close()
        filename=workdir+name+"/system/tolerance"
        f=open(filename,'w')
        f.write("tol "+str(self.tol)+";\n")
        f.write("#inputMode merge\n")
        f.close()
        filename=workdir+name+"/constant/polyMesh/domainsize"
        f=open(filename,'w')
        f.write("xout "+str(-g.xlen)+";\n")  # leave a small empty region at the beginning that will be removed by openfoam later
        f.write("x1 "+str(-0.5*g.xlen)+";\n")  # leave a small empty region at the beginning that will be removed by openfoam later
        f.write("y1 "+str(-0.5*g.ylen)+";\n")
        f.write("z1 "+str(-0.5*g.zlen)+";\n")
        f.write("x2 "+str(0.5*g.xlen)+";\n")
        f.write("y2 "+str(0.5*g.ylen)+";\n")
        f.write("z2 "+str(0.5*g.zlen)+";\n")
        f.write("xgridout "+str(int(self.gridres*g.xlen*1.5))+";\n")
        f.write("xgrid "+str(int(self.gridres*g.xlen))+";\n")
        f.write("ygrid "+str(int(self.gridres*g.ylen))+";\n")
        if dimension<3:
            f.write("zgrid "+str(1)+";\n")
        else:
            f.write("zgrid "+str(int(self.gridres*g.zlen))+";\n")
        f.write("scalegrid "+str(self.scalegrid)+";\n")
        for bc,lab in zip(bcs,["inletbc","outletbc","lateralbc","poresbc"]):
            if bc in (["symmetry","symmetryPlane","empty","cyclic","cyclicAMI"]):
                f.write(lab+"mesh "+bc+";\n")
            elif bc=="fixedValue":
                if "let" in lab:
                    f.write(lab+"mesh patch;\n")
                else:
                    f.write(lab+"mesh wall;\n")
            else:
                f.write(lab+"mesh patch;\n")
        f.write("#inputMode merge\n")
        f.close()
        filename=workdir+name+"/system/ndomains"
        f=open(filename,'w')
        if (self.nprocs>7 and dimension>2):
            self.nprocsz=self.nprocsy=int(self.nprocs**third)
            self.nprocsx=int(self.nprocs/self.nprocsy/self.nprocsz)
            while (self.nprocsx*self.nprocsy*self.nprocsz!=self.nprocs):
                self.nprocsz-=1
                self.nprocsx=int(self.nprocs/self.nprocsy/self.nprocsz)
        else:
            self.nprocsz=1
            self.nprocsy=int(self.nprocs**0.5)
            self.nprocsx=int(self.nprocs/self.nprocsy)
            while (self.nprocsx*self.nprocsy*self.nprocsz!=self.nprocs):
                self.nprocsy-=1
                self.nprocsx=int(self.nprocs/self.nprocsy)
        self.nprocs=self.nprocsx*self.nprocsy*self.nprocsz
        f.write("ndomains "+str(self.nprocs)+";\n")
        f.write("ndomainsx "+str(self.nprocsx)+";\n")
        f.write("ndomainsy "+str(self.nprocsy)+";\n")
        f.write("ndomainsz "+str(self.nprocsz)+";\n")
        f.write("#inputMode merge\n")
        f.close()
        filename=workdir+name+"/pressure"
        f=open(filename,'w')
        f.write("viscosity "+str(self.viscosity)+";\n")
        f.write("sourcep "+str(self.sourcep/self.ylen/self.zlen)+";\n")
        f.write("inletp "+str(p0+self.deltap*self.xlen+capillarypressure)+";\n")
        f.write("outletp "+str(p0)+";\n")
        f.write("temperature "+str(temperature)+";\n")
        for bc,lab in zip(bcs,["inletbc","outletbc","lateralbc","poresbc"]):
            if bc in (["symmetry","symmetryPlane","empty","cyclic","cyclicAMI"]):
                f.write(lab+"prim "+bc+";\n") # primary variable (velocity for NS, T for laplacian
                f.write(lab+"sec "+bc+";\n")  # secondary variable (pressure for NS)
            elif bc=="fixedValue":
                f.write(lab+"prim "+bc+";\n")
                if "2phase" in self.pdeproblem:
                    f.write(lab+"sec fixedFluxPressure;\n")
                else:
                    f.write(lab+"sec zeroGradient;\n")
            else:
                f.write(lab+"prim "+bc+";\n")
                f.write(lab+"sec fixedValue;\n")
        for bcval,lab in zip(bcsvalue,["inletbc","outletbc","lateralbc","poresbc"]):
            f.write(lab+"val "+str(bcval)+";\n")
        # tertiary variable (transport)
        for bc,lab in zip(bcs3,["inletbc","outletbc","lateralbc","poresbc"]):
            f.write(lab+"ter "+bc+";\n")
        for bcval,lab in zip(bcs3value,["inletbc","outletbc","lateralbc","poresbc"]):
            f.write(lab+"valter "+str(bcval)+";\n")
        f.write("#inputMode merge\n")
        f.close()
        filename=workdir+name+"/time"
        f=open(filename,'w')
        f.write("timetot "+str(timetot*self.xlen)+";\n")
        f.write("deltat "+str(deltat*self.xlen)+";\n")
        f.write("timeoutf "+str(timeoutf*self.xlen)+";\n")
        f.write("timeout "+str(timeout*self.xlen)+";\n")
        f.write("timeinj1 "+str(timeinj*self.xlen*0.999999)+";\n")
        f.write("timeinj2 "+str(timeinj*self.xlen)+";\n")
        f.write("#inputMode merge\n")
        f.close()
        filename=workdir+name+"/contactangle"
        f=open(filename,'w')
        f.write("contactangle "+str(contactangle)+";\n")
        f.write("#inputMode merge\n")
        f.close()
        # ------------------------------
        ## stats of output grain size distribution (not scaled)
        #log=log+"Mean="+str(mean(g.r))+"microns\n"
        #log=log+"Dev.Std.="+str(samplestd(g.r))+"\n"
        #log=log+"Number of grains:"+str(g.n_grains_out)+"\n"
        #log=log+"Expected porosity:"+str(g.porosity_out)+"\n"
        # ------------------------------
        # prepare the simulation for the darcy continuum solver
        log=log+"End snappy - written all setup files\n"
        if (self.pdeproblem=="darcy"):
            g.write_field(workdir+name+"/",perm2,perm1)
        return log.encode()
    #-----------------------

    #----------------- COMPUTE THE QUANTITIES OF INTEREST
    def qoi(self):
        filename=workdir+self.name+"/post_out"
        try:
            res=[[j.strip('()') for j in str.split(i)] for i in open(filename,'r').readlines()]
        except:
            print(filename)
            return [nan]
        # Quantity of interest defined in post_functions in the templates
        volume            = read_qoi(res,"domain volume SUM")
        meanvel           = read_qoi(res,"domain U SUM")/volume
        porosity          = volume/self.ylen/self.zlen/self.xlen
        watersaturation   = read_qoi(res,"domain alpha.water SUM")/volume
        eff_diffusion     = -read_qoi(res,"domain T GRADSUM")/self.xlen/self.ylen/self.zlen
        surface           = read_qoi(res,"pores area SUM")/volume
        forcep            = read_qoi(res,"pores p DIRECTEDSUM")
        forceU            = read_qoi(res,"pores U SHEAR")*viscosity
        permeability      = self.viscosity*meanvel[0]/(self.deltap+self.sourcep/self.ylen/self.zlen)
        val=[porosity,                      # 0
        surface,                            # 1
        meanvel,                            # 2
        permeability,                       # 3
        forcep,                             # 4  pressure forces
        forceU,                             # 5  viscous forces
        eff_diffusion,                      # 6
        watersaturation,                    # 7
        read_qoi(res,"domain U NORML1"),    # 8
        read_qoi(res,"domain U NORML2"),    # 9
        read_qoi(res,"domain U GRADNORML1"),# 10
        read_qoi(res,"domain U GRADNORML2"),# 11
        read_qoi(res,"domain p NORML1"),    # 12
        read_qoi(res,"domain p NORML2"),    # 13
        read_qoi(res,"domain p GRADNORML1"),# 14
        read_qoi(res,"domain p GRADNORML2"),# 15
        read_qoi(res,"domain T NORML1"),    # 16
        read_qoi(res,"domain T NORML2"),    # 17
        read_qoi(res,"domain T GRADNORML1"),# 18
        read_qoi(res,"domain T GRADNORML2"),# 29
        0]
        return val

    #---------------------

    def max_procs(self):
        np=(int(self.totcells/4000)) # 4000 is related to the setup of GAMG solver
        return max(1,np)

    # -------------------------
# ------------------------- END SOLVER CLASS

def runfoam(name,parallel=1,coarse_res="None"):
    if dimension<3:
        log=run2d(name,parallel,coarse_res)
    else:
        log=run3d(name,parallel,coarse_res)
    return log


def run2d(name,parallel=False,coarse_res="None"):
    print("2D simulations not implemented")
# REPLACE WITH ANOTHER ALLRUN FILE
#       #print("---------running "+name)
#       runlist=["foamJob", "-screen", "-case", name]
#       if parallel:
#               runlist=runlist+["-parallel"]
#       print("blockmesh 2D")
#       log=subprocess.check_output(["blockMesh", "-case", name])
#       print("snappy")
#       log=log+subprocess.check_output(["snappyHexMesh", "-case", name, "-overwrite"])
#       print("extrude+toposet")
#       shutil.copyfile(name+"/system/extrudeMeshDict1",name+"/system/extrudeMeshDict")
#       log=log+subprocess.check_output(["extrudeMesh", "-case", name])
#       log=log+subprocess.check_output(["topoSet", "-case", name])
#       print("createpatch+extrude")
#       shutil.copyfile(name+"/system/createPatchDict1",name+"/system/createPatchDict")
#       log=log+subprocess.check_output(["createPatch", "-case", name, "-overwrite"])
#       shutil.copyfile(name+"/system/extrudeMeshDict2",name+"/system/extrudeMeshDict")
#       log=log+subprocess.check_output(["extrudeMesh", "-case", name])
#       shutil.copyfile(name+"/system/createPatchDict2",name+"/system/createPatchDict")
#       log=log+subprocess.check_output(["createPatch", "-case", name, "-overwrite"])
#       log=log+subprocess.check_output(["renumberMesh", "-case", name, "-overwrite"])
#       shutil.copytree(name+"/0_2d",name+"/0")
#
# #     print("potential")
# #     log=log+subprocess.check_output(["potentialFoam", "-case", name])
# #     log=log+subprocess.check_output(["patchIntegrate", "-case", name, "U", "minX"])
#
#       if (coarse_res is not None):
#               print("interpolate")
#               log=log+subprocess.check_output(runlist+["mapFields", "../"+coarse_res, "-consistent", "-sourceTime", "latestTime"])
#       if parallel:
#               log=log+subprocess.check_output(["decomposePar", "-case", name])
#       print("simplefoam")
#       log=log+subprocess.check_output(runlist+["$(getApplication)"])
#       print("postprocess")
#       if parallel:
#               log=log+subprocess.check_output(["reconstructPar", "-latestTime", "-case", name])
# #     log=log+subprocess.check_output(["patchIntegrate", "-latestTime", "-case", name, "U", "minX"])
#
#       filename=workdir+name+"/log_all"
#       f=open(filename,'w')
#       f.write(log.decode('ascii'))
#       log=log+subprocess.check_output([name+"/qoi.sh", filename])
#       log=log+subprocess.check_output(["touch", name+"/temp.foam"])
#       f.close()
#       print("--------- end running "+name)
    return log
# -------------------------

def run3d(name,parallel=1,coarse_res=None):
    print("running "+name+" with nprocs=",parallel)
    #if (parallel>1):
    #    log=subprocess.check_output([name+'/Allrun_parallel', coarse_res])
    #else:
    log=subprocess.check_output([workdir+name+'/Allrun', coarse_res])
    filename=workdir+name+"/log_error"
    f=open(filename,'w')
    f.write(log.decode('ascii'))
    f.close()
    #print("---------end running "+name)
    return log

# --- retrieve specific data from the results vector that was stored in post_out file
def read_qoi(res,val,all_iterations=False):
    x=[[float(s) for s in i[3:]] for i in res if i[0:3]==val.split()]
    if len(x)==0:
        return nan
    if len(x[0])==1:
        x=[i[0] for i in x]
    if (all_iterations):
        return array(x)
    if (aitken_extrapolation and len(x)>10):
        import aitken
        x=aitken.aitken(x)
    if isinstance(x[-1],list):
        return array(x[-1])
    else:
        return x[-1]

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
