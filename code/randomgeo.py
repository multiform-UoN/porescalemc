# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module for generating packings with random placement

try:
    random
except:
    from general import *

from randomgeoDefault import *


try:
    exec(open("randomgeoDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file randomgeoDict.py not found or not readable. Loaded default values")


try:
    exec(open("runDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file runDict.py not found or not readable. Using module specific options")

# this works only if random.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("randomgeoDict.py").read(),globals_)

from geomodule import *

# ------------ CLASS DEFINITION FOR A GEOMETRY
class geom:
    # ---------- GRAIN GENERATION - CLASS INITIALIZATION
    def __init__(self, l=0, name=str(time.time()),hierarchy="",tot_level=0):
        self.level      = l
        self.tot_level  = max(tot_level,l)
        self.name       = name
        if (dimension<3):
            self.zlen   = 1
        self.hierarchy  = hierarchy
        self.ellipsoid  = ellips
        self.compute_field = compute_field
        self.resolution = resolution
        self.set_input()


    # ----------- SET GEOMETRY INPUT
    # this function is called to initialize the sample (with update =0)
    # and to keep the same sample and solve it at different levels (update=1)
    def set_input(self,update=None):
        if (update is None):  # initialization (copy global variables)
            self.update=False
            self.cov = coeffvar  # coefficient of variation of grains size
            self.porosity   = minpor        # minimum allowed porosity (estimated)
            # hierarchy of domain size
            if ('d' not in self.hierarchy):
                self.xlen   = xlen
                self.ylen   = ylen
                self.zlen   = zlen
            else:
                self.xlen   = xlen*get_len(self.level)
                self.ylen   = ylen*get_len(self.level)
                self.zlen   = zlen*get_len(self.level)
            if ('s' not in self.hierarchy): # stone size # TO DO - resize the existing grains
                self.mu=mu
            else:
                self.mu=get_mean(self.level) # adapt grain size to the level
            if ('n' not in self.hierarchy):
                self.n_grains=ngrains
            else:
                self.n_grains=get_n(self.level) # adapt the number of grains to the level
        else:
            self.update=True
            self.level=update
            if ('d' in self.hierarchy): # domain size
                self.xlen   = xlen*get_len(self.level)
                self.ylen   = ylen*get_len(self.level)
                self.zlen   = zlen*get_len(self.level)
            if ('s' in self.hierarchy): # stone size
                self.mu     = get_mean(self.level) # adapt grain size to the level
            if ('n' in self.hierarchy): # stone numbers
                self.n_grains     = get_n(self.level) # adapt the number of grains to the level
        if ('p' not in self.hierarchy):
            self.max_tries   = max_tries
        else:
            self.max_tries   = get_max_tries(self.level)
        #self.sample(update)
    # --------------------------

    # ---------- RANDOM SAMPLING
    # when update>0,None it keeps the previous sample and add new grains
    #@jit
    def sample(self):
        if periodic:
            detbc  = 0
        else:
            detbc = detachedbc*2-1
        print(detbc)
        vspace = 2*voidspace
        if (psd=="lognormal"):
            lmu=math.log(self.mu)
            ls=math.sqrt(math.log(self.cov*self.cov+1.0))
            mu_eq=self.mu*exp(ls**2*3/2)
            def sample_mu(): #lognormal distribution
                return random.lognormvariate(lmu,ls)
        elif (psd=="constant"):
            mu_eq=((self.mu*(1+self.cov))**4-(self.mu*(1-self.cov))**4)**(1/3)
            def sample_mu(): #lognormal distribution
                return self.mu
        elif (psd=="uniform"):
            mu_eq=self.mu
            def sample_mu(): #uniform distribution
                return random.uniform(self.mu-self.cov*mu,self.mu+self.cov*self.mu)
        vol0   = ((self.xlen-vspace*mu_eq*4/3)-4/3*detbc*mu_eq)*(self.ylen-4/3*detbc*mu_eq)*(self.zlen-4/3*detbc*mu_eq)   # total volume
        #vol0   = ((self.xlen-voidspace*mu_eq*4/3))*(self.ylen)*(self.zlen)   # total volume

        #@jit
        def sample_x(r=0):
            shift = ((self.xlen*0.5-vspace*r)-detbc*r)
            return (random.random()*2*shift-shift)

        #@jit
        def sample_y(r=0):
            return random.random()*(self.ylen-2*detbc*r)-(self.ylen*0.5-r*detbc)

        #@jit
        def sample_z(r=0):
            return random.random()*(self.zlen-2*detbc*r)-(self.zlen*0.5-r*detbc)

        def multiple_region_update(xyz,r=0): # works only without vspace or detbc or detached
            if xyz[0]<(self.xlen-detbc*r*2)*(pow(region_ratio,2)/(pow(region_ratio,2)+one)-half):
                xyz[0]=-abs(sample_x())*two/(region_ratio+1)-(self.xlen-detbc*r*2)*(half-one/(region_ratio+1))
                return (one+one/pow(region_ratio,2))/(one+one/pow(region_ratio,1))
            else:
                xyz[0]=abs(sample_x())*two*region_ratio/(region_ratio+1)-(self.xlen-detbc*r*2)*(half-one/(region_ratio+1))
                return region_ratio*(one+one/pow(region_ratio,2))/(one+one/pow(region_ratio,1))

        def sample_xyz(eps, p, A, index = [-1], m = [0]):

            L = len(index)
            value = False
            flag = 0

            while value == False:

                mnew=get_ellips_radii(m)
                #generazione x,y,z casuali
                x = sample_x(mnew[0])
                y = sample_y(mnew[1])
                z = sample_z(mnew[2])

                if index == [-1]:
                    value = True

                for i in range(len(index)):
                    j = index[i]
                    # plane
                    if j == 1:
                        distance = dist_point_plane(A[i][0], A[i][1], A[i][2], A[i][3], x, y, z)
                        tmp = disc_prob(distance, eps, p)
                        if tmp == True:
                            flag = flag + 1

                    # line
                    if j == 2:
                        distance = dist_point_line(A[i][0:4], A[i+1][0:4], x, y, z)
                        tmp = cont_prob(distance, eps, p)
                        if tmp == True:
                            flag = flag + 1

                    # point
                    if j == 3:
                        distance = dist_point_point(A[i][0:4], A[i+1][0:4], A[i+2][0:4], x, y, z)
                        tmp = cont_prob(distance, eps, p)
                        if tmp == True:
                            flag = flag + 1

                    #void index
                    if j == 0:
                        flag = flag + 1

                if flag >= L:
                    value = True
                flag = 0

            return array([x,y,z])


        dfactor=1.0
        nfactor=1
        if ('d' in self.hierarchy):
            dfactor = get_len(self.tot_level-1)/get_len(self.level)
            nfactor = dfactor**dimension
        if (self.update):
            g=self.gold
            nstart=len(self.gold)
            volume=vol0-self.porosity_out*vol0
        else:
            nstart=0
            g=[]
            volume=0.0

        #------------------------------------------------
        # Insertion of random oscillations in the fractures for every sample
        for i in range(len(Matrix)):
            Matrix[i]=random_frac(PointX,PointY,PointZ)
        #-------------------------------------------------
        log=""
        ntries=0
        # ---------- GRAIN GENERATION
        g=list(g)
        p=(vol0-volume)/vol0
        self.eig=[]
        self.eif=[]
        #print("pippo")
        for i in range(nstart,self.n_grains*nfactor):
            #print(i)
            if (p<self.porosity):
                #print("Reached desired porosity "+str(p),len(g))
                break
            ntries=0
            if (ellips):
                rx = [sample_mu() for i in range(3)]
                rr = sample_ellips(rx)
                xyz=sample_xyz(EPS,Prob,Matrix,index,rr)*dfactor
                if (regions>1):
                    rr=[y*multiple_region_update(xyz,y) for y in rr]
                if (detached) and len(g)>0:
                    while inside_point(g,concatenate([xyz,rr]))>-1:
                        #rr=sample_ellips(sample_mu)
                        xyz=sample_xyz(EPS,Prob,Matrix,index,rr)*dfactor
                        ntries+=1
                        if (ntries>self.max_tries):
                            break
                vol=abs(math.pi*linalg.det(rr.reshape(3,3))*4.*third)
            else:
                rr=[sample_mu()]
                xyz=sample_xyz(EPS,Prob,Matrix,index,rr)*dfactor
                if (regions>1):
                    rr=[y*multiple_region_update(xyz,y) for y in rr]
                if (detached) and len(g)>0:
                    while inside_point(g,concatenate([xyz,rr]))>-1:
                        #rr=sample_mu()
                        xyz=sample_xyz(EPS,Prob,Matrix,index,rr)*dfactor
                        ntries+=1
                        if (ntries>self.max_tries):
                            break
                vol=math.pi*rr[0]**3*4.*third
            if (ntries>self.max_tries):
                break
            #print("placed grain", i, vol, rr, xyz, p, vol0, volume)
            g.append(concatenate([xyz,rr]))
            ## --- ESTIMATION OF POROSITY
            volume+=vol/nfactor
            #if (detached): # exact
                #volume+=vol/nfactor
            #elif (jodreytory): # statistical approximation assuming JD algorithm perfectly works
                #pcut=(vol0-volume)/(vol0-volume*min_dist**3)
                #volume+=vol*pcut/nfactor
            #else:
                #volume+=vol*p/nfactor # statistical approximation for independent spheres
            p=(vol0-volume)/vol0
        # ---------------------
        n_grains_out = len(g)
        #print(n_grains_out, " placed grains with volume and porosity ", volume, p)

        # ----POST PROCESSING
        g=array(g)
        if (jodreytory):
                ntries=run_jodreytory(g,self.xlen*dfactor,self.ylen*dfactor,self.zlen*dfactor,self.max_tries,nmoves_jd,min_dist,max_dist)
        self.ntries_out=ntries
        if ('d' in self.hierarchy):
            self.g=[]
            for i in g:
                if inside_box([self.xlen*0.5,self.ylen*0.5,self.zlen*0.5],i).all():
                    self.g.append(i)
        else:
            self.g=g
        self.g=array(self.g)

        # ----  SAVE OUTPUT VALUES
        self.gold=g
        self.n_grains_out=len(self.g)  # update the effective number of grains
        self.porosity_out=p # update the effective estimated porosity
        self.volume = volume

        if (periodic):
            self.add_symmetric_grains()

        #d = compute_dist_matrix(g)
        #ind = where(abs(d)<min_dist)
        #self.n_intersections = size(ind)/len(ind)
        # TODO compute volume intersections and estimate porosity

        print("N grains "+str(self.n_grains_out)+" with estimated porosity "+str(p))
    #---------------------------------

    # -------- INCREASE/DECREASE levels
    def levelup(self):
        geom1=copy.deepcopy(self)
        geom1.set_input(self.level+1)
        return geom1

    def leveldown(self):
        geom1=copy.deepcopy(self)
        geom1.set_input(self.level-1)
        return geom1

    # -------- add symmetric grains
    def add_symmetric_grains(self):
        g2=[]
        for i in self.g:
            bcinside = inside_box([self.xlen*0.5,self.ylen*0.5,self.zlen*0.5],i)
            for ix in range(3):
                if bcinside[ix]==False:
                    i2 = i.copy()
                    i2[ix] = i[ix]-sign(i[ix])*self.xlen
                    g2.append(i2)
        self.g=concatenate([self.g,g2])



    # -------- output geometry
    def write(self,geom_format=geom_out,destination="."):
        if (geom_out=="gmsh"):
            import gmshgetdp as gg
        elif (geom_out=="openfoam"):
            import openfoam as of
        elif (geom_out=="openscad"):
            self.write_openscad(destination)
            return
        gg.pdeproblem=""
        outg=gg.solver(self)
        outg.name=destination
        outg.snappy(self)

    # -------- output geometry
    def write_stl(self,filename="./geom"):
        self.write_openscad(filename+".scad")
        log=subprocess.check_output(['openscad', '-o',filename+".stl",filename+".scad"])
        return log.decode()


    # -------- write OpenSCAD format
    def write_openscad(self,filename="./geom.scad"):
        import solid as scad
        scad.use("geodesic_sphere.scad")
        mesh = scad.translate([-xlen/2,-ylen/2,-zlen/2])(scad.cube([xlen,ylen,zlen]))
        sphere = geodesic_sphere()
        for i,circ in enumerate(self.g):
            rx, ax, tx = rotmat_to_euler(self.get_transformation_matrix(i))
            pos = [circ[0], circ[2], circ[1]]
            rx = [rx[0],rx[1],rx[2]]
            if i==0:
                obj = scad.rotate(tx/2/pi*360,ax) \
                    (scad.scale(rx)(scad.translate(pos)(sphere)))
            else:
                obj += scad.rotate(tx/2/pi*360,ax) \
                    (scad.scale(rx)(scad.translate(pos)(sphere)))
        mesh = mesh - obj

        scad.scad_render_to_file(mesh,filename,file_header='$fn = %s;\n' % 10)


    def get_radii(self):
        if self.ellipsoid:
            return array([get_ellips_radii(i) for i in self.g])
        else:
            return array(self.g)[:,[3,3,3]]

    def get_transformation_matrix(self,i):
        return transformation_matrix(self.g[i],1)

    # ------- FUNCTION TO COMPUTE FOURIER TRANSFORM
    @jit
    def fourier(self,kx,ky,kz):
        # TODO extend for ellipsoids
        # now only the first number is used
        f=zeros_like(kx,dtype=complex)
        d=sqrt(kx**2+ky**2+kz**2)
        for i in range(len(self.g)):
            f+=exp(-2j*math.pi*(kx*self.g[i,0]+ky*self.g[i,1]+kz*self.g[i,2]))*self.g[i,3]**3*fball(d*self.g[i,3])
        return f

    @jit
    def fourier_field(self,nx,ny,nz,smooth,sl):
        #print("Begin fourier field")
        dx=self.xlen/(nx-1)
        dkx=1/dx/(nx-1)
        Lkx=dkx*nx/2
        kx=arange(-Lkx,Lkx,dkx)
        dy=self.ylen/(ny-1)
        dky=1/dy/(ny-1)
        Lky=dky*ny/2
        ky=arange(-Lky,Lky,dky)
        dz=self.zlen/(nz-1)
        dkz=1/dz/(nz-1)
        Lkz=dkz*nz/2
        kz=arange(-Lkz,Lkz,dkz)
        if constant_smoothing:
            sx=sl
            sy=sl
            sz=sl
        else:
            sx=dx*sl
            sy=dy*sl
            sz=dz*sl
        (x,y,z)=meshgrid(kx,ky,kz, sparse=False, indexing='ij')
        f=self.fourier(x,y,z)
        #f=fft.fftn(minimum(1,abs(fft.ifftn(fft.ifftshift(f)/(dx*dy*dz)))))*(dx*dy*dz)
        f=f*fourier_smooth(smooth,x,y,z,sx,sy,sz)
        f=fft.ifftshift(f)/(dx*dy*dz)
        self.ff=array(f)

    @jit
    def invert_fourier_field(self):
        self.f=abs(fft.fftshift(fft.ifftn(self.ff)))

    # ------------------------

    # ------- FUNCTION TO CREATE A FIELD FROM GRAIN POSITION
    @jit
    def create_field(self,nx,ny,nz):
        print("Begin create field")
        f=zeros((nx+1,ny+1,nz+1))
        for i in range(nx+1):
            px = -self.xlen/2 + i*self.xlen/(nx-1)
            for j in range(ny+1):
                py = -self.ylen/2 + j*self.ylen/(ny-1)
                for k in range(nz+1):
                    pz = -self.zlen/2 + k*self.zlen/(nz-1)
                    if (inside_point(g,array([px,py,pz,0]))>-1):
                        f[i,j,k]=1;
        self.f=f

    @jit
    def smooth_field(self):
        # smooth in real space, with smoothing_length cells
        print("Begin smooth field")
        for s in range(int(smoothing_length)):
            for i in range(size(self.f,0)):
                ip=min(size(self.f,0)-1,i+1)
                im=max(0,i-1)
                for j in range(size(self.f,1)):
                    jp=min(size(self.f,1)-1,j+1)
                    jm=max(0,j-1)
                    for k in range(size(self.f,2)):
                        kp=min(size(self.f,2)-1,k+1)
                        km=max(0,k-1)
                        self.f[i,j,k]=0.5*(self.f[i,j,k]+(self.f[ip,j,k]+self.f[im,j,k]+
                           self.f[i,jp,k]+self.f[i,jm,k]+self.f[i,j,kp]+self.f[i,j,km])/6)
    # ------------------------

    # ------- FUNCTION TO PRODUCE and WRITE PERMEABILITY or POROSITY FIELD

    def transform_field(self,law=upscaling_type,re=0):
        #print("Begin transform field")
        self.f=transform_fieldv(self.f,law,self.mu,re,self.porosity_out,upscaling_limiter,upscaling_constant)

    def write_field(self,path="./",f0=0,f1=1):
        # f0 and f1 are the linear coefficients of the field (perm1 and perm2)
        try:
            self.f
        except:
            self.fourier_field(int(self.resolution*self.xlen),int(self.resolution*self.ylen),int(self.resolution*self.zlen),smoothing, smoothing_length)
            self.invert_fourier_field()
            self.transform_field()
        filename=path+'porosity.3d'
        f=open(filename,'w')
        nx=size(self.f,0)
        ny=size(self.f,1)
        nz=size(self.f,2)
        f.write(str(nx)+" "+str(ny)+" "+str(nz)+"\n")
        for i in range(nx):
            px = -self.xlen/2 + i*self.xlen/(nx-1)
            for j in range(ny):
                py = -self.ylen/2 + j*self.ylen/(ny-1)
                for k in range(nz):
                    pz = -self.zlen/2 + k*self.zlen/(nz-1)
                    p=f1+(f0-f1)*self.f[i][j][k]
                    f.write(str(px)+" "+str(py)+" "+str(pz)+" ")
                    f.write(str(p)+" "+str(p)+" "+str(p)+"\n")
        f.close()
        #print("End write field")
    # ------------------------

#----------------------------- END CLASS GEOM

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
