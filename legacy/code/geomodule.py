# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module for geometrical operations with grain packings (used by bsand and randomgeo)


try:
    random
except:
    from general import *


from randomgeoDefault import *

try:
    exec(open("randomgeoDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file randomgeoDict.py not found or not readable. Loaded default values")


# this works only if random.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("randomgeoDict.py").read(),globals_)


#from scipy.spatial import KDTree
from scipy.spatial import distance


# --------- transform rotation matrix to deform euler axis and angle
def rotmat_to_euler(m):
    mm = reshape(m,[3,3])
    rr = sqrt(diag(mm.T.dot(mm)))
    a = mm.dot(diag(1/rr))
    theta = arccos(0.5*trace(a)-.5)
    ax = [(a[2,1]-a[1,2])/(2*sin(theta)), (a[0,2]-a[2,0])/(2*sin(theta)), (a[1,0]-a[0,1])/(2*sin(theta)) ]
    return rr, ax, theta


# --------- compute distance matrix
@jit
def compute_dist_matrix(g):
    #if ellips:
    return nan_to_num(distance.pdist(g,ellipsdist))
    #return squeeze(array(distance.pdist(g[:,:3])/compute_radius_matrix(g[:,3])))

# --------- compute distance matrix
@jit
def compute_dist_pair(g1,g2):
    #if ellips:
        return nan_to_num(distance.cdist(g1,g2,ellipsdist))
    #return squeeze(distance.cdist(g1[:,:3],g2[:,:3])/compute_radius_pair(g1[:,3],g2[:,3]))

# ---------- compute radius matrix (sum of radii of grains
@jit
def compute_radius_matrix(r):
    #n=length(r)
    #return outer(ones(n),r))+vectorform(outer(r,ones(n))
    #TODO ellipsoids (see function below)
    return vectorform(matrix(r)+transpose(matrix(r)))

# ---------- compute radius matrix (sum of radii of grains)
@jit
def compute_radius_pair(r1,r2):
    #n=length(r)
    #return outer(ones(n),r))+vectorform(outer(r,ones(n))
    #TODO ellipsoids (see function below)
    return sum(meshgrid(r2,r1),0)

@jit
def transformation_matrix(x,b):
# return the transformation matrix from the vector of generalized radii
    if len(x)<5:
        return diag([x[-1],x[-1],x[-1]]).dot(b)
    #m=zeros([3,3])
    #m[transpose(tril(ones([3,3],dtype=bool)))]=x[-6:]
    #m[tril(ones([3,3],dtype=bool))]=x[-6:]
    return x[-9:].reshape(3,3).dot(b)

@jit
def stretching_length(x,b):
# return the stretching in one direction (b must be unitary)
    if len(x)<4:
        return 0
    if len(x)==4:
        return sqrt(((ones(3)*x[-1]*b)**2).sum(0))
    #m=zeros([3,3])
    #m[transpose(tril(ones([3,3],dtype=bool)))]=x[-6:]
    #m[tril(ones([3,3],dtype=bool))]=x[-6:]
    r=get_ellips_radii(x)
    return sqrt(((x[-9:].reshape(3,3).dot(eye(3)/r).T.dot(b)*r)**2).sum(0))

@jit
def get_ellips_radii(x):
    return sqrt((transformation_matrix(x,eye(3))**2).sum(0))

# ------- FUNCTION TO FIND IF AN ELLIPSOID or POINT TOUCHES A ELLIPSOID
@jit
def inside_point(g,g0):
    # return the index of the closest grain
    ##dist=transpose([a,b,c])-[a0,b0,c0]
    ##dist=distance.cdist([[0,0,0]],dist)[0]/compute_radius_pair(d,d0,dist)
    #dist=distance.cdist(array([g0])[:,:3],array(g)[:,:3])[0]
    #dist=dist/distance.cdist(array([g0]),array(g),ellipsdist_correction)[0]
    #dist=nan_to_num(dist)
    #inside=argmin(dist)
    #if (dist[inside]<=min_dist):
        #return inside
    #return -1
    inside=-1
    for i,k in enumerate(g):
        if ellipsdist(k,g0)<min_dist:
            return i
    return inside
# ------------------------

# Sample a rotation matrix
def sample_ellips(rr):
    # How to generate a random unitary matrix (Maris Ozols, 2009)
    # It has to preserve the Haar measure
    q,r = linalg.qr([[random.normalvariate(0,1) for i in range(3)] for j in range(3)])
    eif = q.dot(diag(diag(r)/abs(diag(r))))
    eif = eif/linalg.det(eif)
    return ravel(eif.dot(diag(array(rr))))

# dummy sample radius
def sample_one():
    return 1.0

@jit
def ellipsdist_eucl(x,y):
    # compute relative distance between ellipsoids x=(c1,m1) y=(c2,m2) with c=centers m=transf matrices
    d=array(x[:3])-array(y[:3])
    n=sqrt((d**2).sum())
    return n/(stretching_length(x,d/n)+stretching_length(y,d/n))

@jit
def ellipsdist_per(x,y):
    # compute relative distance between ellipsoids x=(c1,m1) y=(c2,m2) with c=centers m=transf matrices
    d=array([modsign(x[0]-y[0],xlen/2),modsign(x[1]-y[1],ylen/2),modsign(x[2]-y[2],zlen/2)])
    n=sqrt((d**2).sum())
    return n/(stretching_length(x,d/n)+stretching_length(y,d/n))

if periodic:
    ellipsdist=ellipsdist_per
else:
    ellipsdist=ellipsdist_eucl

@jit
def pointdist_eucl(x,y):
    # compute distance vector between points
    return array(x[:3])-array(y[:3])

@jit
def pointdist_per(x,y):
    return array([modsign(x[0]-y[0],xlen/2),modsign(x[1]-y[1],ylen/2),modsign(x[2]-y[2],zlen/2)])

if periodic:
    pointdist=pointdist_per
else:
    pointdist=pointdist_eucl

@jit
def ellipsdist_correction(x,y):
    # compute relative distance between ellipsoids x=(c1,m1) y=(c2,m2) with c=centers m=transf matrices
    d=array(x[:3])-array(y[:3])
    n=sqrt((d**2).sum())
    return stretching_length(x,d/n)+stretching_length(y,d/n)


# ------- FUNCTION TO FIND IF AN ELLIPSOID or SPHERES TOUCHES THE CONTAINER
# exact for spheres
# ellipsoid only for non-rotated ones
#@jit
def inside_box(box,m):
    rx=[stretching_length(m,[1,0,0]),stretching_length(m,[0,1,0]),stretching_length(m,[0,0,1])]
    inside=[]
    for i in range(3):
        inside+=[abs(m[i])+abs(rx[i])<box[i]]
    return array(inside)
# ------------------------

# ------- FUNCTION TO COMPUTE VOLUME OF THE SPHERE INTERSECTIONS
# exact for spheres, approximation for ellipsoids
@jit
def sphere_intersection(dist,rr,dd):
    return math.pi*(dd+rr-dist)**2*(dist**2+2*dist*rr-3*rr**2+2*dist*dd+6*rr*dd-3*dd**2)/(12*dist)
# ------------------------

## --------- compute overlapping volume (only using pairs, volume is overestimated)
# TODO transformation matrix can be used for each distance and , assuming spherical cap, an estimate can be obtained also for ellipsoids
#def compute_overlapping_volume(g):
    #n_grains_out = len(self.x)
    #vol=0
    #for i in range(n_grains_out-1):
            #dist=inside_point(self.x[i+1:],self.y[i+1:],self.z[i+1:],self.r[i+1:],self.x[i],self.y[i],self.z[i],self.r[i],'all')
            #for j in where(dist<1)[0]:
                #vol+=sum(sphere_intersection(self.x[j+i+1],self.y[j+i+1],self.z[j+i+1],self.r[j+i+1],self.x[i],self.y[i],self.z[i],self.r[i],dist[j]))
    #return vol


# ------ FUNCTIONS FOR FRACTURES (Nathan)
# Functions to create random fractures -> Output: Matrix and index vector
# Main idea: from three points passes one and one only plane
#@jit
def plane(xlen = 2,ylen = 1,zlen = 1):
    x = []
    y = []
    z = []
    #for i in range(3):
        #x.append(xlen*random.random()-xlen/2)
        #y.append(ylen*random.random()-ylen/2)
        #z.append(zlen*random.random()-zlen/2)
    # TODO  ----- DEBUG FOR FIXED PLANE Z=0
    x = [ 0, 0.1, 0.2]
    y = [ 0, 0.1, 0]
    z = [ 0, 0, 0]
    A = (y[1]-y[0])*(z[2]-z[0])-(z[1]-z[0])*(y[2]-y[0])
    B = -(x[1]-x[0])*(z[2]-z[0])+(z[1]-z[0])*(x[2]-x[0])
    C = (x[1]-x[0])*(y[2]-y[0])-(y[1]-y[0])*(x[2]-x[0])
    a = A
    b = B
    c = C
    d = -(A*x[0]+B*y[0]+C*z[0])
    return a,b,c,d,x,y,z

#@jit
def build_matrix(nplane = 1,nline = 0,npoint = 0,xlen = 2,ylen = 1,zlen = 1):
    A = []
    index = []
    x = []; y = []; z = [];
    #generate planes
    for i in range(nplane):
        Plane = plane()
        A.append(Plane[0:4])
        index.append(1)
        x.append(Plane[4])
        y.append(Plane[5])
        z.append(Plane[6])
    #generate lines
    for i in range(nline):
        A.append(plane())
        A.append(plane())
        index.append(2)
        index.append(0)
    #generate points
    for i in range(npoint):
        A.append(plane())
        A.append(plane())
        A.append(plane())
        index.append(3)
        index.append(0)
        index.append(0)
    return A, index, x, y, z

# Function to generate random oscillations in the fracture every time the
# generates a sample
@jit
def random_frac(x,y,z,zlen = 1):
    x = x[0]; y = y[0]; z = z[0]
    x = array(x); y = array(y); z = array(z);
    x = x + x*(randomize_frac*random.random()-randomize_frac/2)
    y = y + y*(randomize_frac*random.random()-randomize_frac/2)
    z = z + z*(randomize_frac*random.random()-randomize_frac/2)
    A = (y[1]-y[0])*(z[2]-z[0])-(z[1]-z[0])*(y[2]-y[0])
    B = -(x[1]-x[0])*(z[2]-z[0])+(z[1]-z[0])*(x[2]-x[0])
    C = (x[1]-x[0])*(y[2]-y[0])-(y[1]-y[0])*(x[2]-x[0])
    a = A
    b = B
    c = C
    d = -(A*x[0]+B*y[0]+C*z[0])
    return a,b,c,d



# The following functions compute the distance from the random point (x,y,z) to the element of the
# matrix (plane, line, point)
@jit
def dist_point_plane(a, b, c, d, x, y, z):
    num = a*x+b*y+c*z+d
    num = abs(num)
    den = math.sqrt(a**2+b**2+c**2)
    distance = num/den
    return distance

@jit
def dist_point_line(A, B, x, y, z):
    # I find the plane perpendicular to the line and passing for the point
    a = A[1]*B[2]-A[2]*B[1]
    b = -(A[0]*B[2]-A[2]*B[0])
    c = A[0]*B[1]-A[1]*B[0]
    d = -(a*x+b*y+c*z)
    # I find the projection of the point in the line
    M = [A[0:3],B[0:3],[a,b,c]]
    n = [-A[3],-B[3],-d]
    point = linalg.solve(M,n)
    # I compute the distance point-line
    distance = math.sqrt((point[0]-x)**2+(point[1]-y)**2+(point[2]-z)**2)
    return distance

@jit
def dist_point_point(A, B, C, x, y, z):
    M=[A[0:3],B[0:3],C[0:3]]
    n=[-A[3],-B[3],-C[3]]
    point = linalg.solve(M,n)
    # I compute the distance point-point
    distance = math.sqrt((point[0]-x)**2+(point[1]-y)**2+(point[2]-z)**2)
    return distance

#--------------------------------------------------------------------------------------
# Functions to transform porosity to the inverse of permeability
def transform_fieldv(field,law,mu,re=0,porglob=0.5,invkg=invkglob,upscaling_limiter="none",upscaling_constant=1):
    eps=1-maximum(minimum(field,1-1e-3),1e-3)
    if (law=="kozeny"): # valid only for porosity < 0.4
        if (upscaling_limiter in ("maxmin","minmax")):
            eps0=max(0.08,min(1-porglob,porglob))
            eps0=1.3*eps0**2+0.26*eps0-0.0042 # minimum porosity to ensure consistent homogenization (harmonic average of inverse permeability of 0-1 field equal to effective inverse permeability computed with smooth field)
            eps1=1-eps0
        elif (upscaling_limiter=="min"):
            eps1=min(max(0.05,porglob),1)
            eps1=-0.79*eps1**3+1.5*eps1**2+0.3*eps1-0.012
            eps0=1e-5
        else:
            eps0=1e-5
            eps1=1-eps0
        eps=minimum(eps1,maximum(eps,eps0))
        return (180*(1-eps)**2)/(eps**2*mu**2)/upscaling_constant
    elif (law=="poiseuille"):
        return (32*(1-eps))/(eps*mu**2)/upscaling_constant
    elif (law=="ergun"):
        return (150*(1-eps)**2/(eps**2*mu**2)+1.75*re*(1-eps)/(eps*mu**2))/upscaling_constant
    elif (law=="wenyu"):
        if (upscaling_limiter in ("maxmin","minmax")):
            eps0=max(0.08,min(porglob,1-porglob))
            eps0=0.77*eps0**2+0.53*eps0-0.0061 # minimum porosity to ensure consistent homogenization (harmonic average of inverse permeability of 0-1 field equal to effective inverse permeability computed with smooth field)
            eps1=1-eps0
        elif (upscaling_limiter=="min"):
            eps1=min(max(0.05,porglob),1)
            eps1=-0.38*eps1**3+0.84*eps1**2+0.55*eps1-0.0089
            eps0=1e-5
        else:
            eps0=1e-5
            eps1=1-eps0
        eps=minimum(eps1,maximum(eps,eps0))
        return (18*(1-eps))/eps**3.65*mu**2/upscaling_constant
    elif (law=="hill"):
        return 1/upscaling_constant; #TODO
    elif (law=="aposteriori"): # linear fit of 0 and 1/(K*porglob). these 2 values homogenize to 1/K
        return eps*invkg/porglob/upscaling_constant
    elif (law=="custom"):
        return custom_porosity_transform(eps,mu,re)/upscaling_constant
    elif (law=="mixed"): # mixed kozeny wen-yu
        if (upscaling_limiter in ("maxmin","minmax")):
            eps0=max(0.05,min(porglob,1-porglob))
            eps0=-0.85*eps0**3+2.1*eps0**2+0.09*eps0+0.0045 # minimum porosity to ensure consistent homogenization (harmonic average of inverse permeability of 0-1 field equal to effective inverse permeability computed with smooth field)
            eps1=1-eps0
        elif (upscaling_limiter=="min"):
            eps1=min(max(0.05,porglob),1)
            eps1=1.2*eps1**4-3.2*eps1**3+2.9*eps1**2+0.034*eps1-0.0017
            eps0=1e-5
        else:
            eps0=1e-5
            eps1=1-eps0
        eps=minimum(eps1,maximum(eps,eps0))
        koz=(180*(1-eps)**2)/(eps**2*mu**2)/upscaling_constant
        wen=(18*(1-eps))/eps**3.65*mu**2/upscaling_constant
        return (eps<0.3102)*koz+(eps>=0.3102)*wen
    else:
        # no upscaling, using permeabilities from dictionaries
        return eps

#--------------------------------------------------------------------------------------
# Functions to compute spherical fourier transform
@jit
def fball(rh):
#       if (rh<1e-4):
#               return 4*math.pi/3
    rh=maximum(1e-5,rh)
    rr=2*math.pi*rh
    return 2*(sin(rr)/rr**2-cos(rr)/rr)/rh
#       return 2*scipy.special.sph_jn(1,2*math.pi*rho)[1][0]/rho

# Return the fourier transform of the averaging kernels
@jit
def fourier_smooth(smooth,x,y,z,dx,dy,dz):
    if(smooth=="cell"):
        return sinc(x*dx)*sinc(y*dy)*sinc(z*dz)
    elif (smooth=="gaussian"):
        return 0.5*exp(-0.5*(math.pi**2)*((dx*x)**2+(dy*y)**2+(dz*z)**2))
    else:
        return 1

#--------------------------------------------------------------------------------------
# Functions that define a discrete and a continuos probability

@jit
def disc_prob(distance, eps, p):
    ## Discontinuos probability
    tmp = False
    if distance < eps:
        k = random.random()
        if k < p:
            tmp = True
    else:
        tmp = True
    return tmp

@jit
def cont_prob(distance, eps, p):
    ## Continuos probability
    tmp = False
    k=random.random()
    if distance <= eps/4:
        if k < p:
            tmp = True
    elif distance <= eps and distance > eps/4:
        phat = min(1,p+distance/eps)
        if k < phat:
            tmp = True
    elif distance > eps:
        tmp = True
    return tmp

    # check if the matrix of planes exists
try:
    PointX
except: # otherwise create random fractures
    tmp=build_matrix(nfractures,0,0)
    Matrix=tmp[0]
    index=tmp[1]
    PointX = tmp[2]
    PointY = tmp[3]
    PointZ = tmp[4]
else:
    # TODO -
    # Matrix= # creare la matrice dai punti
    index = ones(len(Matrix))  # TODO


#----------------------------------------------
#--- JODREY-TORY algorithm

@jit
def vectorform(m):
    n=shape(m)[0]
    return m[transpose(tril(ones([n,n],dtype=bool),-1))]

@jit
def squareform(v,dtype):
    n=int(ceil(sqrt(len(v) * 2)))
    m=zeros([n,n],dtype=dtype)
    m[transpose(tril(ones([n,n],dtype=bool),-1))]=v
    return m

@jit
def minmax(x,y):
    return min(max(x,-y),y)

if (periodic):
    move = modsign
    detbc = 0
else:
    move=minmax
    detbc=detachedbc*2-1

def run_jodreytory(g,xl,yl,zl,max_t,nmoves,mindist,maxdist):
# nmoves = number of intersections to move per each try
# xl, yl, zl length of domain
# TODO reduce the matrix from the beginning only to the couples close enough
    print("Running Jodrey-Tory algorithm with n=", len(g))
    n_grains_out = len(g)
    if n_grains_out<2:
        return 0
    ntries=0
    overlap_ind=[]
    vect2square=vectorform(transpose(mgrid[0:n_grains_out,0:n_grains_out]))
    square2vect=squareform(mgrid[0:len(vect2square)],int)
    square2vect=square2vect+transpose(square2vect)
    t0=time.time()
    dist_matrix=compute_dist_matrix(g)
    update=[]
    for ntries in range(int(max_t)):
        print("JD iter ", ntries,dist_matrix[overlap_ind])
        dist_matrix=compute_dist_matrix(g)
        detach_factor=1.
        #if len(update)>0:
            #t0=time.time()
            #update=unique(update)
            #X=g[update,:]
            #dist_matrix[ravel(square2vect[:,update])]=ravel(compute_dist_pair(g,X))
            #dist_matrix[0]=distance.pdist(g[0:2,:])
        overlap_ind=argpartition(dist_matrix,nmoves-1)[:nmoves]
        if (nanmin(dist_matrix[overlap_ind])>mindist):
            if (maxdist<0):
                break
            overlap_ind=[]
            detach_factor=-.5
            for ii in range(n_grains_out):
                mindi=square2vect[ii,argpartition(dist_matrix[square2vect[ii,:]],cluster-1)[:cluster]]
                [overlap_ind.append(i) for i in mindi[dist_matrix[mindi]>maxdist]]
                if len(overlap_ind)>=nmoves:
                    break
        overlap_ind=unique(overlap_ind)
        if len(overlap_ind)==0:
            #print("JD: converged, len(overlap_ind)=0")
            break
        update=[]
        for kk in overlap_ind:
            i,k=vect2square[kk]
            d=dist_matrix[kk]
            if d>mindist:
                continue
            rx1=stretching_length(g[i,:],[1,0,0])
            ry1=stretching_length(g[i,:],[0,1,0])
            rz1=stretching_length(g[i,:],[0,0,1])
            rx2=stretching_length(g[k,:],[1,0,0])
            ry2=stretching_length(g[k,:],[0,1,0])
            rz2=stretching_length(g[k,:],[0,0,1])
            #rand_pert=0.5+array([random.random(), random.random(), random.random()])
            dxyz = detach_factor*eps_jd*pointdist(g[i],g[k])/d
            g[i,0]=move(g[i,0]+dxyz[0],((xl*0.5-2*voidspace*rx1)-detbc*rx1))
            g[k,0]=move(g[k,0]-dxyz[0],((xl*0.5-2*voidspace*rx2)-detbc*rx2))
            g[i,1]=move(g[i,1]+dxyz[1],(yl*0.5-detbc*ry1))
            g[k,1]=move(g[k,1]-dxyz[1],(yl*0.5-detbc*ry2))
            g[i,2]=move(g[i,2]+dxyz[2],(zl*0.5-detbc*rz1))
            g[k,2]=move(g[k,2]-dxyz[2],(zl*0.5-detbc*rz2))
            update.append(i)
            update.append(k)
        if int(max_t)==ntries+1:
            print("JD: max iter reached", ntries+1)
    return ntries+1


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
