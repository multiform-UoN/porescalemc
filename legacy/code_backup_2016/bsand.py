# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Gianluca Boccardo & Matteo Icardi
# March 2015
#--------------------------------------------------------------------------------------------
# module for generating random geometries with Bsand (Blender + Bullet Physics)

from randomgeo import *
geom_base = geom # rename base class

from bsandDefault import *

try:
    exec(open("bsandDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file bsandDict.py not found or not readable. Loaded default values")

try:
    exec(open("runDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file runDict.py not found or not readable. Using module specific options")

# this works only if bsand.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("bsandDict.py").read(),globals_)


import bpy
from mesh_volume_tools import *

if margin<=0:
    usemargin=False
else:
    usemargin=True

####################BLENDER FUNCTIONS####################
def fix_context():
    #bpy.context.area.type = 'VIEW_3D'
##Fix bpy.context if some command (like .blend import) changed/emptied it"""
    #for window in bpy.context.window_manager.windows:
        #print('win',window)
        #screen = window.screen
        #for area in screen.areas:
            #print('area',area)
            #if area.type == 'VIEW_3D':
                #for region in area.regions:
                    #print('region',region)
                    #if region.type == 'WINDOW':
                        #override = {'window': window, 'screen': screen, 'area': area, 'region': region}
                        #print('override')
                        #bpy.ops.screen.screen_full_area(override)
                        #break
    return


def deselect_all():
    for i in bpy.context.selected_objects:
        bpy.context.selected_objects[0].select=False
        #it looks weird with respect to select_all below, as in it refers to [0] every time,
        #but remember each loop cycle [0] gets DEselected, so each cycle [0] is actually a new object


def select_all(pattern=""):
    for i in bpy.context.selectable_objects:
        if pattern in i.name:
            i.select=True
        

####################CONTAINER DEFINITION####################
def box_create(sob,remove_top=True,name="Container",shift=0):
    if not remove_top:
        bpy.ops.mesh.primitive_cube_add()
        bpy.context.active_object.scale[0]=sob[0]*.5
        bpy.context.active_object.scale[1]=sob[1]*.5
        bpy.context.active_object.scale[2]=sob[2]*.5
        obj=bpy.context.object
        obj.name=name+"Cube"
    else:
        bpy.ops.mesh.primitive_plane_add(location=(0,0,-sob[2]*.5+shift), rotation=(0,0,0))
        bpy.ops.rigidbody.object_add(type='PASSIVE')
        obj=bpy.context.object
        obj.name=name+"MinZ"
        obj.scale[0]=(sob[0]*.5)
        obj.scale[1]=(sob[1]*.5)
        obj.rigid_body.restitution=restitution
        #obj.rigid_body.use_margin=usemargin
        #obj.rigid_body.collision_margin=margin
        obj.rigid_body.friction=1.0
        obj.rigid_body.collision_shape='BOX'
        #MinX side
        bpy.ops.mesh.primitive_plane_add(location=(-sob[0]*.5,0,shift), rotation=(0,pi*.5,0))
        bpy.ops.rigidbody.object_add(type='PASSIVE')
        obj=bpy.context.object
        obj.name=name+"MinX"
        obj.scale[0]=(sob[2]*.5)
        obj.scale[1]=(sob[1]*.5)
        obj.rigid_body.restitution=restitution
        #obj.rigid_body.use_margin=usemargin
        #obj.rigid_body.collision_margin=margin
        obj.rigid_body.friction=1.0
        obj.rigid_body.collision_shape='BOX'
        #MaxX side
        bpy.ops.mesh.primitive_plane_add(location=(sob[0]*.5,0,shift), rotation=(0,pi*.5,0))
        bpy.ops.rigidbody.object_add(type='PASSIVE')
        obj=bpy.context.object
        obj.name=name+"MaxX"
        obj.scale[0]=(sob[2]*.5)
        obj.scale[1]=(sob[1]*.5)
        obj.rigid_body.restitution=restitution
        #obj.rigid_body.use_margin=usemargin
        #obj.rigid_body.collision_margin=margin
        obj.rigid_body.friction=1.0
        obj.rigid_body.collision_shape='BOX'
        #MinY side
        bpy.ops.mesh.primitive_plane_add(location=(0,-sob[1]*.5,shift), rotation=(pi*.5,0,0))
        bpy.ops.rigidbody.object_add(type='PASSIVE')
        obj=bpy.context.object
        obj.name=name+"MinY"
        obj.scale[0]=(sob[0]*.5)
        obj.scale[1]=(sob[2]*.5)
        obj.rigid_body.restitution=restitution
        #obj.rigid_body.use_margin=usemargin
        #obj.rigid_body.collision_margin=margin
        obj.rigid_body.friction=1.0
        obj.rigid_body.collision_shape='BOX'
        #MaxY side
        bpy.ops.mesh.primitive_plane_add(location=(0,sob[1]*.5,shift), rotation=(pi*.5,0,0))
        bpy.ops.rigidbody.object_add(type='PASSIVE')
        obj=bpy.context.object
        obj.name=name+"MaxY"
        obj.scale[0]=(sob[0]*.5)
        obj.scale[1]=(sob[2]*.5)
        obj.rigid_body.restitution=restitution
        #obj.rigid_body.use_margin=usemargin
        #obj.rigid_body.collision_margin=margin
        obj.rigid_body.friction=1.0
        obj.rigid_body.collision_shape='BOX'

def cylinder_create(d,h,remove_top=True,remove_bottom=False,name="Container",shift=0):
    # TODO removing both top and bottom does not work. it makes everything disappear
    bpy.ops.mesh.primitive_cylinder_add(radius=d*.5,depth=h,location=(0,0,shift))
    bpy.context.object.name=name
    if not remove_top and not remove_bottom:
        return
    obj=bpy.context.active_object #pointer to active object
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)      # set to Edit Mode
    #EDIT MODE ON 
    bpy.ops.mesh.select_all(action='DESELECT') #deselect all, just to be sure
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    #OBJECT MODE ON 
    if remove_top:
        obj.data.polygons[-4].select = True #Select only those to delete (-4 is top face, -1 is bottom face)
    elif remove_bottom:
        obj.data.polygons[-1].select = True #Select only those to delete (-4 is top face, -1 is bottom face)
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)      # set to Edit Mode AGAIN
    #EDIT MODE ON
    bpy.ops.mesh.delete(type='FACE') #Finally delete the faces
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    #OBJECT MODE ON 
    bpy.ops.rigidbody.object_add(type='PASSIVE')
    obj.rigid_body.collision_shape='MESH'
    #obj.rigid_body.use_margin=usemargin
    #obj.rigid_body.collision_margin=margin
    obj.rigid_body.friction=1.0
####################CONTAINER DEFINITION END####################


# ------------ CLASS DEFINITION FOR A GEOMETRY
class geom(geom_base):
    # ---------- GRAIN GENERATION - CLASS INITIALIZATION
    # see randomgeo parent class


    # ----------- SET GEOMETRY INPUT
    # this function is called to initialize the sample (with update =0)
    # and to keep the same sample and solve it at different levels (update=1)
    def set_input(self,update=None):
        #bpy.ops.wm.read_homefile()
        self.stl = True # Tell openfoam to use STL file
        geom_base.set_input(self,update) # call parent's method
        if "cyl" in container:
            #self.cyl_diam   = sqrt(self.ylen**2+self.zlen**2) # cylinder outside box
            self.cyl_diam   = min(self.ylen,self.zlen)  # cylinder inside box
            self.cyl_height = self.xlen+cyl_extra_height
        elif "box" in container:
            self.box_side   = [self.ylen+box_extra_length, self.zlen+box_extra_length, self.xlen+box_extra_length]
        if 'r' in self.hierarchy:
            self.refine_level = refine_level + self.level
        else:
            self.refine_level  = refine_level
        self.xlen=self.xlen*deposition_height
        self.grain_rescale=set_grain_rescale(self.level)

        ############### REMOVE ALL EXISTING OBJECTS and RESET
        bpy.context.scene.frame_set(0)
        deselect_all()
        for obj in bpy.data.objects:
            obj.select=True
            bpy.ops.object.delete(use_global=True)

        ############### CONTAINER CREATION ###############
        if "box" in container:
            box_create(self.box_side)
            # second "safety" layer
            box_create([self.box_side[0]+self.mu,self.box_side[1]+self.mu,self.box_side[2]*3+self.mu],name="Container2",shift=self.box_side[2])
        else: #cylinder
            cylinder_create(self.cyl_diam,self.cyl_height)
            # second "safety" layer
            cylinder_create(self.cyl_diam+self.mu,self.cyl_height*3+self.mu,name="Container2",shift=self.cyl_height)


        ############### SET WORLD PARAMETERS
        bpy.data.scenes["Scene"].gravity[2]=-abs(gravity)
        bpy.data.scenes["Scene"].rigidbody_world.steps_per_second=set_nsteps(self.level)
        bpy.data.scenes["Scene"].rigidbody_world.solver_iterations=set_niter(self.level)
    # --------------------------


    # ---------- RANDOM SAMPLING
    # see randomgeo parent class
    def sample(self):
        geom_base.sample(self)
        self.xlen=self.xlen/deposition_height
        self.blender_run()
        if postprocess:
            self.blender_post()
        if write_final_file:
            self.write("blender") # no need to write it always

    # ---------- BLENDER SETUP AND SIMULATION
    def blender_run(self):
        ####################NONPRIMITIVE MODEL IMPORT####################
        if non_primitive_model==0:
            bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=self.refine_level, size=1., location=(0,0,0), rotation=(0, 0, 0))
            model_name="Icosphere"
            primitive=1
        else:
            model_name=non_primitive_model
            cwd=bpy.path.abspath("//")
            structure=cwd+"Library/"
            bpy.ops.wm.link_append(filepath=structure+model_name+".blend/Object/"+model_name,directory=structure+model_name+".blend/Object/",filename=model_name,link=False)
            bpy.ops.object.make_local(type='SELECT_OBJECT')
            bpy.data.objects[model_name].select=True
            bpy.ops.transform.translate(value=(0,0,-10))
            for i in range(len(list(bpy.data.objects))):
                if bpy.data.objects[i].name == model_name:
                    i_index=i
        obj=bpy.data.objects[model_name]
        bpy.ops.rigidbody.object_add(type='ACTIVE')
        obj.game.physics_type = 'RIGID_BODY'
        obj.rigid_body.mass=mean_mass*self.mu**3
        obj.rigid_body.friction=friction
        obj.rigid_body.linear_damping=linear_damping
        obj.rigid_body.angular_damping=angular_damping
        obj.rigid_body.restitution=restitution
        obj.rigid_body.use_margin=usemargin
        obj.rigid_body.collision_margin=margin
        obj.rigid_body.use_deactivation = True
        #obj.rigid_body.use_start_deactivated = True
        obj.rigid_body.collision_shape='SPHERE' #TODO change for objects other than spheres

        #################### MODEL PLACEMENT ####################
        xshift=(deposition_shift)*self.xlen
        sce = bpy.context.scene
        obs = []
        for namer,grain in enumerate(self.g):
            #Add the icosphere
            copy = obj.copy()
            copy.location = [grain[1],grain[2],grain[0]+xshift]
            copy.scale = [grain[3],grain[3],grain[3]]
            #copy.parent_set(obj)
            #copy.data = copy.data.copy() # also duplicate mesh, remove for linked duplicate
            name='Grain'+str(namer)
            copy.name=name
            #bpy.ops.object.transform_apply(location=True, rotation=False, scale=True)
            obs.append(copy)
        for ob in obs:
            sce.objects.link(ob)
        sce.update()
        deselect_all()
        select_all("Grain")
        bpy.ops.rigidbody.objects_add()
                    #bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
                    #bpy.ops.transform.rotate(value=random.randint(-2,2), axis=(1.0, 1.0, 0.5))
                    #obj.rigid_body.mass=mean_mass*obj.scale[2]**3
                    #bpy.context.object.rigid_body.use_deactivation = True
                    #bpy.context.object.rigid_body.use_start_deactivated = True
        deselect_all()
        bpy.data.objects[model_name].select=True #Select original model
        bpy.ops.object.delete()
        ####################MODEL PLACEMENT END####################

        ####################BLENDER ANIMATION SETTINGS####################
        if autotime:
            if (linear_damping<0.5):
                fall_time = 2./abs(gravity)*sqrt((deposition_shift+deposition_height)*self.xlen)
                fin_vel = fall_time*abs(gravity)
            else:
                fin_vel = abs(gravity)/linear_damping
                fall_time = (deposition_shift+deposition_height)*self.xlen/fin_vel
            tot_time  = btime*(fall_time + fin_vel*mean_mass*relax_time)
        else:
            tot_time = max_time
        #print(fin_vel, fall_time, tot_time/btime)
        bpy.context.scene.frame_end = tot_time
        bpy.data.scenes['Scene'].rigidbody_world.point_cache.frame_end= tot_time
        ####################BLENDER ANIMATION SETTINGS END####################

        ####################COMMAND LINE AUTOMATION####################
        if write_blender_file:
            bpy.ops.wm.save_as_mainfile(filepath=workdir+self.name+"b1-start.blend")
        cnt=1
        while (cnt<=bpy.context.scene.frame_end and deposition):
            bpy.context.scene.frame_set(cnt)
            cnt=cnt+1
        #bpy.ops.ptcache.bake_all()
        #bpy.context.scene.frame_set(tot_time)
        if write_blender_file:
            bpy.ops.wm.save_as_mainfile(filepath=workdir+self.name+"b2-endSimulation.blend")
        ####################COMMAND LINE AUTOMATION END####################

        #################### EXPORT and RESCALE
        deselect_all()
        for obj in bpy.data.objects:
            if "Container" in obj.name:
                obj.select=True
                bpy.ops.object.delete()
        for i in bpy.data.objects:
            if "Grain" in i.name:
                i.scale=i.scale*self.grain_rescale
        if write_blender_file:
            bpy.ops.wm.save_as_mainfile(filepath=workdir+self.name+"b3-forExport.blend")
        #################### EXPORT END
    #---------------------------------

    def blender_post(self):
        deselect_all()
        maxd=self.xlen/2+box_extra_length*(container=="box")+cyl_extra_height*(container=="cyl")
        i=0
        iname=0
        g=list(self.g)
        while i<len(g):
            gname='Grain'+str(iname)
            locx=bpy.data.objects[gname].matrix_world.translation[0]
            locy=bpy.data.objects[gname].matrix_world.translation[1]
            locz=bpy.data.objects[gname].matrix_world.translation[2]
            if update_grain:
                g[i][0]=locz
                g[i][1]=locx
                g[i][2]=locy
            #print(gname, locz, self.xlen/2)
            #Deleting grains whose z-location is lower than the base of the container
            if abs(locz)>maxd and remove_grain:
                bpy.data.objects[gname].select=True
                bpy.ops.object.delete()
                #print("removed grain",gname,i)
                deselect_all()
                if update_grain:
                    del g[i]
                else:
                    i+=1
                iname+=1
                continue
            i+=1
        self.g=array(g)
        # Here ends the loop over all the original grains. But some have been deleted! (the Fallen Ones)
        # Let's create an array of remaining grains
        grainArray=[]
        for item in bpy.data.objects:
            if item.name.find('Grain')!=-1:
                grainArray.append(item)
        maxHeight=max([obj.matrix_world.translation[2] for obj in bpy.data.objects])
        minHeight=min([obj.matrix_world.translation[2] for obj in bpy.data.objects])
        deselect_all()
        ### WHY CREATING A NEW CYLINDER? BETTER TO USE THE SAME CYLINDER USED FOR THE SIMULATION (even if it's not totally filled. A new "measuring" REV can be defined in the parameters #TODO
        ############### CONTAINER CREATION(again) ###############
        if "box" in container:
            box_create([self.ylen*measure_volume,self.zlen*measure_volume,self.xlen*measure_volume],False)
        else: #cylinder
            cylinder_create(min(self.zlen,self.ylen)*measure_volume,self.xlen*measure_volume,False,False)
        #bpy.ops.mesh.primitive_cylinder_add(radius=self.cyl_diam*0.5,depth=(maxHeight-minHeight),location=(0,0,((maxHeight+minHeight)/2)))
        justCreated=bpy.context.active_object
        justCreated.name='Bulk'   
        # POROSITY CALC PART 1 @@@@@@@@@@@@@@@@
        bpy.context.scene.objects.active = bpy.data.objects['Bulk']
        bpy.ops.object.volume()
        totVol=bpy.data.objects['Bulk']['volume']
        #print('Total volume is',totVol)
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        # Select all grains!
        deselect_all()
        sumVol=0
        sumVol0=0
        for item in grainArray:
            item.select=True
            bpy.context.scene.objects.active = item
            bpy.ops.object.volume()
            sumVol+=item["volume"]
            sumVol0+=prod(item.scale)*4./3.*pi
        # The obj selected just now becomes active
        # That's because if the active object isn't one of the selected objects to join, the join operation fails.
        # Ah, the joys of Blender scripting!
        bpy.context.scene.objects.active = item
        # Join all selected, but first we need to override the context to join also without X-server!
        override=bpy.context.copy()
        for i in bpy.context.window_manager.windows:
            override['window']=i
            break
        scene = bpy.context.scene
        # one of the objects to join
        override['active_object'] = grainArray[0]
        override['selected_objects'] = grainArray
        # we need the scene bases as well for joining
        override['selected_editable_bases'] = [scene.object_bases[ob.name] for ob in grainArray]
        if write_blender_file:
            bpy.ops.wm.save_as_mainfile(override,filepath=workdir+self.name+"b4-beforeJoin.blend")
        bpy.ops.object.join(override)
        justCreated=grainArray[0]#override.active_object
        justCreated.name='GrainBed' 
        if write_blender_file:
            bpy.ops.wm.save_as_mainfile(filepath=workdir+self.name+"b5-beforeBoolean.blend")
        #Next two lines: big cyl objA has a boolean modifier mymod, whose item is grain bed objB
        objA=bpy.data.objects['Bulk']
        objB=bpy.data.objects['GrainBed']
        mymod=objA.modifiers.new("Moddy",'BOOLEAN')
        mymod.object=objB
        #And whose operation is DIFFERENCE
        mymod.operation='DIFFERENCE'
        #Now, objA is made active because if not fuck if Blender will work, as usual.
        bpy.context.scene.objects.active = objA
        #And now Moddy is applied (Moddy is the name of mymod) (while we have objA active)
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Moddy")
        # POROSITY CALC PART 2 @@@@@@@@@@@@@@@@
        bpy.context.scene.objects.active = bpy.data.objects['Bulk']
        bpy.ops.object.volume()
        voidVol=bpy.data.objects['Bulk']['volume']
        varepsilon=(voidVol/totVol)
        if write_blender_file:
            bpy.ops.wm.save_as_mainfile(filepath=workdir+self.name+"b6-postProcessed.blend")
        # Delete the cylinder to run openfoam simulation
        deselect_all()
        obj=bpy.data.objects['Bulk']
        obj.select=True
        bpy.ops.object.delete() 
        deselect_all()

        #print('Porosity is',varepsilon,voidVol,totVol-sumVol,totVol-sumVol0)
        self.gold=self.g
        self.n_grains_out=len(self.g)  # update the effective number of grains
        self.porosity_out=varepsilon # update the effective estimated porosity
        self.volume=totVol
        return

    # -------- output geometry
    def write(self,geom_out="blender"):
        geom_base.write(self,geom_out)
        if (geom_out=="blender"):
            bpy.ops.wm.save_as_mainfile(filepath=workdir+self.name+".blend")

    # -------- output geometry
    def write_stl(self,filename):
        deselect_all()
        for obj in bpy.data.objects:
            if "Container" in obj.name:
                obj.select=True
                bpy.ops.object.delete()
        if container=="cyl":
            cylinder_create(self.cyl_diam,self.cyl_height*3,False,False)
        select_all()
        bpy.ops.export_mesh.stl(filepath=filename+"grain.stl")
        return ""


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
