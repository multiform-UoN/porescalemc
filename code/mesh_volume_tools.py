# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####


# Copyright 2012 Tom SF Haines
# thaines@gmail.com
#
# Provides 4 operators, that may be accessed via the 'space' menu.

# *Volume:
# Calculates the volume of a well defined mesh - the output appears momentarily on the top bar, but is also written to the custom property of the selected object, in the variable 'volume'
# A well defined mesh is one with no holes, and no stray vertices/edges where the normals all point outwards and the geometry does not intersect with itself. Multiple separate chunks are acceptable as long as they do not intersect.
# Most likely users of this are people who use Blender to create meshes for 3D printers, to estimate the material required.

# *Record Volume:
# Does exactly what volume does, but calculates it for every frame in the animation, and creates an animated custom property on the object with the name 'volume' - you can then see how the property changes as the animation is run, or use it to drive something.
# Added this because I could - not really sure what it could be used for. Possibly useful to drive a gauge on some complex machine? Might have a debugging use when animating an object that should not change volume.

# *Split Volume:
# You select two objects - first a well defined mesh, then an arbitrary object that represents a plane. The plane can be any actual object - it assumes a plane in the x-y axes, which happens to match up with the plane provided by Blender. Note however that this is a true plane, i.e. infinite.
# The plane can then split the mesh in half, defining a volume above and a volume below the mesh (above is the positive z direction, below the negative z direction.). You define one of 3 custom properties on the plane object:
## below - Contains your desired volume below the plane.
## above - Contains your desired volume above the plane.
## ratio - Defines how much of the volume should be below the plane in terms of the total volume, i.e. below = ratio * total_volume.
# This operator then moves the plane the shortest distance to satisfy the constraint.
# In the event it is not possible to satisfy it the plane will end up at the extremity of the object that is closest to satisfying it.
# Main use is in testing values to find good ones to use for the below operator on the current frame.

# *Record Split Volume:
# Does the exact same thing as split volume, except it makes the adjustment for every frame of the animation range and records the objects location in each position.
# Note that if you animate the custom properties of the plane (below, above or ratio.) it will satisfy the constraint for the correct animated term.
# Basically good if you want to model something related to volume, typically where the volume is held constant as something deforms. Most obvious example is water in a deforming container (Which can simply be a rotating container.).
# Every possible deformation is supported. The animation means you could fill an object with water at a constant rate accounting for the strange shape of the object - not the same as a real water simulation, but potentially useful in cases when that would be overkill.



import math

import bpy
from mathutils import Vector, Matrix
from mathutils.geometry import area_tri, intersect_line_plane



bl_info = {
    'name': 'Volume Tools',
    'description': 'Provides an operator to provide the volume of a mesh. Also provides an operator to position a plane to satisfy a volume-based invariant relative to a mesh, i.e. the volume of the mesh above/below the plane or the ratio between the two. Animation is supported, so you can get the volume during an animation or animate a plane to maintain an invariant. This is typically used to position water planes in containers, for when a full simulation is overkill/unwanted (e.g. saline bag, which deforms when being picked up.).',
    'author': 'Tom SF Haines',
    'version': (1, 0),
    'blender': (2, 63, 0),
    'location': 'Operators: Volume, Record Volume, Split Volume, Record Split Volume',
    'category': 'Mesh'}



def dir_range(normal, mesh, mesh_to_world):
    """Given an object and a normal defining a plane at the origin returns the lowest and highest distances from the plane, e.g. the range of values for which dist in volumeSplitPlane will not have 0 in one of its return values."""
    # Initialise with the first vertex..
    dist0 = normal.dot(mesh_to_world * mesh.vertices[0].co)
    low = dist0
    high = dist0
    
    for i in range(1, len(mesh.vertices)):
        d = normal.dot(mesh_to_world * mesh.vertices[i].co)
        if d<low: low = d
        if d>high: high = d
    
    return (low, high)



def volumeSplitPlane(normal, dist, mesh, mesh_to_world):
    """Given a plane and a (clean, triangulated) mesh returns a tuple (below, above) of the volume of the plane above and below the mesh. Plane is defined in world space with a normal and distance from the origin."""
    below = 0.0
    above = 0.0
    
    def height_to_vol(h1, h2, h3, area, sign):
        # Given 3 heights, all positive or negative, updates below/above accordingly, making use of the projected triangle area...
        nonlocal below, above
        change = sign * area * (h1 + h2 + h3) / 3.0
        
        if (h1+h2+h3) > 0.0: above += change
        else: below += change
    
    def do_face(v1, v2, v3, norm):
        # Process a triangle - start by getting its vertices...
        v1 = mesh_to_world * mesh.vertices[v1].co
        v2 = mesh_to_world * mesh.vertices[v2].co
        v3 = mesh_to_world * mesh.vertices[v3].co
        
        # Calculate distances to the plane...
        h1 = normal.dot(v1) - dist
        h2 = normal.dot(v2) - dist
        h3 = normal.dot(v3) - dist
        
        # Project them back to the plane...
        p1 = v1 - normal*h1
        p2 = v2 - normal*h2
        p3 = v3 - normal*h3
        
        # Check if the face is all on one side of the plane - the simple case...
        ph1 = h1>0.0
        ph2 = h2>0.0
        ph3 = h3>0.0
        
        sign = 1.0 if normal.dot(norm)>0.0 else -1.0
        
        if ph1==ph2 and ph1==ph3:
            # Calculate the area of the projected triangle...
            area = area_tri(p1, p2, p3)
            
            # Make the update...
            height_to_vol(h1, h2, h3, area, sign)
        else:
            # The complex scenario - need to chop the triangle into 3, where each is entirly on one side of the plane, and handle each in turn...
            ## Switch around the vertices so that both v2 and v3 are on the same side...
            if ph2!=ph3:
                if ph1==ph2:
                    v1, v3 = v3, v1
                    h1, h3 = h3, h1
                    p1, p3 = p3, p1
                else: # phi1==phi3
                    v1, v2 = v2, v1
                    h1, h2 = h2, h1
                    p1, p2 = p2, p1
            
            ## Intercept the lines v1 - v2 and v1 - v3 with the plane...
            i12 = intersect_line_plane(v1, v2, normal * dist, normal)
            i13 = intersect_line_plane(v1, v3, normal * dist, normal)
            assert(i12!=None and i13!=None)
            
            ## We have defined 3 triangles - pass them each through to the height_to_vol method in turn...
            area = area_tri(p1, i12, i13)
            height_to_vol(h1, 0.0, 0.0, area, sign)
            
            area = area_tri(i12, p2, p3)
            height_to_vol(0.0, h2, h3, area, sign)
            
            area = area_tri(i12, p3, i13)
            height_to_vol(0.0, h3, 0.0, area, sign)
            
    
    # Process each triangle in the mesh in turn...
    for face in mesh.tessfaces:
        norm = mesh_to_world.to_3x3() * face.normal
        norm.normalize()
        
        do_face(face.vertices[0], face.vertices[1], face.vertices[2], norm)
        if len(face.vertices)==4:
            do_face(face.vertices[0], face.vertices[2], face.vertices[3], norm)
    
    return (below, above)



def findSplitPoint(normal, mesh, mesh_to_world, task, goal, lam = 1e-3):
    """Finds and returns the distance from the origin required by a plane to acheive a certain goal. There are 3 possible goals, set by the task parameter - 'below', 'above' and 'ratio'. If 'below' then the plane is set so the volume below it is equal to goal, if above then for the volume above it. If set to ratio then goal is the ratio of the total volume that should go below the plane. Uses recursive subdivision, so can be quite slow."""
    assert(task in ['below', 'above', 'ratio'])
    
    # Get the search range...
    low, high = dir_range(normal, mesh, mesh_to_world)
    
    # Do a partial first split - some special casing is required...
    first = True
    half = 0.5 * (low + high)
    below, above = volumeSplitPlane(normal, half, mesh, mesh_to_world)
    
    vol = below + above
    high_below = vol
    high_above = 0.0
    low_below = 0.0
    low_above = vol
    
    # Keep splitting until we are close enough...
    while (high - low) > lam:
        if not first:
            # Choose a split point - use linear interpolation to try and reduce the iteration count...
            if task=='below':
                t = (goal - low_below) / max(high_below - low_below, 1e-3)
            elif task=='above':
                t = (goal - low_above) / max(high_above - low_above, 1e-3)
            else: # task=='ratio'
                lr = low_below / (low_below + low_above)
                hr = high_below / (high_below + high_above)
                t = (goal - lr) / max(hr - lr, 1e-5)
            
            if t<0.1: t = 0.1
            elif t>0.9: t = 0.9
            
            half = (1.0-t) * low + t * high
        
            # Do the split...
            below, above = volumeSplitPlane(normal, half, mesh, mesh_to_world)
        else:
            first = False
        
        # Update the ranges...
        if task=='below':
            goHigh = goal < below
        elif task=='above':
            goHigh = goal > above
        else: # task=='ratio'
            goHigh = (goal * (below + above)) < below
        
        if goHigh:
            high = half
            high_below = below
            high_above = above
        else:
            low = half
            low_below = below
            low_above = above
    
    print('Below =', low_below, high_below)
    print('Above =', low_above, high_above)
    print('Selected offset from origin =', 0.5*(low + high))
    return 0.5 * (low + high)



class VolumeOperator(bpy.types.Operator):
    bl_idname = 'object.volume'
    bl_label = 'Volume'
    
    @classmethod
    def poll(cls, context):
        return (context.active_object is not None) and context.active_object.type=='MESH'

    def execute(self, context):
        mesh = context.active_object
        final_mesh = mesh.to_mesh(bpy.context.scene, True, 'RENDER')
        below, above = volumeSplitPlane(Vector((0.0, 0.0, 1.0)), 0.0, final_mesh, mesh.matrix_world)
        
        vol = below + above
        #self.report({'INFO'}, 'Volume = %.3f BU^3' % vol)
        mesh['volume'] = vol
        
        return {'FINISHED'}



class VolumeRecordOperator(bpy.types.Operator):
    bl_idname = 'object.volume_record'
    bl_label = 'Record Volume'
    
    @classmethod
    def poll(cls, context):
        return (context.active_object is not None) and context.active_object.type=='MESH'
    
    def execute(self, context):
        # Get the current frame range...
        start = context.scene.frame_start
        end = context.scene.frame_end
        
        # Iterate it, and for each frame calculate the volume and record it into the volume property...
        orig_frame = context.scene.frame_current
        mesh = context.active_object
        
        for frame in range(start, end+1):
            context.scene.frame_set(frame)
            print('Frame %i' % frame)
            
            final_mesh = mesh.to_mesh(bpy.context.scene, True, 'RENDER')
            below, above = volumeSplitPlane(Vector((0.0, 0.0, 1.0)), 0.0, final_mesh, mesh.matrix_world)
        
            vol = below + above
            mesh['volume'] = vol
            mesh.keyframe_insert('["volume"]')

        context.scene.frame_set(orig_frame)
        
        return {'FINISHED'}



class SplitVolumeOperator(bpy.types.Operator):
    bl_idname = 'object.split_volume'
    bl_label = 'Split Volume'
    
    @classmethod
    def poll(cls, context):
        ok = context.active_object is not None
        ok = ok and len(context.selected_objects)==2
        for object in context.selected_objects:
            if object!=context.active_object:
                ok = ok and object.type=='MESH'
        return ok
    
    def execute(self, context):
        # Get the plane and the mesh, including the planes normal...
        plane = context.active_object
        normal = plane.matrix_world.to_3x3() * Vector((0.0, 0.0, 1.0))
        normal.normalize()
        
        for object in context.selected_objects:
            if object!=plane: mesh = object
        
        # Determine which mode we need to use, and get the relevent value...
        if 'ratio' in plane.keys():
            task = 'ratio'
            goal = plane['ratio']
        elif 'above' in plane.keys():
            task = 'above'
            goal = plane['above']
        elif 'below' in plane.keys():
            task = 'below'
            goal = plane['below']
        else:
            self.report({'ERROR'}, 'Active object needs a custom property')
            return {'CANCELLED'}
        
        # Find the relevent offset of the plane from the origin...
        mesh_to_world = mesh.matrix_world
        m = mesh.to_mesh(bpy.context.scene, True, 'RENDER')
        final_dist = findSplitPoint(normal, m, mesh_to_world, task, goal)
        
        # Move it the required distance...
        start_dist = normal.dot(plane.matrix_world * Vector((0.0, 0.0, 0.0, 1.0)))
        offset = normal * (final_dist - start_dist)
        plane.matrix_world = Matrix.Translation(offset) * plane.matrix_world
        
        return {'FINISHED'}



class SplitVolumeRecordOperator(bpy.types.Operator):
    bl_idname = 'object.split_volume_record'
    bl_label = 'Record Split Volume'
    
    @classmethod
    def poll(cls, context):
        ok = context.active_object is not None
        ok = ok and len(context.selected_objects)==2
        for object in context.selected_objects:
            if object!=context.active_object:
                ok = ok and object.type=='MESH'
        return ok
    
    def execute(self, context):
        # Get the plane and the mesh...
        plane = context.active_object
        
        for object in context.selected_objects:
            if object!=plane: mesh = object
        
        # Determine which mode we need to use, and get the relevent value...
        if 'ratio' in plane.keys():
            task = 'ratio'
        elif 'above' in plane.keys():
            task = 'above'
        elif 'below' in plane.keys():
            task = 'below'
        else:
            self.report({'ERROR'}, 'Active object needs a custom property')
            return {'CANCELLED'}
        
        # Get the current frame range...
        start = context.scene.frame_start
        end = context.scene.frame_end
        orig_frame = context.scene.frame_current
        
        # Loop and process for each frame...
        for frame in range(start, end+1):
            context.scene.frame_set(frame)
            print('Frame %i' % frame)
            
            # Get the planes normal...
            normal = plane.matrix_world.to_3x3() * Vector((0.0, 0.0, 1.0))
            normal.normalize()
            
            # Find the relevent offset of the plane from the origin...
            mesh_to_world = mesh.matrix_world
            m = mesh.to_mesh(bpy.context.scene, True, 'RENDER')
            goal = plane[task]
            final_dist = findSplitPoint(normal, m, mesh_to_world, task, goal)
        
            # Move it the required distance...
            start_dist = normal.dot(plane.matrix_world * Vector((0.0, 0.0, 0.0, 1.0)))
            offset = normal * (final_dist - start_dist)
            plane.matrix_world = Matrix.Translation(offset) * plane.matrix_world
            
            # Add the keyframe...
            plane.keyframe_insert('location')

        context.scene.frame_set(orig_frame)
        return {'FINISHED'}



def register():
    bpy.utils.register_class(VolumeOperator)
    bpy.utils.register_class(VolumeRecordOperator)
    bpy.utils.register_class(SplitVolumeOperator)
    bpy.utils.register_class(SplitVolumeRecordOperator)


def unregister():
    bpy.utils.unregister_class(VolumeOperator)
    bpy.utils.unregister_class(VolumeRecordOperator)
    bpy.utils.unregister_class(SplitVolumeOperator)
    bpy.utils.unregister_class(SplitVolumeRecordOperator)


register()
#if __name__ == "__main__":
#    register()
