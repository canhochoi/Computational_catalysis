from ase import Atoms
from ase.visualize import view
import numpy as np
from ase.io import read, write
from ase.build import fcc111, rotate
from ase.build import root_surface
from ase.build import root_surface_analysis
import math

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def point_in_cell_2d(point, cell, eps=1e-8):
    """This function takes a 2D slice of the cell in the XY plane and calculates
    if a point should lie in it.  This is used as a more accurate method of
    ensuring we find all of the correct cell repetitions in the root surface
    code.  The Z axis is totally ignored but for most uses this should be fine.
    """
    # Define area of a triangle
    def tri_area(t1, t2, t3):
        t1x, t1y = t1[0:2]
        t2x, t2y = t2[0:2]
        t3x, t3y = t3[0:2]
        return abs(t1x * (t2y - t3y) + t2x * (t3y - t1y) + t3x * (t1y - t2y)) / 2

    # c0, c1, c2, c3 define a parallelogram
    c0 = (0, 0)
    c1 = cell[0, 0:2]
    c2 = cell[1, 0:2]
    c3 = c1 + c2

    # Get area of parallelogram
    cA = tri_area(c0, c1, c2) + tri_area(c1, c2, c3)

    # Get area of triangles formed from adjacent vertices of parallelogram and
    # point in question.
    pA = tri_area(point, c0, c1) + tri_area(point, c1, c2) + tri_area(point, c2, c3) + tri_area(point, c3, c0)

    # If combined area of triangles from point is larger than area of
    # parallelogram, point is not inside parallelogram.
    return pA <= cA + eps

#for limiting in 2D plane 
eps = 1e-8

#values to take root
#need to input  
root = 13
#root value
root_value = np.sqrt(root)

#index of atoms to get vectors and angle rotation (from ntimes repeating structure)
#need to input
index_0 = 1
index_1 = 64
#print(pos_mol[index_1])

#import structure
slab = read('POSCAR')
ntimes = 5 #number of times to repeat
slab_repeat = slab.repeat((ntimes,ntimes,1))
#view(slab_repeat)
cell_vector = slab.cell
vector0_x = cell_vector[0,:]
vector0_y = cell_vector[1,:]
z = cell_vector[2,2]
#print(cell_vector)
#print(vector0_y)
pos_mol = slab_repeat.get_positions()
#print(pos_mol)

#index of atoms to get vectors and angle rotation
#index_0 = 0
#index_1 = 42
#print(pos_mol[index_1])

#rotation angle
a_i = cell_vector[0,:]   
a_f = pos_mol[index_1] - pos_mol[index_0]
rotation_angle = angle_between(a_i,a_f)
rotation_angle_deg = math.degrees(rotation_angle)
#print(rotation_angle)


#rotation of x-y axis
# minus sign to rotate counter-clockwise
rotation_angle = -rotation_angle 
rotation_dir = [[np.cos(rotation_angle), -np.sin(rotation_angle)],
                [np.sin(rotation_angle), np.cos(rotation_angle)]]


#new unit cell after rotation
cell_vectorxy = cell_vector[0:2,0:2]
#print(cell_vectorxy)
#vector length
vector_length = np.sqrt(np.sum(np.square(cell_vectorxy[1])))
#vector_length = round(vector_length, 2) #round to 2 decimals
#print(vector_length)

cell_length = vector_length*root_value

#normalize
cell_vectorxy = cell_vectorxy/np.linalg.norm(cell_vectorxy[0])

cell_new = np.array([np.dot(x,rotation_dir)*cell_length for x in cell_vectorxy])

#STEP0 - set up new unit cell for root structures
cell_new3D = np.append(cell_new, [[0,0]], 0)
cell_new3D = np.append(cell_new3D, [[0],[0],[z]], 1)
print(cell_new3D)
print(vector0_x)
print(vector0_y)

#STEP1 - duplicate atoms
atoms = slab.copy()   ##IMPORTANT to avoid changing of vector0_x into vector of new cell unit
atoms *= (2*root, 2*root, 1)
view(atoms)
#translate backward so that new unit cell covers all atoms
shift_vector = [-2.51*vector_length,0, 0]
atoms.translate(shift_vector)

#STEP2 - set new unit cell
atoms.cell[0:2, 0:2] = cell_new
view(atoms)
#STEP3 - remove extra atoms
del atoms[[atom.index for atom in atoms if not point_in_cell_2d(atom.position, atoms.cell, eps=eps)]]

atoms.cell = cell_new3D
view(atoms)
print(vector0_x)
#STEP4 - rotate both atoms and new unit cell to original orientation of slab (111)
#rotate(atoms,cell_new3D[0,:],[1,0,0],cell_new3D[1,:],[1/2,np.sqrt(3)/2,0])
#rotate(atoms,cell_new3D[0,:],vector0_x,cell_new3D[1,:],vector0_y)
atoms.rotate(-13.90,'z',rotate_cell = True )
atoms = atoms.repeat((1,1,1))
#print(atoms.get_tags())
atoms.center(vacuum = 8, axis = 2)
view(atoms)

atoms.write('POSCAR_root')

#print(cell_new3D)
#print(vector0_x)
