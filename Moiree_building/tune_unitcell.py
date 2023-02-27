#!/usr/local/intel/2020/intelpython3/bin/python

from ase import Atoms
from ase.io import read, write
from ase.visualize import view
import sys
import os

#syntax as follows: tune_unitcell.py input_1 input_2
#input_1 is the POSCAR of the targeted unit cell length (e.g POSCAR_slab)
#input_2 is the POSCAR of the unit cell for compressing/stretching (e.g POSCAR_Mol)
#factor = sys.argv[1]

file_name_slab = sys.argv[1]
file_name_Mol = sys.argv[2]
file_name_Mol_new = file_name_Mol + '_scaled'
print(file_name_slab)
print(file_name_Mol)
print(file_name_Mol_new)

slab = read(file_name_slab)
len_x, len_y, len_z, angle_yz, angle_xz, angle_xy = slab.get_cell_lengths_and_angles()
cell = slab.get_cell()
#print(cell)

Mol = read(file_name_Mol)
cell_Mol = Mol.get_cell()
len_x_Mol, len_y_Mol, len_z_Mol, angle_yz_Mol, angle_xz_Mol, angle_xy_Mol = Mol.get_cell_lengths_and_angles()

factor = len_x/len_x_Mol

cell_Mol[0,:] = factor*cell_Mol[0,:]
cell_Mol[1,:] = factor*cell_Mol[1,:]

#Increase the unit cell size slightly to create a distortion

#We scale the atoms to keep the internal symmetries consistent
Mol.set_cell(cell_Mol, scale_atoms = True)      #set the unit cell to the atoms object called slab
view(Mol)               #view the optimized structure
Mol.write(file_name_Mol_new)    #write optimized structure as CONTCAR
