#!/usr/local/intel/2020/intelpython3/bin/python

from __future__ import division
from ase import Atom, Atoms
from ase.io import read,write
import numpy as np
import sys
import os
from ase.visualize import view
from ase.constraints import *

#Tej Choksi, 10/15/2017 - edited 07/24/2020

#USAGE:
#for top position
#python addLayer.py POSCAR_adsorbate POSCAR_surface n_ad top n_sufatom

#addLayer.py is the file, corresponds to sys.argv[0]
#POSCAR_adsorbate is POSCAR file of adsorbate, corresponds to sys.argv[1]
#POSCAR_surface is POSCAR file of surface --> sys.argv[2]
#n_ad is adsorbate atom to be put on top of surface --> sys.argv[3]
#top is the top position
#n_sufatom is surface atom where adsorbate atom is placed on top 

#for bridge position
#python ... bridge n_sufatom1 n_sufatom2

#for hcp or fcc 
#python ... hcp/fcc n_sufatom1 n_sufatom2 n_sufatom3

#for user input
#python ... user height z_pos x_pos y_pos 


# python addLayer.py 3 top 37 
# python addLayer.py 0 user 1.3 40 pos_x pos_y
# 3 - 3rd atom from old structure  is the reference
# top - site where it is to be placed, can also be bridge, fcc, hcp
# 37  etc one, two or three positions that the reference atom sits on

#Read the oxide layer to be added
slab_ads = read(sys.argv[1])

#Read positions from the adsorbate geometry
slab_at = slab_ads.get_number_of_atoms()
pos_slab = slab_ads.get_positions()

#Get Reference coordinate from the old adsorbate
print(sys.argv[3])
r_C1 = pos_slab[int(sys.argv[3])]

#Read the new metal slab being considered
#slab = read('POSCAR')
slab = read(sys.argv[2])

#Select the adsorption site
a = sys.argv[4]

pos = slab.get_positions()

#Get position to which the adsorbate is being shifted

c1 = np.zeros(3)

if a == 'hcp' or a == 'fcc':
	print(a)

	c1[2] = 1.8+slab[int(sys.argv[5])].z

	r_1 = pos[int(sys.argv[5]),:]
	r_2 = pos[int(sys.argv[6]),:]
	r_3 = pos[int(sys.argv[7]),:]

	#This position is an fcc site
	for i in range(len(c1)-1):
		c1[i] = (1/3)*(r_1[i] + r_2[i]+ r_3[i])

elif a == 'bridge':
	print(a)

	c1[2] = 2.277+slab[int(sys.argv[5])].z

	r_1 = pos[int(sys.argv[5]),:]
	r_2 = pos[int(sys.argv[6]),:]

	for i in range(len(c1)-1):
		c1[i] = 0.5*(r_1[i] + r_2[i])

elif a == 'top':
	print(a)

	c1[2] = 2 + slab[int(sys.argv[5])].z
	print(slab[int(sys.argv[5])].z)
	r_1 = pos[int(sys.argv[5]),:]

	for i in range(len(c1)-1):
		c1[i] = r_1[i]

elif a == 'user':
	print(a)
	c1[2] = float(sys.argv[5])+slab[int(sys.argv[6])].z
	c1[0] =  float(sys.argv[7]) 
	c1[1] =  float(sys.argv[8])
else:
	print('Enter either hcp, bridge, fcc, or top')

#Get unit vector connecting old and new adsorbate positions

dc1_C1 = np.zeros(3)

for i in range(len(r_C1)):
	dc1_C1[i] = c1[i] - r_C1[i]

magc1_C1 = np.linalg.norm(dc1_C1)

unitc1_C1 = np.divide(dc1_C1, magc1_C1)

#Shift adsorbate to new position
pos_slab.resize((slab_at,3))

for i in range(pos_slab.shape[0]):
	for j in range(pos_slab.shape[1]):
		pos_slab[i,j] = pos_slab[i,j] + magc1_C1*unitc1_C1[j]

slab_ads.set_positions(pos_slab)

for atom in slab_ads:
	slab.append(Atom(atom.symbol, atom.position))

slab.center(vacuum = 8, axis = 2)
view(slab)
slab.write('POSCAR_new')
os.system("cp POSCAR_new ../POSCAR")
