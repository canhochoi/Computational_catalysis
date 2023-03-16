from ase.io import read, write
from ase.visualize import view
from ase import Atom, Atoms
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.neb import NEBTools
from matplotlib import gridspec
import pandas as pd
from pandas.plotting import table

nimages = 21

images = read('AIDNEB.traj@{}:'.format(-nimages))

nebtools = NEBTools(images)

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = nebtools.get_barrier()

# Get the actual maximum force at this point in the simulation.
max_force = nebtools.get_fmax()
print('Max force: {} ev/A'.format(max_force))

energy_list = np.nan*np.ones((nimages,1))
TS_energy = -1e8
index = 0
for i in range(nimages):
	energy_list[i] = images[i].get_calculator().results['energy']
	if energy_list[i] > TS_energy:
		TS_energy = energy_list[i]
		index = i
#TS_energy = np.amax(energy_list)
#index = np.where(energy_list == TS_energy)
#print(energy_list)
#print(index, TS_energy)
#print(images[0:-1].get_calculator().results['energy'])

fig_neb = nebtools.plot_band()
fig_neb.savefig('spline.png')

#save images
idx_list = [0, index, nimages-1]

#taken from matplotlib in ase
#https://wiki.fysik.dtu.dk/ase/ase/visualize/visualize.html#matplotlib-plotting
for i in idx_list:
	image = images[i]
	#write('image%d.png' % i, image)
	fig, ax = plt.subplots()
	plot_atoms(image, ax)
	fig.savefig('image%d.png' % i)

# create a figure
fig = plt.figure()

# to change size of subplot's
# set height of each subplot as 8
fig.set_figheight(6)

# set width of each subplot as 8
fig.set_figwidth(10)

# create grid for different subplots
spec = gridspec.GridSpec(ncols=3, nrows=2, width_ratios=[1, 1, 1], wspace=0.5, hspace=0.5, height_ratios=[1, 1])
title_name = ['IS','TS','FS']
for idx, val in enumerate(idx_list):
	ax = fig.add_subplot(spec[0,idx])
	image = images[val]
	plot_atoms(image, ax)
	ax.set_title(title_name[idx])
	ax.xaxis.set_visible(False)  # hide the x axis
	ax.yaxis.set_visible(False)  # hide the y axis

#plot spline traj
ax = fig.add_subplot(spec[1,0:2])
nebtools.plot_band(ax = ax)


fig.savefig('combine.png')

#print(Ef)
quantities = ['Ea (eV)', 'Edelta (eV)', 'Max_f (eV/A)', 'IS (eV)', 'TS (eV)', 'FS (eV)']
val_list = [Ef, dE, max_force, energy_list[0][0], TS_energy[0], energy_list[-1][0]]
#val_round = [ round(elem, 2) for elem in val_list]
data_dict = {'Quantity': quantities,
	     'Value': val_list}
df_data = pd.DataFrame(data_dict)
output_file = 'Barrier_info.csv'
df_data.to_csv(output_file)

ax = fig.add_subplot(spec[1,2], frame_on = False)
ax.xaxis.set_visible(False)  # hide the x axis
ax.yaxis.set_visible(False)  # hide the y axis
#tbl = table(ax, df_data, rowLabels = ['']*df_data.shape[0],loc = 'center')

#ref link
#https://stackoverflow.com/questions/19726663/how-to-save-the-pandas-dataframe-series-data-as-a-figure

tbl = ax.table(cellText = df_data.values, colLabels = df_data.columns, loc = 'center')
tbl.auto_set_column_width(col=list(range(len(df_data['Value'])))) 
tbl.auto_set_font_size(False)
tbl.set_fontsize(12)
fig.savefig('combine.png')
print(df_data.values[0])

#display(fig)
#fig.savefig('combine.png')
#with pd.ExcelWriter(output_file) as writer:
#	df_data.to_excel(writer) 
