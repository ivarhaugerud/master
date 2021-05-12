import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata
import matplotlib.ticker as tick
import matplotlib

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
root = "../master_latex/results/"


filenames = ["data_square/U_Re0.0_b0.4.h5", 
			 "data_square/U_Re0.0_b0.8.h5",
			 "data_square/U_Re0.0_b1.6.h5",
			 #"data_square/U_Re3.16289307602155_b0.4.h5",
			 #"data_square/U_Re3.1599334969768544_b0.8.h5",
			 #"data_square/U_Re3.1616049232683583_b1.6.h5",
			 "data_square/U_Re31.624510620418054_b0.4.h5",
			 "data_square/U_Re31.555168251015914_b0.8.h5",
			 "data_square/U_Re30.182298532693334_b1.6.h5"]

rough = ["0.4", "0.8", "1.6", "0.4", "0.8", "1.6", "0.4", "0.8", "1.6"]
Reyno = ["0", "0", "0", "32", "32", "32"] #"3.16", "3.16", "3.16",

b = 0.3
b_max = 1.65
epsilon = 0.0
max_vel = 0

kwargs = {'format': '%.1f'}

Map = matplotlib.cm.get_cmap('Spectral_r')

for k in range(len(filenames)):
	fig = plt.figure(k)
	ax = fig.add_subplot(111)
	ax.set_aspect('equal')
	filename = filenames[k]
	f = h5py.File(filename, 'r')
	vector = np.array(list(f["VisualisationVector"]["0"]))
	geometry = np.array(list(f["Mesh"]["mesh"]["geometry"]))
	x_axis = geometry[:,0]
	y_axis = geometry[:,1]

	Nx = int(1e3)
	Ny = int(1e3)
	x = np.linspace(min(x_axis), max(x_axis), Nx)
	y = np.linspace(min(y_axis), max(y_axis), Ny)
	X, Y = np.meshgrid(x,y)

	# Interpolate using three different methods and plot
	ux = griddata((x_axis, y_axis), vector[:,0], (X, Y), method='cubic')
	uy = griddata((x_axis, y_axis), vector[:,1], (X, Y), method='cubic')
	b = 0.5*max(x)

	for i in range(len(x)):
		for j in range(len(y)):
			if b/2 < x[i] < 3*b/2 and abs(y[j]) > 0.5*(2-b):
				ux[j,i] = 0
				uy[j,i] = 0

	ux = np.roll(ux, shift=int(Nx/2), axis=1)
	uy = np.roll(uy, shift=int(Nx/2), axis=1)

	speed = np.sqrt(np.square(ux) + np.square(uy))
	plt.gca().set_aspect('equal', adjustable='box')

	if k % 3 == 0:
		for i in [0, 1, 2, 3]:
			plt.streamplot(X+2*b*i, Y, ux, uy, density=0.6+0.3, color='k', arrowsize=0.6)
			CS = plt.contourf(X+2*b*i, Y, speed, 10, cmap=Map)
			if i == 0:
				cbar = plt.colorbar(CS, fraction=0.02725, pad=0.02)
				cbar.ax.set_ylabel('Velocity $[U]$', fontsize=8)
				cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
			plt.fill_between([0+2*b*i, b/2+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],     color="k")
			plt.fill_between([0+2*b*i, b/2+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],     color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2], color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2], color="k")
			plt.fill_between([0, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1-b_max/2, -1-b_max/2],   color="k")
			plt.fill_between([0, 2*b+2*b*i], [1+b/2, 1+b/2],   [1+b_max/2, 1+b_max/2],     color="k")
			plt.axis([0, i*2*b+2*b, -1-b_max/2, 0])


	if int(k - 1) % 3 == 0:
		for i in [0, 1]:
			plt.streamplot(X+2*b*i, Y, ux, uy, density=1.3+0.15, color='k', arrowsize=0.6)
			CS = plt.contourf(X+2*b*i, Y, speed, 10, cmap=Map)
			if i == 0:
				cbar = plt.colorbar(CS,  fraction=0.02725, pad=0.02)
				cbar.ax.set_ylabel('Velocity $[U]$', fontsize=8)
				cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
			plt.fill_between([0+2*b*i, b/2+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],     color="k")
			plt.fill_between([0+2*b*i, b/2+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],     color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2], color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2], color="k")
			plt.fill_between([0, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1-b_max/2, -1-b_max/2],   color="k")
			plt.fill_between([0, 2*b+2*b*i], [1+b/2, 1+b/2],   [1+b_max/2, 1+b_max/2],     color="k")
			plt.axis([0, i*2*b+2*b, -1-b_max/2, 0])

	if int(k - 2) % 3 == 0:
		plt.streamplot(X, Y, ux, uy, density=2.0+0.2, color='k', arrowsize=0.6)
		CS = plt.contourf(X, Y, speed, 10, cmap=Map)
		cbar = plt.colorbar(CS,  fraction=0.02725, pad=0.02)
		cbar.ax.set_ylabel('Velocity $[U]$', fontsize=8)
		cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
		plt.fill_between([0, b/2], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],         color="k")
		plt.fill_between([0, b/2], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],         color="k")
		plt.fill_between([3*b/2, 2*b], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],     color="k")
		plt.fill_between([3*b/2, 2*b], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],     color="k")
		plt.fill_between([0, 2*b], [-1-b/2, -1-b/2], [-1-b_max/2, -1-b_max/2], color="k")
		plt.fill_between([0, 2*b], [1+b/2, 1+b/2],   [1+b_max/2, 1+b_max/2],   color="k")
		plt.axis([0, 2*b, -1-b_max/2, 0])

	filename = "figures/streamlines_b=" + rough[k] + "_Re="+ Reyno[k]
	filename = filename.replace(".", "_")
	filename += ".png"
	plt.savefig(root+filename, bbox_inches="tight", dpi=750)
	os.system('pdfcrop %s %s &> /dev/null &'%(root+filename, root+filename))
	plt.tick_params(axis='both', which='major', labelsize=8)
	plt.tick_params(axis='both', which='minor', labelsize=8)
plt.show()



