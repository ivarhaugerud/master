import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata
import matplotlib.ticker as tick
import matplotlib

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

filenames = ["data_square/U_Re0.0_b0.4.h5", 
			 "data_square/U_Re0.0_b0.8.h5",
			 "data_square/U_Re0.0_b1.6.h5",
			 "data_square/U_Re31.624510620418054_b0.4.h5",
			 "data_square/U_Re31.555168251015914_b0.8.h5",
			 "data_square/U_Re30.182298532693334_b1.6.h5"]

rough = ["0.4", "0.8", "1.6", "0.4", "0.8", "1.6"]
Reyno = ["0", "0", "0", "31", "31", "31"]

b = 0.3
b_max = 1.65
epsilon = 0.0
max_vel = 0

kwargs = {'format': '%.1f'}

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
			plt.streamplot(X+2*b*i, Y, ux, uy, density=1.0, color='k')
			CS = plt.contourf(X+2*b*i, Y, speed, 100, cmap="Spectral")
			if i == 0:
				cbar = plt.colorbar(CS, fraction=0.02725, pad=0.02)
				cbar.ax.set_ylabel('Velocity $[U]$', fontsize=12)
				cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
			plt.fill_between([0+2*b*i, b/2+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],     color="k")
			plt.fill_between([0+2*b*i, b/2+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],     color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2], color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2], color="k")
			plt.fill_between([0, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1-b_max/2, -1-b_max/2],   color="k")
			plt.fill_between([0, 2*b+2*b*i], [1+b/2, 1+b/2],   [1+b_max/2, 1+b_max/2],     color="k")
			#plt.axis([0, i*2*b+2*b, -1-b_max/2, 0])


	if int(k - 1) % 3 == 0:
		for i in [0, 1]:
			plt.streamplot(X+2*b*i, Y, ux, uy, density=1.3, color='k')
			CS = plt.contourf(X+2*b*i, Y, speed, 100, cmap="Spectral")
			if i == 0:
				cbar = plt.colorbar(CS,  fraction=0.02725, pad=0.02)
				cbar.ax.set_ylabel('Velocity $[U]$', fontsize=12)
				cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
			plt.fill_between([0+2*b*i, b/2+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],     color="k")
			plt.fill_between([0+2*b*i, b/2+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],     color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2], color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2], color="k")
			plt.fill_between([0, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1-b_max/2, -1-b_max/2],   color="k")
			plt.fill_between([0, 2*b+2*b*i], [1+b/2, 1+b/2],   [1+b_max/2, 1+b_max/2],     color="k")
			#plt.axis([0, i*2*b+2*b, -1-b_max/2, 0])

	if int(k - 2) % 3 == 0:
		plt.streamplot(X, Y, ux, uy, density=2.0, color='k')
		CS = plt.contourf(X, Y, speed, 100, cmap="Spectral")
		cbar = plt.colorbar(CS,  fraction=0.02725, pad=0.02)
		cbar.ax.set_ylabel('Velocity $[U]$', fontsize=12)
		cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
		plt.fill_between([0, b/2], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],         color="k")
		plt.fill_between([0, b/2], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],         color="k")
		plt.fill_between([3*b/2, 2*b], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],     color="k")
		plt.fill_between([3*b/2, 2*b], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],     color="k")
		plt.fill_between([0, 2*b], [-1-b/2, -1-b/2], [-1-b_max/2, -1-b_max/2], color="k")
		plt.fill_between([0, 2*b], [1+b/2, 1+b/2],   [1+b_max/2, 1+b_max/2],   color="k")
		#plt.axis([0, 2*b, -1-b_max/2, 0])

	filename = "figures/streamlines_b=" + rough[k] + "_Re="+ Reyno[k]
	filename = filename.replace(".", "_")
	filename += ".pdf"
	plt.savefig(filename, bbox_inches="tight")
	os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
	plt.show()



