import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata


filenames = ["data_square/U_Re0.0_b0.2.h5", "data_square/U_Re0.0_b0.6.h5", "data_square/U_Re0.0_b1.4.h5"]#, "data_square/U_Re99.1296244009446_b1.4.h5", "data_square/U_Re0.0_b1.4.h5"]

for k in range(len(filenames)):
	filename = filenames[k]
	print(filename)
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

	ux = ux.roll(shift=int(Nx/2), axis=0)
	uy = uy.roll(shift=int(Ny/2), axis=1)

	b = 0.5*max(x)

	for i in range(len(x)):
		for j in range(len(y)):
			if b/2 < x[i] < 3*b/2 and abs(y[j]) > 0.5*(2-b):
				ux[j,i] = 0
				uy[j,i] = 0

	plt.figure(k)
	speed = np.sqrt(np.square(ux) + np.square(uy)) # / speed.max()
	plt.streamplot(X, Y, ux, uy, density=1.5, color='k')
	CS = plt.contourf(X, Y, speed, 100)
	cbar = plt.colorbar(CS)
	#plt.savefig("figures/test_"+str(i))
plt.show()



