import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata

#general variabels for all runs
U  = 1
dx = 1e-4
dt = dx/U
Nr_p_approx = 1000
Re = 0
datafiles = {}

datafiles["Re0"] = ["data_square/U_Re0.0_b0.0.h5",
                    "data_square/U_Re0.0_b0.1.h5",
                    "data_square/U_Re0.0_b0.2.h5",
                    "data_square/U_Re0.0_b0.3.h5",
                    "data_square/U_Re0.0_b0.4.h5",
                    "data_square/U_Re0.0_b0.5.h5",
                    "data_square/U_Re0.0_b0.6.h5",
                    "data_square/U_Re0.0_b0.7.h5",
                    "data_square/U_Re0.0_b0.8.h5",
                    "data_square/U_Re0.0_b0.9.h5",
                    "data_square/U_Re0.0_b1.0.h5",
                    "data_square/U_Re0.0_b1.1.h5",
                    "data_square/U_Re0.0_b1.2.h5",
                    "data_square/U_Re0.0_b1.4.h5",
                    "data_square/U_Re0.0_b1.5.h5"]

for a in range(len(datafiles["Re0"])):
	f = h5py.File(datafiles["Re0"][a], 'r')
	vector = np.array(list(f["VisualisationVector"]["0"]))
	geometry = np.array(list(f["Mesh"]["mesh"]["geometry"]))

	# export grid 
	Nx = int(2.5*1e3)
	Ny = int(2.5*1e3)
	x = np.linspace(min(geometry[:,0]), max(geometry[:,0]), Nx)
	y = np.linspace(min(geometry[:,1]), max(geometry[:,1]), Ny)
	X, Y = np.meshgrid(x,y)
	b = 0.5*max(x)

	# Interpolate uneven grid onto an even grid
	ux = griddata((geometry[:,0], geometry[:,1]), vector[:,0], (X, Y), method='nearest')
	uy = griddata((geometry[:,0], geometry[:,1]), vector[:,1], (X, Y), method='nearest')

	# we are only interested in spreading horisontally, take mean velocity in x-direction
	mean_U = np.mean(vector[:,0])

	# get function to call for random walk simulation
	u_x = sci.RectBivariateSpline(y, x, ux[:,:]/mean_U)
	u_y = sci.RectBivariateSpline(y, x, uy[:,:]/mean_U)

	#visualize field to test if everything worked
	"""
	fig1, ax2 = plt.subplots()
	CS = ax2.contourf(X, Y, ux, 100)
	cbar = fig1.colorbar(CS)
	plt.show()

	fig1, ax2 = plt.subplots()
	CS = ax2.contourf(X, Y, uy, 100)
	cbar = fig1.colorbar(CS)
	plt.show()

	plt.plot(y, u_y(y,    b*0.05, grid=False))
	plt.plot(y, u_y(y,    b, grid=False))
	plt.plot(y, u_y(y,    7*b/2, grid=False))
	plt.show()

	plt.plot(y, u_x(y,    b*0.05, grid=False))
	plt.plot(y, u_x(y,    7*b/2, grid=False))
	plt.plot(y, u_x(y,    b, grid=False))
	plt.show()
	"""
	#run 
	l  = 2*b
	Nt = int(20*l/(U*dt))

	slicing = int(len(geometry)/Nr_p_approx)
	number_of_points = int(len(geometry[::slicing, 0]))

	pos = np.zeros((Nt, number_of_points, 2))
	pos[0, :, 0] = geometry[::slicing, 0]
	pos[0, :, 1] = geometry[::slicing, 1]
	Np = int(len(pos[0,:,0]))

	# if b=0, l=0.05, so just for correct saving
	if abs(b - 0.05) < 1e-4:
		b = 0

	print("Running parameters", "Re=",Re, "b=", b, "Nt=", Nt, "Np=", Np)

	# Running simulation
	for i in range(1, Nt):
	    pos[i, :, 0] = pos[i-1, :, 0] + u_x(pos[i-1,:,1], (pos[i-1,:,0]+l)%l, grid=False)*dt
	    pos[i, :, 1] = pos[i-1, :, 1] + u_y(pos[i-1,:,1], (pos[i-1,:,0]+l)%l, grid=False)*dt

	"""
	print("Finished simulation")
	nr_of_converged = 0
	plt.scatter(pos[0,:,0], pos[0,:,1])
	for j in range(Np):
		print(pos[-1, j, 0] - pos[0,j,0])
		if abs(pos[-1, j, 0] - pos[0,j,0]) < l:
			plt.plot(pos[0,j,0], pos[0,j,1], "ro")
			nr_of_converged += 1
	"""
	np.save("analyzed_data/RZ_b="+str(b)+"_Re="+str(Re), pos)

