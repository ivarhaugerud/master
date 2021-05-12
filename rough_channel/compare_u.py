import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata
import matplotlib.ticker as tick 

name = "results_oscwavychannel/Lx2.0_tau1.0_eps0.2_nu0.8_D0.1_fzero0.0_fone3.0_res50_dt0.02/u.h5"

f = h5py.File(name, 'r')
#index of all lattice-points, 0 could be any timestep as this is indep of time
geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))
timesteps = 100#len(list(f["VisualisationVector"]))
x, y = geometry[:,0], geometry[:,1]
speed = np.zeros(timesteps)

for i in range(timesteps):
	u = np.array(list(f["VisualisationVector"][str(i)]))
	speed[i] = np.mean(np.sqrt(np.square(u[:,0]) + np.square(u[:,1])))
plt.plot(speed)
plt.show()



epsilon = 0.2
kappa = np.pi
eta = np.linspace(0, 2, 100)

for i in range(timesteps):
	plt.clf()
	u = np.array(list(f["VisualisationVector"][str(i)]))
	Nx = int(2.5*1e3)
	Ny = int(2.5*1e3)
	x = np.linspace(min(geometry[:,0]), max(geometry[:,0]), Nx)
	y = np.linspace(min(geometry[:,1]), max(geometry[:,1]), Ny)
	X, Y = np.meshgrid(x,y)
	
	# Interpolate uneven grid onto an even grid
	ux = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='nearest')
	uy = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='nearest')
	speed = np.sqrt(np.square(ux) + np.square(uy))

	ux[np.where(np.abs(uy)<1e-5)] = 0
	uy[np.where(np.abs(uy)<1e-5)] = 0

	plt.streamplot(X, Y, ux, uy, density=1.0, color='k')
	CS = plt.contourf(X, Y, speed, levels=np.linspace(0, 0.4, 15))
	cbar = plt.colorbar(CS)
	cbar.set_label(r"Velocity $[U]$", fontsize=12)
	cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter("%.1f"))
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	#plt.fill_between(x[0,:], (-1-epsilon)*np.ones(len(x[0,:]))-0.05, y[0,:], color="k")
	#plt.fill_between(x[0,:], ( 1+epsilon)*np.ones(len(x[0,:]))+0.05, -y[0,:], color="k")
	plt.axis([min(eta), max(eta), -1-epsilon-0.05, 1+epsilon+0.05])
	#plt.axis("equal")
	plt.draw()
	plt.pause(0.05)
plt.show()
