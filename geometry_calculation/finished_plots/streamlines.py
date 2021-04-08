import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata
import matplotlib.ticker as tick
import matplotlib
Map = matplotlib.cm.get_cmap('Spectral_r')

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

base = "../data_test/plot_u/"
filenames = [#"Lx6.28_tau3.0_eps0.3_nu1.2_D1.0_fzero0.0_fone12.0_res100_dt0.006/u.h5",
			 "Lx6.28_tau3.0_eps0.7_nu1.2_D1.0_fzero0.0_fone12.0_res100_dt0.006/u.h5",
			 "Lx6.28_tau6.0_eps0.3_nu1.2_D1.0_fzero0.0_fone12.0_res100_dt0.012/u.h5"]

kappa = 1
pi = np.pi
epsilon = np.array([0.7, 0.3, 0.7, 0.3])
timesteps = int(3/0.006)
frames = 18
skip = int(timesteps/frames)
skip = 1
periods_of_flow = 8/3

for k in range(len(filenames)):
	fig = plt.figure(k)
	ax = fig.add_subplot(111)
	ax.set_aspect('equal')
	geometry = np.array(list(h5py.File(base+filenames[k], 'r')["Mesh"]["0"]["mesh"]["geometry"]))
	x_axis = geometry[:,0]
	y_axis = geometry[:,1]

	Nx = int(1e3)
	Ny = int(1e3)
	x = np.linspace(min(x_axis), max(x_axis), Nx)
	y = np.linspace(min(y_axis), max(y_axis), Ny)
	stream_points    = np.array(list(zip( np.zeros(20), np.linspace(-1-epsilon[k], 1+epsilon[k], 10))))
	stream_points2   = np.array(list(zip( max(x)*np.ones(20), np.linspace(-1-epsilon[k], 1+epsilon[k], 10))))
	X, Y = np.meshgrid(x,y)

	for j in range(int(17)):
		plt.clf()
		#u = np.array(list(h5py.File(base+filenames[k], 'r')["VisualisationVector"][str(int(-250+periods_of_flow*timesteps)-j*skip-1)])) 
		u = np.array(list(h5py.File(base+filenames[k], 'r')["VisualisationVector"][str(int(900)-j*skip-1)])) 

		# Interpolate using three different methods and plot
		ux = griddata((x_axis, y_axis), u[:,0], (X, Y), method='cubic')
		uy = griddata((x_axis, y_axis), u[:,1], (X, Y), method='cubic')

		#ux = np.roll(ux, shift=int(Nx/2), axis=1)
		#uy = np.roll(uy, shift=int(Nx/2), axis=1)

		speed = np.sqrt(np.square(ux) + np.square(uy))
		plt.gca().set_aspect('equal', adjustable='box')
		#plt.title(str(int(str(int(920)-j*skip-1))))

		#plt.streamplot(X, Y, ux, uy, density=0.6, color='k')
		plt.streamplot(X, Y, ux, uy, color='k', start_points=stream_points,  density=5)
		plt.streamplot(X, Y, ux, uy, color='k', start_points=stream_points2, density=5)
		CS = plt.contourf(X, Y, speed, 10, cmap="Spectral", map=Map)
		cbar = plt.colorbar(CS)#, fraction=0.02725, pad=0.02)
		cbar.ax.set_ylabel('Velocity $[U]$', fontsize=8)
		cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
		plt.tick_params(axis='both', which='major', labelsize=8)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		plt.xlabel(r"Horizontal position $x$   $[a]$  ", fontsize=8)
		plt.ylabel(r"Vertical position   $y$   $[a]$  ", fontsize=8)
		plt.fill_between(x, (-(1+epsilon[k]+0.1))*np.ones(len(x)), -(1+epsilon[k]*np.cos(kappa*x)), color="k")
		plt.fill_between(x, ( (1+epsilon[k]+0.1))*np.ones(len(x)),  (1+epsilon[k]*np.cos(kappa*x)), color="k")
		#plt.fill_between(x+2*np.pi/kappa, (-(1+epsilon[k]))*np.ones(len(x)), y[0,:], color="k")
		#plt.fill_between(x+2*np.pi/kappa, ( (1+epsilon[k]))*np.ones(len(x)), -y[0,:], color="k")
		plt.pause(0.01)
	plt.show()
	"""


	if int(k - 1) % 3 == 0:
		for i in [0, 1]:
			plt.streamplot(X+2*b*i, Y, ux, uy, density=1.3, color='k')
			CS = plt.contourf(X+2*b*i, Y, speed, 10, cmap="Spectral")
			if i == 0:
				cbar = plt.colorbar(CS,  fraction=0.02725, pad=0.02)
				cbar.ax.set_ylabel('Velocity $[U]$', fontsize=16.5)
				cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
			plt.tick_params(axis='both', which='major', labelsize=16.5)
			plt.tick_params(axis='both', which='minor', labelsize=16.5)
			plt.fill_between([0+2*b*i, b/2+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2],     color="k")
			plt.fill_between([0+2*b*i, b/2+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2],     color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [+1+b/2, +1+b/2], [+1-b/2, +1-b/2], color="k")
			plt.fill_between([3*b/2+2*b*i, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1+b/2, -1+b/2], color="k")
			plt.fill_between([0, 2*b+2*b*i], [-1-b/2, -1-b/2], [-1-b_max/2, -1-b_max/2],   color="k")
			plt.fill_between([0, 2*b+2*b*i], [1+b/2, 1+b/2],   [1+b_max/2, 1+b_max/2],     color="k")
			plt.axis([0, i*2*b+2*b, -1-b_max/2, 0])

	if int(k - 2) % 3 == 0:
		plt.streamplot(X, Y, ux, uy, density=2.0, color='k')
		CS = plt.contourf(X, Y, speed, 10, cmap="Spectral")
		cbar = plt.colorbar(CS,  fraction=0.02725, pad=0.02)
		cbar.ax.set_ylabel('Velocity $[U]$', fontsize=16.5)
		cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
		plt.tick_params(axis='both', which='major', labelsize=16.5)
		plt.tick_params(axis='both', which='minor', labelsize=16.5)
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
	plt.savefig(filename, bbox_inches="tight", dpi=750)
	#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
	"""
	#plt.show()



