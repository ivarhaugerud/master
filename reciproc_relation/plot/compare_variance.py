import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import matplotlib
import pandas as pd 

plt.style.use("bmh")
sns.color_palette("hls", 1)

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

Nx = 140
Ny = 64

Dx = 25
Dy = 25

datafiles = 300
t = np.linspace(0, 1, datafiles)
C = np.zeros((Nx, Ny, datafiles, 6))
x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

for i in range(datafiles):
	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_1_0.txt")
	C[:, :, i, 0] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_5.txt")
	C[:, :, i, 1] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_8.txt")
	C[:, :, i, 2] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_9.txt")
	C[:, :, i, 3] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_single_maxima.txt")
	C[:, :, i, 4] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_maxima.txt")
	C[:, :, i, 5] = (np.reshape(data, (Nx, Ny)))

C_m    = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
C_var  = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
C_x    = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
C_xx   = np.zeros((datafiles, Ny, len(C[0,0,0,:])))

C_m_int    = np.zeros((datafiles, len(C[0,0,0,:])))
C_var_int  = np.zeros((datafiles, len(C[0,0,0,:])))
C_x_int    = np.zeros((datafiles, len(C[0,0,0,:])))
C_xx_int   = np.zeros((datafiles, len(C[0,0,0,:])))

import scipy.integrate as sci 
C[:, :, :11+84, 5] = 0 #some noise in the beginning
plt.figure(2)
start = 40
x_axis[-25:] = -np.linspace(25, 1, 25)
print(x_axis)
for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		if i > 190:
			x_axis[-start:] = -np.linspace(start, 1, start)
		else:
			x_axis = np.linspace(0, Nx-1, Nx)
		for j in range(Ny):
			C_m[i,  j, f]	= np.trapz(C[:, j, i, f], x_axis)
			C_x[i,  j, f]   = np.trapz(C[:, j, i, f]*x_axis, x_axis)
			C_xx[i, j, f]   = np.trapz(C[:, j, i, f]*x_axis*x_axis, x_axis)
			#C_var[:,j, f]   = C_xx[:,j,f]/C_m[:,j,f] - (C_x[:,j,f]/C_m[:,j,f])**2

		C_m_int[i, f] = np.trapz(C_m[i,:,f], y_axis)
		C_x_int[i, f] = np.trapz(C_x[i,:,f], y_axis)
		C_xx_int[i, f] = np.trapz(C_xx[i,:,f], y_axis)

		C_var_int[i, f] = C_xx_int[i,f]/C_m_int[i,f] - (C_x_int[i,f]/C_m_int[i,f])**2

plt.plot(t, C_var_int[:, 0], label="No cutoff")
plt.plot(t, C_var_int[:, 1], label="0.5 cutoff")
plt.plot(t, C_var_int[:, 2], label="0.8 cutoff")
plt.plot(t, C_var_int[:, 3], label="0.9 cutoff")
plt.plot(t, C_var_int[:, 4], label="Single maximum")
plt.plot(t, C_var_int[:, 5], label="Max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Variance of concentration", fontsize=14)
plt.legend(loc="best", fontsize=14)
plt.savefig("../powerpoint/figures/variance_high_D.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/variance_high_D.pdf", "../powerpoint/figures/variance_high_D.pdf"))
plt.show()

for i in range(datafiles):
	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_1_0.txt")
	C[:, :, i, 0] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_0_5.txt")
	C[:, :, i, 1] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_0_8.txt")
	C[:, :, i, 2] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_0_9.txt")
	C[:, :, i, 3] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_single_maxima.txt")
	C[:, :, i, 4] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_maxima.txt")
	C[:, :, i, 5] = (np.reshape(data, (Nx, Ny)))

C_m    = np.zeros((datafiles, len(C[0,0,0,:])))
C_var  = np.zeros((datafiles, len(C[0,0,0,:])))
C_x    = np.zeros((datafiles, len(C[0,0,0,:])))
C_xx   = np.zeros((datafiles, len(C[0,0,0,:])))
import scipy.integrate as sci 
C[:, :, :120, 4] = 0 #some noise in the beginning
C[:, :, :120, 5] = 0 #some noise in the beginning

plt.figure(2)
start = 40
x_axis[-25:] = -np.linspace(25, 1, 25)
print(x_axis)
for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		if i > 190:
			x_axis[-start:] = -np.linspace(start, 1, start)
		else:
			x_axis = np.linspace(0, Nx-1, Nx)
		for j in range(Ny):
			C_m[i, f]   += np.trapz(C[:, j, i, f], x_axis)
			C_x[i, f]   += np.trapz(C[:, j, i, f]*x_axis, x_axis)
			C_xx[i, f]  += np.trapz((C[:, j, i, f]*x_axis*x_axis), x_axis)
		C_var[:, f]  = C_xx[:,f]/C_m[:,f] - (C_x[:,f]/C_m[:,f])**2

plt.plot(t, C_var[:, 0], label="no cutoff")
plt.plot(t, C_var[:, 1], label="0.5 cutoff")
plt.plot(t, C_var[:, 2], label="0.8 cutoff")
plt.plot(t, C_var[:, 3], label="0.9 cutoff")
plt.plot(t, C_var[:, 4], label="single maximum")
plt.plot(t, C_var[:, 5], label="max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Variance of concentration", fontsize=14)
plt.legend(loc="best", fontsize=14)
plt.savefig("../powerpoint/figures/variance_low_D.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/variance_low_D.pdf", "../powerpoint/figures/variance_low_D.pdf"))
plt.show()
