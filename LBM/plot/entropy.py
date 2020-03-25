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

u = np.loadtxt("../data/peak_1104_peak_1104_u.txt")
u_x = np.reshape(u[0, :], (Nx, Ny))

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
	

#normalizing
for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		if np.sum(np.sum(C[:,:,i,f])) != 0:
			C[:, :, i, f] = abs(C[:, :, i, f])/np.trapz(np.trapz(abs(C[:,:,i,f]), x_axis, axis=0), y_axis)

S_each_y    = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
mass_each_y = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
mass = np.zeros((datafiles, len(C[0,0,0,:])))
S    = np.zeros((datafiles, len(C[0,0,0,:])))

import scipy.integrate as sci 
C[(np.where(C[:, :, :, :] == 0))] = 1
for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		for j in range(Ny):
			S_each_y[i, j, f] = -np.trapz((abs(C[:, j, i, f]))*np.log((abs(C[:, j, i, f]))), x_axis)
			mass_each_y[i, j, f] = np.trapz(C[:, j, i, f], x_axis)

		S[i, f] = np.trapz(S_each_y[i, :, f], y_axis)

plt.plot(t, S[:, 0], label="No cutoff")
plt.plot(t, S[:, 1], label="0.5 cutoff")
plt.plot(t, S[:, 2], label="0.8 cutoff")
plt.plot(t, S[:, 3], label="0.9 cutoff")
plt.plot(t, S[:, 4], label="Single maximum")
plt.plot(t, S[:, 5], label="Max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Entropy $S = -\int\rho\ln{\rho}$", fontsize=14)
plt.legend(loc="best", fontsize=12, ncol=1)
plt.savefig("../powerpoint/figures/entropy_high_D.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/entropy_high_D.pdf", "../powerpoint/figures/entropy_high_D.pdf"))
plt.show()

S_each_y = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
S    = np.zeros((datafiles, len(C[0,0,0,:])))

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

for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		if np.sum(np.sum(C[:,:,i,f])) != 0:
				C[:, :, i, f] = abs(C[:, :, i, f])/np.trapz(np.trapz(abs(C[:,:,i,f]), x_axis, axis=0), y_axis)

C[(np.where(C[:, :, :, :] == 0))] = 1

for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		for j in range(Ny):
			S_each_y[i, j, f] = -np.trapz((abs(C[:, j, i, f]))*np.log((abs(C[:, j, i, f]))), x_axis)
			mass_each_y[i, j, f] = np.trapz(C[:, j, i, f], x_axis)

		S[i, f] = np.trapz(S_each_y[i, :, f], y_axis)

plt.plot(t, S[:, 0], label="No cutoff")
plt.plot(t, S[:, 1], label="0.5 cutoff")
plt.plot(t, S[:, 2], label="0.8 cutoff")
plt.plot(t, S[:, 3], label="0.9 cutoff")
plt.plot(t, S[:, 4], label="Single maximum")
plt.plot(t, S[:, 5], label="Max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
#plt.yscale("log")
plt.ylabel(r"Entropy $S = -\int\rho\ln{\rho}$", fontsize=14)
plt.legend(loc="best", fontsize=12, ncol=1)
plt.savefig("../powerpoint/figures/entropy_low_D.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/entropy_low_D.pdf", "../powerpoint/figures/entropy_low_D.pdf"))
plt.show()
