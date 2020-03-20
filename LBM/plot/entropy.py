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
	C[:, :, i, 0] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*1e2

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_5.txt")
	C[:, :, i, 1] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*1e2

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_8.txt")
	C[:, :, i, 2] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*1e2

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_9.txt")
	C[:, :, i, 3] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*1e2

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_single_maxima.txt")
	C[:, :, i, 4] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*1e2

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_maxima.txt")
	C[:, :, i, 5] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*1e2

for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		C[:, :, i, f] /= np.sum(np.sum(C[:,:,i,f]))
	print(np.sum(np.sum(np.sum(C[:,:,:,f], axis=0),axis=0), axis=0))

S_each_y    = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
mass_each_y = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
mass = np.zeros((datafiles, len(C[0,0,0,:])))
S    = np.zeros((datafiles, len(C[0,0,0,:])))

import scipy.integrate as sci 
plt.figure(2)

for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		for j in range(Ny):
			#C[np.where( abs((u_x[:, :])) < 1e-8), i, f] = 0

			S_each_y[i, j, f] = np.trapz((C[:, j, i, f])*np.log((C[:, j, i, f])), x_axis)
			mass_each_y[i, j, f] = np.trapz(C[:, j, i, f], x_axis)

		S[i, f] = np.trapz(S_each_y[i, :, f], y_axis)
		mass[i, f] = np.trapz(mass_each_y[i, :, f], y_axis)
	S[:, f] - S[0,f]+1

plt.plot(mass)
plt.show()

plt.plot(t, S[:, 0], label="No cutoff")
plt.plot(t, S[:, 1], label="0.5 cutoff")
plt.plot(t, S[:, 2], label="0.8 cutoff")
plt.plot(t, S[:, 3], label="0.9 cutoff")
plt.plot(t, S[:, 4], label="Single maximum")
plt.plot(t, S[:, 5], label="Max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
#plt.yscale("log")
plt.ylabel(r"Entropy $S = \rho\ln{\rho}$", fontsize=14)
plt.legend(loc="best", fontsize=12, ncol=1)
plt.savefig("../powerpoint/figures/entropy_high_D.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/entropy_high_D.pdf", "../powerpoint/figures/entropy_high_D.pdf"))
plt.show()

S_each_y = np.zeros((datafiles, Ny, len(C[0,0,0,:])))
S    = np.zeros((datafiles, len(C[0,0,0,:])))

for i in range(datafiles):
	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_1_0.txt")
	C[:, :, i, 0] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*5*1e2

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_0_5.txt")
	C[:, :, i, 1] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*5*1e2

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_0_8.txt")
	C[:, :, i, 2] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*5*1e2

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_0_9.txt")
	C[:, :, i, 3] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*5*1e2

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_single_maxima.txt")
	C[:, :, i, 4] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*5*1e2

	data = np.loadtxt("../data/peak_1104_C_"+str(i)+"_maxima.txt")
	C[:, :, i, 5] = (np.reshape(data, (Nx, Ny))) + np.ones((Nx, Ny))*5*1e2

for i in range(datafiles):
	C[:, :, i, 0] /= np.sum(np.sum(C[:,:,i,0]))
	C[:, :, i, 1] /= np.sum(np.sum(C[:,:,i,1]))
	C[:, :, i, 2] /= np.sum(np.sum(C[:,:,i,2]))
	C[:, :, i, 3] /= np.sum(np.sum(C[:,:,i,3]))
	C[:, :, i, 4] /= np.sum(np.sum(C[:,:,i,4]))
	C[:, :, i, 5] /= np.sum(np.sum(C[:,:,i,5]))

for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		for j in range(Ny):
			S_each_y[i, j, f] = np.trapz((C[:, j, i, f])*np.log((C[:, j, i, f])), x_axis)
		S[i, f] = np.trapz(S_each_y[i, :, f], y_axis)
	S[:, f] - S[0,f]+1

plt.plot(t, S[:, 0], label="No cutoff")
plt.plot(t, S[:, 1], label="0.5 cutoff")
plt.plot(t, S[:, 2], label="0.8 cutoff")
plt.plot(t, S[:, 3], label="0.9 cutoff")
plt.plot(t, S[:, 4], label="Single maximum")
plt.plot(t, S[:, 5], label="Max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
#plt.yscale("log")
plt.ylabel(r"Entropy $S = \rho\ln{\rho}$", fontsize=14)
plt.legend(loc="best", fontsize=12, ncol=1)
plt.savefig("../powerpoint/figures/entropy_low_D.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/entropy_low_D.pdf", "../powerpoint/figures/entropy_low_D.pdf"))
plt.show()

