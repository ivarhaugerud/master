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

C_max    = np.zeros((datafiles, len(C[0,0,0,:])))

for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		#C[:,:,i,f] /= np.sum(np.sum(C[:,:,i,f]))
		C_max[i, f] = np.max(np.max(C[:,:,i,f]))

plt.plot(t, C_max[:, 0], label="No cutoff")
plt.plot(t, C_max[:, 1], label="0.5 cutoff")
plt.plot(t, C_max[:, 2], label="0.8 cutoff")
plt.plot(t, C_max[:, 3], label="0.9 cutoff")
plt.plot(t, C_max[:, 4], label="Single maximum")
plt.plot(t, C_max[:, 5], label="Max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Variance of concentration", fontsize=14)
plt.legend(loc="best", fontsize=14)
#plt.savefig("../powerpoint/figures/variance_high_D.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/variance_high_D.pdf", "../powerpoint/figures/variance_high_D.pdf"))
plt.show()