import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib
import pandas as pd
import matplotlib.animation as animation
import scipy.integrate as sci

plt.style.use("bmh")
sns.color_palette("hls", 1)

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

datafiles = 100
var_in = np.zeros(datafiles)
Nx = 256
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

C = np.zeros((Nx, Ny, datafiles))
C_r  = np.zeros((datafiles))

for i in range(datafiles):
	data = np.loadtxt("../data/C_"+str(i)+"_front.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	summed = np.sum(C[:,:,i], axis=1)
	C_x = x_axis*summed
	C_xx = x_axis*C_x
	m = np.trapz(summed)
	C_r[i] = np.trapz(C_xx)/m - (np.trapz(C_x)/m)**2


x_inject = np.argmax(np.sum(C[:,:,0], axis=0))
y_inject = np.argmax(np.sum(C[:,:,0], axis=1))

print(x_inject, y_inject)
at_injection_point = np.zeros(datafiles)

plt.plot(C_r, label="x")
plt.legend(loc="best")
plt.show()

for i in range(datafiles):
	data = np.loadtxt("../data/C_"+str(i)+"_back.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	summed = np.sum(C[:,:,i], axis=1)
	C_x = x_axis*summed
	C_xx = x_axis*C_x
	m = np.trapz(summed)
	C_r[i] = np.trapz(C_xx)/m - (np.trapz(C_x)/m)**2

	at_injection_point[i] = C[x_inject, y_inject, i] + C[x_inject+1, y_inject, i] + C[x_inject-1, y_inject, i] + C[x_inject, y_inject+1, i] + C[x_inject, y_inject-1, i]

plt.plot(C_r, label="x")
plt.legend(loc="best")
plt.show()

plt.plot(at_injection_point)
plt.show()
