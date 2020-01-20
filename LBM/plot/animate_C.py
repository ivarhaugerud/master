import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import matplotlib
import pandas as pd 
import matplotlib.animation as animation

plt.style.use("bmh")
sns.color_palette("hls", 1)

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

datafiles = 10
Nx = 256
Ny = 64
C = np.zeros((Nx, Ny, datafiles))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

for i in range(datafiles):
	data = np.loadtxt("../data/C_"+str(i)+".txt")
	C[:, :, i] = np.reshape(data, (Nx, Ny))
	#C[np.where(C[:,:,i] < 1e-3), i] = 0

fig,ax = plt.subplots()

def animate(i):
	ax.clear()
	ax.contourf(C[:,:,i])
	ax.axis("equal")

interval = 200 #in seconds     
ani = animation.FuncAnimation(fig,animate,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
