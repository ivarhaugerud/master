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

datafiles = 1000
Nx = 160
Ny = 64
C = np.zeros((Nx, Ny, datafiles))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)


max_C = np.max(np.max(np.max(C)))

fig,ax = plt.subplots()

def animate(i):
	data = np.loadtxt("../data/C_"+str(i+100)+"_back.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	ax.clear()
	ax.contourf(C[:,:,i])#, levels=np.linspace(-100, np.max(np.max(C[:,:,i])), 50))
	ax.axis("equal")


interval = 5 #in ms     
ani = animation.FuncAnimation(fig,animate,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
