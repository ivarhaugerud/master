import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import matplotlib
import pandas as pd 
import matplotlib.animation as animation
from scipy import ndimage

plt.style.use("bmh")
sns.color_palette("hls", 1)

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

datafiles = 300
var_in = np.zeros(datafiles)
Nx = 140
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)
C = np.zeros((Nx, Ny, datafiles))

u = np.loadtxt("../data/0303heat_heat_u.txt")
u_x = np.reshape(u[0, :], (Nx, Ny))

def animate_forward(i):
	fig.suptitle(str(i))
	data = np.loadtxt("../data/peak_C_"+str(i)+"_maxima.txt")

	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	print(np.sum(np.sum(C[100,:,i])))
	print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), C.flatten()[np.argmax(C[:,:,i])])

	#C[np.where( (C[:,:,i]) < 10) ] = 0
	#C[np.where( abs(u_x[:, :]) < 1e-8)] = -1
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(C[:,:,i])#, levels=np.linspace(-0.00001, 0.0004, 35))
	ax.axis("equal")

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig,ax = plt.subplots()
interval = 1 #in ms     
ani = animation.FuncAnimation(fig,animate_forward,datafiles, interval=interval, blit=False, repeat=False)
#ani.save('im.mp3')
plt.show()

def animate_back(i):
	fig.suptitle(str(i))
	data = np.loadtxt("../data/peak_C_"+str(i)+"_1_0.txt")

	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), np.argmax(C[:,:,i]))

	C[np.where( np.abs(C[:,:,i]) < 0.1) ] = 0
	C[np.where( abs(u_x[:, :]) < 1e-8)] = -1
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(C[:,:,i], levels=np.linspace(-0.01, 5500, 25))
	ax.axis("equal")

fig,ax = plt.subplots()

interval = 1 #in ms     
ani = animation.FuncAnimation(fig,animate_back,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
