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

datafiles = 100
var_in = np.zeros(datafiles)
Nx = 256
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)
C = np.zeros((Nx, Ny, datafiles))

def animate_forward(i):
	fig.suptitle(str(i))
	data = np.loadtxt("../data/C_"+str(i)+"_front.txt")
	u = np.loadtxt("../data/final_vel.txt")
	u_x = np.reshape(u[0, :], (Nx, Ny))

	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	C[np.where( np.abs(C[:,:,i]) < 0.01) ] = 0
	C[np.where( abs(u_x[:, :]) < 1e-8)] = -100
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(C[:,:,i], levels=np.linspace(-0.1, 0.5, 25))
	ax.axis("equal")
	print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), np.argmax(C[:,:,i]))

fig,ax = plt.subplots()
interval = 5 #in ms     
ani = animation.FuncAnimation(fig,animate_forward,datafiles, interval=interval, blit=False, repeat=False)
plt.show()


def animate_back(i):
	fig.suptitle(str(i))
	data = np.loadtxt("../data/C_"+str(i)+"_back.txt")
	u = np.loadtxt("../data/final_vel.txt")
	u_x = np.reshape(u[0, :], (Nx, Ny))

	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	C[np.where( np.abs(C[:,:,i]) < 10) ] = 0
	C[np.where( abs(u_x[:, :]) < 1e-8)] = -100
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(C[:,:,i], levels=np.linspace(-0.1, 500, 25))
	ax.axis("equal")
	print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), np.argmax(C[:,:,i]))

fig,ax = plt.subplots()

interval = 5 #in ms     
ani = animation.FuncAnimation(fig,animate_back,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
