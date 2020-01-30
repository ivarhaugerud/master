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
#var = np.zeros(datafiles)
var_in = np.zeros(datafiles)
Nx = 256
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

C = np.zeros((Nx, Ny, datafiles))
C_r  = np.zeros(datafiles)

for i in range(datafiles):
	data = np.loadtxt("../data/C_"+str(i)+"_back.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))

	mean_x, mean_y = (ndimage.measurements.center_of_mass(C[:,:,i]))

	for x in range(Nx):
		for y in range(Ny):
			C_r[i] += C[x, y, i]*np.sqrt((x-mean_x)**2 + (y-mean_y)**2)

plt.plot(C_r)
plt.show()



max_C = np.max(np.max(np.max(C)))

fig,ax = plt.subplots()

def animate_forward(i):
	fig.suptitle(str(i))	
	data = np.loadtxt("../data/C_"+str(i)+"_front.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(C[:,:,i], levels=np.linspace(1e-5, 0.5, 25))
	ax.axis("equal")
	print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), np.argmax(C[:,:,i]))

fig,ax = plt.subplots()

interval = 5 #in ms     
ani = animation.FuncAnimation(fig,animate_forward,datafiles, interval=interval, blit=False, repeat=False)
plt.show()

def animate_back(i):
	fig.suptitle(str(i))
	data = np.loadtxt("../data/C_"+str(i)+"_back.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(C[:,:,i], levels=np.linspace(0, 100, 25))
	ax.axis("equal")
	print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), np.argmax(C[:,:,i]))

fig,ax = plt.subplots()

interval = 5 #in ms     
ani = animation.FuncAnimation(fig,animate_back,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
