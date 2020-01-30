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
C_r  = np.zeros((datafiles, 2))
"""
for i in range(datafiles):
	data = np.loadtxt("../data/C_"+str(i)+"_back.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	m = np.sum(np.sum(C[:,:,i]))
	C_r[i, 0] = np.sqrt((np.sum( np.dot(np.square(x_axis), C[:,:,i]))/m)**2 - (np.sum( np.dot(x_axis, C[:,:,i]))/m)**2)
	C_r[i, 1] = np.sqrt((np.sum( np.dot(C[:,:,i], np.square(y_axis)))/m)**2 - (np.sum( np.dot(C[:,:,i], y_axis))/m)**2)

plt.plot(C_r[:, 1], label="y")
plt.plot(C_r[:, 0], label="x")
plt.legend(loc="best")
plt.show()


for i in range(datafiles):
	data = np.loadtxt("../data/C_"+str(i)+"_front.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	m = np.sum(np.sum(C[:,:,i]))
	C_r[i, 0] = np.sqrt((np.sum( np.dot(np.square(x_axis), C[:,:,i]))/m)**2 - (np.sum( np.dot(x_axis, C[:,:,i]))/m)**2)
	C_r[i, 1] = np.sqrt((np.sum( np.dot(C[:,:,i], np.square(y_axis)))/m)**2 - (np.sum( np.dot(C[:,:,i], y_axis))/m)**2)

plt.plot(C_r[:, 1], label="y")
plt.plot(C_r[:, 0], label="x")
plt.legend(loc="best")
plt.show()
"""
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
	ax.contourf(C[:,:,i], levels=np.linspace(-0.1, 100, 25))
	ax.axis("equal")
	print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), np.argmax(C[:,:,i]))

fig,ax = plt.subplots()

interval = 5 #in ms     
ani = animation.FuncAnimation(fig,animate_back,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
