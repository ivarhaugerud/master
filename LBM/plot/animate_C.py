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

datafiles = 900
var_in = np.zeros(datafiles)
Nx = 102
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)
C = np.zeros((Nx, Ny, datafiles))

#u = np.loadtxt("../data/0303heat_heat_u.txt")
#u_x = np.reshape(u[0, :], (Nx, Ny))

def animate_forward(i):
	fig.suptitle(str(i))
	data = np.loadtxt("../data/1307_reciproc_5_oscls_C_"+str(i)+"_front.txt")

	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	C[:, :, i] /= np.sum(np.sum(C[:,:,i]))
	#print(np.sum(np.sum(C[100,:,i])))
	#print(np.max(np.max(np.log10(abs(C[:,:,i]+1)))), np.min(np.min(np.log10(abs(C[:,:,i]+1)))))

	#C[np.where( (C[:,:,i]) < 1e-7) ] = 1e-4
	#C[np.where( abs(u_x[:, :]) < 1e-8)] = -1
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(np.log10(((C[:,:,i]+1e-4))), levels=np.linspace(-4.5, -1.4, 35))
	ax.axis("equal")
	#ax.set_xlabel(r"$x$-position", fontsize=14)
	#ax.set_ylabel(r"$y$-position", fontsize=14)

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig,ax = plt.subplots(figsize=(3,6))

interval = 1 #in ms   
#Writer = animation.writers["ffmpeg"]
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

ani = animation.FuncAnimation(fig,animate_forward, datafiles, interval=interval, blit=False, repeat=False)
#ani.save('../powerpoint/figures/im.mp4')#, writer=writer)
# writer = animation.FFMpegWriter(
#     fps=15, metadata=dict(artist='Me'), bitrate=1800)
#ani.save("../powerpoint/figures/animation_single_maxima.gif")#, Writer=writer)
plt.show()

def animate_back(i):
	fig.suptitle(str(i))
	data = np.loadtxt("../data/1307_reciproc_5_oscls_C_"+str(i)+"_back.txt")

	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
	C[:, :, i] /= np.sum(np.sum(C[:,:,i]))
	#print(np.sum(np.sum(C[100,:,i])))
	#print(np.max(np.max(np.log10(abs(C[:,:,i]+1)))), np.min(np.min(np.log10(abs(C[:,:,i]+1)))))

	#C[np.where( (C[:,:,i]) < 1e-7) ] = 1e-4
	#C[np.where( abs(u_x[:, :]) < 1e-8)] = -1
	ax.clear()
	#ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
	ax.contourf(np.log10(((C[:,:,i]+1e-4))), levels=np.linspace(-4.5, -1.4, 35))
	ax.axis("equal")
	#ax.set_xlabel(r"$x$-position", fontsize=14)
	#ax.set_ylabel(r"$y$-position", fontsize=14)

fig,ax = plt.subplots()

interval = 1 #in ms     
ani = animation.FuncAnimation(fig,animate_back,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
