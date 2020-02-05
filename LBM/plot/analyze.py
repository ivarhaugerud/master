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

Nx = 256
Ny = 64

Sx = 14
Sy = 32

Dx = 40
Dy = 18

datafiles = 200

C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
	data_back = np.loadtxt("../data/050220_C_"+str(i)+"_back.txt")
	C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

	data_front = np.loadtxt("../data/050220_C_"+str(i)+"_front.txt")
	C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

C_back  /= np.sum(np.sum(C_back[:,:, 0]))
C_front /= np.sum(np.sum(C_front[:,:, 0]))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.plot( C_back[Dx, Dy, :])
plt.plot(C_front[Sx, Sy, :], "o")
plt.xlabel(r"Time", fontsize=14)
plt.ylabel(r"Normalized concentration", fontsize=14)
#plt.savefig("../figures/reciprocal_symmetry.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/reciprocal_symmetry.pdf", "../figures/reciprocal_symmetry.pdf"))
plt.show()
#print("Argument of global maxima: ", np.unravel_index(np.argmax(C, axis=None), C.shape))
#print("Center of mass: ", np.mean(np.dot(x_axis, C)), np.mean(np.dot(C, y_axis)))

#plt.imshow(y_axis, x_axis, C)
#plt.show()

x_, y_ = np.meshgrid( y_axis, x_axis)
fig = plt.figure()
plt.ylabel(r"$x$", fontsize=14)
plt.xlabel(r"$y$", fontsize=14)

ax1 = plt.contourf(y_,x_, C_back[:,:,-1])#, levels=np.linspace(0, np.max(np.max(C)), 30))
cbar = fig.colorbar(ax1)
cbar.ax.set_ylabel(r'Concentration $C$', fontsize=14)
plt.show()


x_, y_ = np.meshgrid( y_axis, x_axis)
fig = plt.figure()
plt.ylabel(r"$x$", fontsize=14)
plt.xlabel(r"$y$", fontsize=14)

ax1 = plt.contourf(y_,x_, C_front[:,:,-1])#, levels=np.linspace(0, np.max(np.max(C)), 30))
cbar = fig.colorbar(ax1)
cbar.ax.set_ylabel(r'Concentration $C$', fontsize=14)
plt.show()
