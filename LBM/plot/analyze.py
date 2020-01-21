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

#C_drain = np.loadtxt("../data/sourcefront.txt")

Nx = 140
Ny = 64

Sx = 34
Sy = 34

Dx = 32
Dy = 32

datafiles = 100

C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
	data_back = np.loadtxt("../data/C_"+str(i)+"_back.txt")
	C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

	data_front = np.loadtxt("../data/C_"+str(i)+"_front.txt")
	C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))
	print(np.sum(np.sum(C_back[:,:, i])), np.sum(np.sum(C_front[:,:,i])))
x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.plot( C_back[Dx, Dy, :])
plt.plot(C_front[Sx, Sy, :], "o")
plt.show()
#print("Argument of global maxima: ", np.unravel_index(np.argmax(C, axis=None), C.shape))
#print("Center of mass: ", np.mean(np.dot(x_axis, C)), np.mean(np.dot(C, y_axis)))

#plt.imshow(y_axis, x_axis, C)
#plt.show()

x_, y_ = np.meshgrid( y_axis, x_axis)
fig = plt.figure()
plt.ylabel(r"$x$", fontsize=14)
plt.xlabel(r"$y$", fontsize=14)

print("mass: ", np.sum(np.sum(C)))
ax1 = plt.contourf(y_,x_, C)#, levels=np.linspace(0, np.max(np.max(C)), 30))
cbar = fig.colorbar(ax1)
cbar.ax.set_ylabel(r'Concentration $C$', fontsize=14)
#plt.xscale('log')
#plt.savefig("../figures/relative_diff.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/relative_diff.pdf", "../figures/relative_diff.pdf"))
plt.show()