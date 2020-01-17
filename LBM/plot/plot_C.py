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

C = np.loadtxt("../data/final_C.txt")

Nx = 50
Ny = 50


x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

C = np.reshape(C, (Nx, Ny))

#print("Argument of global maxima: ", np.unravel_index(np.argmax(C, axis=None), C.shape))
#print("Center of mass: ", np.mean(np.dot(x_axis, C)), np.mean(np.dot(C, y_axis)))

#plt.imshow(y_axis, x_axis, C)
#plt.show()

x_, y_ = np.meshgrid( x_axis, y_axis)
fig = plt.figure()
plt.ylabel(r"$x$", fontsize=14)
plt.xlabel(r"$y$", fontsize=14)

ax1 = plt.contourf(y_,x_, C, levels=np.linspace(0, np.max(np.max(C)), 30))
cbar = fig.colorbar(ax1)
cbar.ax.set_ylabel(r'Concentration $C$', fontsize=14)
#plt.xscale('log')
#plt.savefig("../figures/relative_diff.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/relative_diff.pdf", "../figures/relative_diff.pdf"))
plt.show()