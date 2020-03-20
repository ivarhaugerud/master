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

Nx = 140
Ny = 64

Sx = 100
Sy = 36

Dx = 32
Dy = 32

datafiles = 300
t = np.linspace(0, 1, datafiles)

C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
	data_back = np.loadtxt("../data/0303heat_F10-7_C_"+str(i)+"_from.txt")
	C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

	data_front = np.loadtxt("../data/0303heat_F10-7_C_"+str(i)+"_to.txt")
	C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

C_back  /= np.sum(np.sum(C_back[:,:, 0]))
C_front /= np.sum(np.sum(C_front[:,:, 0]))

C_back  -= C_back[Dx, Dy, 0]
C_front -= C_front[Sx, Sy, 0]

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.figure(1)
plt.plot(t, C_front[Sx, Sy, :]*1e3, label=r"$C_A(\mathbf{x}_B, t)$")
plt.plot(t, C_back[Dx, Dy, :]*1e3, "--", label=r"$C_B(\mathbf{x}_A, t)$")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Concentration $\times 10^3$", fontsize=14)
plt.axis([0.05, 1.05, -0.003, 0.06])
plt.legend(loc="best", fontsize=14)
plt.savefig("../powerpoint/figures/heat_both_ways.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/heat_both_ways.pdf", "../powerpoint/figures/heat_both_ways.pdf"))

plt.show()
u = np.loadtxt("../data/0303heat_F10-7_heat_u.txt")
u_x = u[0, :]
u_y = u[1, :]

Nx = 140
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

u_x = np.reshape(u_x, (Nx, Ny))
u_y = np.reshape(u_y, (Nx, Ny))
length = np.sqrt(u_y*u_y + u_x*u_x)

plt.figure(3)
plt.quiver(y_axis, x_axis, np.divide(u_y, 1), np.divide(u_x, 1), length)
plt.axis("equal")
plt.plot([Sy], [Sx], "o", color="k", label="Drain")
plt.plot([Dy], [Dx], "o", color=sns.color_palette()[1], label="Source")
plt.legend(loc="best", fontsize=12)
plt.savefig("../figures/reciprocal_symmetry3.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/reciprocal_symmetry3.pdf", "../figures/reciprocal_symmetry3.pdf"))

boundary_size = np.sum(np.where(abs(u) < 1e-8)[0])
print(boundary_size)
U = np.sum(length)/(Nx*Ny - boundary_size)
average_disc_diameter = 8
visc = (2-0.5)/3

Re = U*average_disc_diameter/visc

tau_g  = 0.50 + 6*pow(10,-5)
D = (tau_g-0.5)/3
Pe = average_disc_diameter*U/D

print("Reynolds number: ", Re)
print("Peclet number: ", Pe)
plt.show()
