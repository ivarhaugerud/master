import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import matplotlib
import pandas as pd 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
root = "../../../master_latex/results/"

Nx = 140
Ny = 64

Sx = 100
Sy = 42

Dx = 34
Dy = 32

datafiles = 300
t = np.linspace(0, 1, datafiles)
C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
	data_back = np.loadtxt("../data/1102reciproc_2_C_"+str(i)+"_back.txt")
	C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

	data_front = np.loadtxt("../data/1102reciproc_2_C_"+str(i)+"_front.txt")
	C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

#C_back  /= np.sum(np.sum(C_back[:,:, 0]))
#C_front /= np.sum(np.sum(C_front[:,:, 0]))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.figure(2)
plt.plot(t,  C_back[Dx, Dy, :]/max(C_back[Dx,Dy, 30:]), label=r"$C_{A}(\mathbf{x}_B, t)$")
plt.plot(t, C_front[Sx, Sy, :]/max(C_back[Dx,Dy, 30:]), color="C3", label=r"$C_{B}(\mathbf{x}_A, t)$")
plt.axis([0.10, 1.025, -0.025, 1.05])
plt.xlabel(r"Time [$T_{max}$]", fontsize=8)
plt.ylabel(r"Normalized Concentration", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=3)
filename = root + "same_Pe.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


plt.figure(2)
plt.plot(t*2, C_back[Dx, Dy, :]/max(C_back[Dx,Dy,30:]), label=r"$C_{A}(\mathbf{x}_B, t/\alpha )$")
plt.plot(t,  C_front[Sx, Sy, :]/max(C_back[Dx,Dy,30:]),  "--", color="C3", label=r"$C_{B}(\mathbf{x}_A, t)$")
plt.axis([0.10, 2.05, -0.025, 1.05])
plt.xlabel(r"Time [$T_{max}$]", fontsize=8)
plt.ylabel(r"Normalized Concentration", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=3)
filename = root + "same_Pe_rescaled.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
"""

u = np.loadtxt("../data/1102reciproc_2_to_u.txt")
u_x = u[0, :]
u_y = u[1, :]

u_x = np.reshape(u_x, (Nx, Ny))
u_y = np.reshape(u_y, (Nx, Ny))
length = np.sqrt(u_y*u_y + u_x*u_x)


boundary_size = np.sum(np.where(abs(u) < 1e-8)[0])
U = np.sum(length)/(Nx*Ny - boundary_size)
average_disc_diameter = 8
visc = (2-0.5)/3

Re = U*average_disc_diameter/visc

tau_g  = 0.50 + 6*pow(10,-5)
D = (tau_g-0.5)/3
Pe = average_disc_diameter*U/D

print("Reynolds number: ", Re)
print("Peclet number: ", Pe)

"""