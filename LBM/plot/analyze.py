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

Nx = 102
Ny = 64

Sx = 10
Sy = 25

Dx = 60
Dy = 25

datafiles = 900
t = np.linspace(0, 1, datafiles)

C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
	data_back = np.loadtxt("../data/1307_reciproc_5_oscls_C_"+str(i)+"_front.txt")
	C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

	data_front = np.loadtxt("../data/1307_reciproc_5_oscls_C_"+str(i)+"_back.txt")
	C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

C_back  /= np.sum(np.sum(C_back[:,:, 0]))
C_front /= np.sum(np.sum(C_front[:,:, 0]))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

C_back[Dx, Dy, :20] = 0 
C_front[Sx, Sy, :20] = 0 

correct_points = 8
dx = 0.01

plt.figure(1)
#plt.title("Change of factor 1.05", fontsize=16)
plt.plot(t, C_back[Dx, Dy, :], label=r"$C_A(\mathbf{x}_B, t)$")
plt.plot(t, C_front[Sx, Sy, :], label="$C_B(\mathbf{x}_A, t)$")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
for i in range(correct_points):
	if  i % 2 != 0:
		plt.fill_between([i/correct_points-dx, i/correct_points+dx], [0, 0], [max(C_back[Dx, Dy, :]), max(C_back[Dx, Dy, :])], color="green", alpha=0.3)

plt.ylabel(r"Normalized concentration", fontsize=14)
#plt.axis([0.19, 0.61, -0.00002, 0.00065])
plt.legend(loc="best", fontsize=14)
plt.savefig("../powerpoint/figures/reciprocal_relation.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/reciprocal_relation.pdf", "../powerpoint/figures/reciprocal_relation.pdf"))
plt.show()

"""
for i in range(datafiles):
	data_back = np.loadtxt("../data/1102reciproc_2_C_"+str(i)+"_back.txt")
	C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

	data_front = np.loadtxt("../data/1102reciproc_2_C_"+str(i)+"_front.txt")
	C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

C_back  /= np.sum(np.sum(C_back[:,:, 0]))
C_front /= np.sum(np.sum(C_front[:,:, 0]))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.figure(2)
plt.title("Change of factor 2.00", fontsize=16)
plt.plot(t, C_back[Dx, Dy, :]*1e3, label="Return")
plt.plot(t, C_front[Sx, Sy, :]*1e3, label="Original")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Normalized concentration $\times 10^3$", fontsize=14)
#plt.axis([0.2, 1.05, -0.02, 0.4])
plt.legend(loc="best", fontsize=12)
plt.show()
#plt.savefig("../figures/reciprocal_symmetry2.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/reciprocal_symmetry2.pdf", "../figures/reciprocal_symmetry2.pdf"))

u = np.loadtxt("../data/0602reciproc_u.txt")
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
U = np.sum(length)/(Nx*Ny - boundary_size)
average_disc_diameter = 8
visc = (2-0.5)/3

Re = U*average_disc_diameter/visc

tau_g  = 0.50 + 6*pow(10,-5)
D = (tau_g-0.5)/3
Pe = average_disc_diameter*U/D

print("Reynolds number: ", Re)
print("Peclet number: ", Pe)


tol = 1e-10
u = np.loadtxt("../data/0602reciproc_105_to_u.txt")
u_x = u[0, :]
u_y = u[1, :]

u1 = np.loadtxt("../data/0602reciproc_105_back_u.txt")
u_x1 = u1[0, :]
u_y1 = u1[1, :]
"""
