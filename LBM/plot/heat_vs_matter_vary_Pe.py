import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import matplotlib
import pandas as pd 
import scipy.integrate as scp 

u_x = np.loadtxt("../data/heat_vary_Pe_factor1_31_03_1___ux.txt")
u_y = np.loadtxt("../data/heat_vary_Pe_factor1_31_03_1___uy.txt")

Nx = 140
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

boundary_size = len(np.where(abs(u_x) < 1e-10)[0])

u_x = np.reshape(u_x, (Nx, Ny))
u_y = np.reshape(u_y, (Nx, Ny))
length = np.sqrt(u_y*u_y + u_x*u_x)
U = abs(np.sum(length)/(Nx*Ny - boundary_size))

average_disc_diameter = 8
visc = (2-0.5)/3
Re = U*average_disc_diameter/visc

tau_g  = 0.50 + 6*pow(10,-5)
D = (tau_g-0.5)/3
Pe = average_disc_diameter*U/D
background_T = 0.0001

M = 0.1-background_T
print("Reynolds number: ", Re)
print("Peclet number: ", Pe)

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
root = "../../../master_latex/results/"


Sx = 100
Sy = 36

Dx = 32
Dy = 32


factors = np.array([1, 2, 4, 8, 16, 32, 64])
datafiles = 10000
skip = 100
datafiles = int(datafiles/skip)
T = 600000
t = np.linspace(0, 1, datafiles)
t = np.linspace(0, T, datafiles)

C_front = np.zeros((len(factors), datafiles))
C_back  = np.zeros((len(factors), datafiles))

#for i in range(datafiles):
#	C_back[i] = np.reshape(np.loadtxt("../data/heat_vary_Pe_factor1_31_03_1_C_"+  str(i*skip) + "_matter.txt"), (Nx,Ny))[Dx, Dy]

for i in range(len(factors)):
	for j in range(datafiles):
		C_front[i, j] = np.reshape(np.loadtxt("../data/heat_vary_Pe_factor1_31_03_" + str(factors[i]) + "_C_" + str(j*skip) + "_heat.txt"), (Nx,Ny))[Sx, Sy]
		C_back[i, j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_factor1_31_03_" + str(factors[i]) + "_C_"+  str(j*skip) + "_matter.txt"), (Nx,Ny))[Dx, Dy]

C_front -= background_T

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)
dif = np.zeros(len(factors))
Pe = Pe/factors
hw = int(datafiles*0.1)

plt.figure(1)
for i in range(len(factors)):
	plt.plot(t/T, C_front[i, :]/max(C_back[i, int(datafiles/2):]), color="C3", label=r"$C_{B_{heat}}(\mathbf{x}_A, t)$")
	plt.plot(t/T, C_back[i,:]/max(C_back[i, int(datafiles/2):]), label=r"$C_{A_{matter}}(\mathbf{x}_B, t)$")
	plt.plot(t[hw], 0, "o")

	plt.xlabel(r"Time [$T_{max}$]",         fontsize=8)
	plt.ylabel(r"Normalized Concentration", fontsize=8)
	plt.tick_params(axis='both', which='major', labelsize=8)
	plt.tick_params(axis='both', which='minor', labelsize=8)
	plt.axis([0.2, 1.05, -0.05, 1.03])
	plt.legend(loc="best", fontsize=8)
	filename = root + "heat_and_matter.pdf"
	#plt.savefig(filename, bbox_inches="tight")
	#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
	#plt.show()
	print(Pe)
	#plt.plot(t, C_front[i,:]-C_back[i,:])
	dif[i] = scp.trapz(abs(C_front[i,hw:]-C_back[i,hw:]), t[hw:])/scp.trapz(abs(C_back[i,hw:]), t[hw:])
	plt.show()

plt.yscale("log")
plt.plot(Pe, dif, "o", markersize=3)
plt.xscale("log")
plt.xlabel(r"Peclet number [Pe]", fontsize=8)
plt.ylabel(r"Integrated consentration difference ", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=1)
filename = root + "vary_Pe_heat.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()