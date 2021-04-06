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

Sx = 100
Sy = 36

Dx = 32
Dy = 32

Sx2 = 13
Sy2 = 26

Sx3 = 135
Sy3 = 39

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


plt.figure(3)
plt.streamplot(y_axis, x_axis, np.divide(u_y, 1), np.divide(u_x, 1))
plt.axis("equal")
plt.plot([Sy], [Sx], "o", label="1")
plt.plot([Sy2], [Sx2], "o", label="2")
plt.plot([Sy3], [Sx3], "o", label="3")

plt.plot([Dy], [Dx], "o", color=sns.color_palette()[1], label="Source")
plt.legend(loc="best", fontsize=12)
plt.savefig("../figures/reciprocal_symmetry3.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../figures/reciprocal_symmetry3.pdf", "../figures/reciprocal_symmetry3.pdf"))
plt.show()


plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
root = "../../../master_latex/results/figures/"



all_factors = np.array([1, 2, 4, 8, 16, 32, 64, 128])
factors = np.array([1, 2, 4, 8, 16, 32, 64, 128])
#f_factors = np.array([0.125, 0.25, 0.5,])

datafiles = 10000
skip = 5
datafiles = int(datafiles/skip)
T = 600000
t = np.linspace(0, 1, datafiles)
t = np.linspace(0, T, datafiles)

C_front = np.zeros((len(all_factors), datafiles))
C_back  = np.zeros((len(all_factors), datafiles))

C_front2 = np.zeros((len(all_factors), datafiles))
C_back2  = np.zeros((len(all_factors), datafiles))

C_front3 = np.zeros((len(all_factors), datafiles))
C_back3  = np.zeros((len(all_factors), datafiles))
#f_factors_file = np.array([4,3,2])
#for i in range(datafiles):
#	C_back[i] = np.reshape(np.loadtxt("../data/heat_vary_Pe_factor1_31_03_1_C_"+  str(i*skip) + "_matter.txt"), (Nx,Ny))[Dx, Dy]

for i in range(len(factors)):
	for j in range(datafiles):
		C_front[i, j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_01_04_" + str(factors[i]) + "_C_" + str(j*skip) + "_heat.txt"), (Nx,Ny))[Sx, Sy]
		C_back[i,  j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_01_04_" + str(factors[i]) + "_C_"+  str(j*skip) + "_matter.txt"), (Nx,Ny))[Dx, Dy]

		C_front2[i, j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_02_04_" + str(factors[i]) + "_C_" + str(j*skip) + "_heat.txt"), (Nx,Ny))[Sx2, Sy2]
		C_back2[i,  j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_02_04_" + str(factors[i]) + "_C_"+  str(j*skip) + "_matter.txt"), (Nx,Ny))[Dx, Dy]

		C_front3[i, j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_03_04_" + str(factors[i]) + "_C_" + str(j*skip) + "_heat.txt"), (Nx,Ny))[Sx3, Sy3]
		C_back3[i,  j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_03_04_" + str(factors[i]) + "_C_"+  str(j*skip) + "_matter.txt"), (Nx,Ny))[Dx, Dy]

#for i in range(len(f_factors)):
#	for j in range(datafiles):
#		C_front[i, j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_01_04_f" + str(f_factors_file[i]) + "_C_" + str(j*skip) + "_heat.txt"), (Nx,Ny))[Sx, Sy]
#		C_back[i,  j]  = np.reshape(np.loadtxt("../data/heat_vary_Pe_01_04_f" + str(f_factors_file[i]) + "_C_"+  str(j*skip) + "_matter.txt"), (Nx,Ny))[Dx, Dy]
C_front  -= background_T
C_front2 -= background_T
C_front3 -= background_T

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)
dif  = np.zeros(len(all_factors))
dif2 = np.zeros(len(all_factors))
dif3  = np.zeros(len(all_factors))

Pe = Pe/all_factors
hw = int(datafiles*0.3)

plt.figure(1)
for i in range(len(all_factors)):
	#plt.plot(t/T, C_front2[i, :], color="C3", label=r"$C_{B_{heat}}(\mathbf{x}_A, t)$")
	#plt.plot(t/T, C_back2[i,:], label=r"$C_{A_{matter}}(\mathbf{x}_B, t)$")

	#plt.plot(t/T, C_front3[i, :], color="C3", label=r"$C_{B_{heat}}(\mathbf{x}_A, t)$")
	#plt.plot(t/T, C_back3[i,:], label=r"$C_{A_{matter}}(\mathbf{x}_B, t)$")

	#plt.plot(t[hw]/T, 0, "o")

	plt.xlabel(r"Time [$T_{max}$]",         fontsize=8)
	plt.ylabel(r"Normalized Concentration", fontsize=8)
	plt.tick_params(axis='both', which='major', labelsize=8)
	plt.tick_params(axis='both', which='minor', labelsize=8)
	#plt.axis([0.2, 1.05, -0.05, 1.03])
	plt.legend(loc="best", fontsize=8)
	filename = root + "heat_and_matter.pdf"
	#plt.savefig(filename, bbox_inches="tight")
	#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
	#plt.show()
	print(Pe[i])
	#plt.plot(t, C_front[i,:]-C_back[i,:])
	dif3[i] = scp.trapz(abs(C_front3[i,hw:]-C_back3[i,hw:]), t[hw:])/scp.trapz(abs(C_back3[i,hw:]), t[hw:])
	dif2[i] = scp.trapz(abs(C_front2[i,hw:]-C_back2[i,hw:]), t[hw:])/scp.trapz(abs(C_back2[i,hw:]), t[hw:])
	dif[i]  = scp.trapz(abs(C_front[i,hw:]-C_back[i,hw:]), t[hw:])/scp.trapz(abs(C_back[i,hw:]), t[hw:])
#plt.show()


plt.figure(2)
#plt.yscale("log")
plt.plot(Pe[:-1], dif[:-1],   "o", markersize=3, label="short")
plt.plot(Pe[:-1], dif3[:-1],  "o", markersize=3, label="medium")
plt.plot(Pe[1:-1], dif2[1:-1],  "o", markersize=3, label="long")
plt.xscale("log")
plt.xlabel(r"Peclet number [Pe]", fontsize=8)
plt.ylabel(r"Integrated consentration difference [Error]", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=1)
filename = root + "vary_Pe_heat.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()