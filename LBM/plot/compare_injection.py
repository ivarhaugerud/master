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

Dx = 25
Dy = 25
"""
C = np.zeros((Nx, Ny, 300))
for i in range(300):
	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_front.txt")
	C[:, :, i] = (np.reshape(data, (Nx, Ny)))
plt.plot(np.max(C[100,:,:], axis=0))
plt.show()
"""
datafiles = 300
t = np.linspace(0, 1, datafiles)
C = np.zeros((Nx, Ny, datafiles, 6))
x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

for i in range(datafiles):
	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_1_0.txt")
	C[:, :, i, 0] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_5.txt")
	C[:, :, i, 1] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_8.txt")
	C[:, :, i, 2] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_0_9.txt")
	C[:, :, i, 3] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_single_maxima.txt")
	C[:, :, i, 4] = (np.reshape(data, (Nx, Ny)))

	data = np.loadtxt("../data/peak_1103_C_"+str(i)+"_maxima.txt")
	C[:, :, i, 5] = (np.reshape(data, (Nx, Ny)))

C[:, :, :, 0] /= np.sum(np.sum(C[:,:,-1,0]))
C[:, :, :, 1] /= np.sum(np.sum(C[:,:,-1,1]))
C[:, :, :, 2] /= np.sum(np.sum(C[:,:,-1,2]))
C[:, :, :, 3] /= np.sum(np.sum(C[:,:,-1,3]))
C[:, :, :, 4] /= np.sum(np.sum(C[:,:,-1,4]))
C[:, :, :, 5] /= np.sum(np.sum(C[:,:,-1,5]))

plt.figure(1)
plt.plot(t, C[Dx, Dy, :, 0]*100, label="no cutoff")
plt.plot(t, C[Dx, Dy, :, 1]*100, label="0.5 cutoff")
plt.plot(t, C[Dx, Dy, :, 2]*100, label="0.8 cutoff")
plt.plot(t, C[Dx, Dy, :, 3]*100, label="0.9 cutoff")
plt.plot(t, C[Dx, Dy, :, 4]*100, label="singe maximum")
plt.plot(t, C[Dx, Dy, :, 5]*100, label="max for each pos")

plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Concentration %", fontsize=14)
#plt.axis([0.03, 1.05, -0.02, 0.4])
plt.legend(loc="best", fontsize=13)
plt.show()

C_var  = np.zeros((datafiles, len(C[0,0,0,:])))
"""
plt.figure(2)
for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		summed = np.sum(C[:, :, i, f], axis=1)
		C_x = x_axis*summed
		C_xx = x_axis*C_x
		m = np.trapz(summed)
		C_var[i, f] = np.trapz(C_xx)/m - (np.trapz(C_x)/m)**2
"""
import scipy.integrate as sci 

plt.figure(2)
for f in range(len(C[0,0,0,:])):
	for i in range(datafiles):
		C_var[i, f] = sci.dblquad(C[:,:,i,f], Nx, Ny, x_axis, y_axis)
		#C_x = x_axis*summed
		#C_xx = x_axis*C_x
		#m = np.trapz(summed)
		#C_var[i, f] = np.trapz(C_xx)/m - (np.trapz(C_x)/m)**2

plt.plot(t, C_var[:, 0], label="no cutoff")
plt.plot(t, C_var[:, 1], label="0.5 cutoff")
plt.plot(t, C_var[:, 2], label="0.8 cutoff")
plt.plot(t, C_var[:, 3], label="0.9 cutoff")
plt.plot(t, C_var[:, 4], label="single maximum")
plt.plot(t, C_var[:, 5], label="max for each pos")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Variance of concentration", fontsize=14)
plt.legend(loc="best", fontsize=13)
plt.show()
#plt.savefig("../figures/reciprocal_symmetry1.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/reciprocal_symmetry1.pdf", "../figures/reciprocal_symmetry1.pdf"))

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

