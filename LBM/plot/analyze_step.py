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
root = "../../../master_latex/results/figures/rough"


Nx = 140
Ny = 64

Sx = 108
Sy = 45

Dx = 10
Dy = 25

datafiles = 300
t = np.linspace(0, 1, datafiles)

C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
  data_back = np.loadtxt("../data/22_04_step_C_"+str(i)+"_step.txt")
  C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

  data_front = np.loadtxt("../data/22_04_step_C_"+str(i)+"_front.txt")
  C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.figure(1)
plt.plot(t, C_back[Dx, Dy, :]*1e3, label="Step ")
plt.plot(t, C_front[Sx, Sy, :]*1e3, color="C3", label="Pulse")
plt.axis([0.2, 1.05, -0.02, 0.75])
plt.xlabel(r"Time [$T_{max}$]", fontsize=8)
plt.ylabel(r"Concentration $[m_a\times 10^3]$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=3)
filename = root + "step_injection.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

plt.figure(2)
deriv = np.gradient(C_back[Dx, Dy, :], t)*1e3
plt.plot(t, deriv, label="Derivative of step")
plt.plot(t, C_front[Sx, Sy, :]*1e3, "--", color="C3", label="Pulse")
plt.axis([0.2, 1.05, -0.02, 0.75])
plt.xlabel(r"Time [$T_{max}$]", fontsize=8)
plt.ylabel(r"Concentration $[m_a\times 10^3]$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=1)
filename = root + "step_injection_integral.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

import scipy.integrate as sci 

integrated = sci.cumtrapz(C_front[Sx,Sy, int(0.1*datafiles):], t[int(0.1*datafiles):])
plt.figure(2)
plt.plot(t, C_back[Dx, Dy, :]*1e3, label="Step")
plt.plot(t[int(0.1*datafiles)+1:], integrated*1e3, "--", color="C3", label="Integrated pulse")
plt.axis([0.1, 1.05, -0.005, 0.14])
plt.xlabel(r"Time [$T_{max}$]", fontsize=8)
plt.ylabel(r"Concentration $[m_a\times 10^3]$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=3)
filename = root + "step_injection_derivative.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


u = np.loadtxt("../data/step_0403_step_u.txt")
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