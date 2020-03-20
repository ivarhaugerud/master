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

Sx = 108
Sy = 45

Dx = 10
Dy = 25

datafiles = 300
t = np.linspace(0, 1, datafiles)

C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
  data_back = np.loadtxt("../data/step_0403_C_"+str(i)+"_step.txt")
  C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

  data_front = np.loadtxt("../data/step_0403_C_"+str(i)+"_front.txt")
  C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

#print(np.sum(C_back[:,:, 0]), np.sum(C_front[:,:, 0]))
C_back  /= np.sum(np.sum(C_back[ :,:, -1]))
C_front /= np.sum(np.sum(C_front[:,:, -1]))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.figure(1)
plt.plot(t, C_back[Dx, Dy, :]/max(C_back[Dx, Dy, -100:]), label="Step ")
plt.plot(t, C_front[Sx, Sy, :]/max(C_front[Sx, Sy, -100:]), label="Pulse")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Normalized concentration", fontsize=14)
plt.axis([0.6, 1.01, -0.05, 1.05])
plt.legend(loc="best", fontsize=12)
plt.savefig("../powerpoint/figures/step_1.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/step_1.pdf", "../powerpoint/figures/step_1.pdf"))
plt.show()

plt.figure(2)
plt.plot(t[:-1], np.diff(C_back[Dx, Dy, :])/max(np.diff(C_back[Dx, Dy, :])), label="Derivative of step")
plt.plot(t, C_front[Sx, Sy, :]/max(C_front[Sx, Sy, int(0.5*datafiles):]), "--", label="Pulse")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Normalized concentration", fontsize=14)
plt.axis([0.6, 1.01, -0.04, 1.05])
plt.legend(loc="best", fontsize=12)
plt.savefig("../powerpoint/figures/step_2.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/step_2.pdf", "../powerpoint/figures/step_2.pdf"))
plt.show()

import scipy.integrate as sci 

integrated = sci.cumtrapz(C_front[Sx,Sy, int(0.5*datafiles):])
plt.figure(2)
plt.plot(t, C_back[Dx, Dy, :]/C_back[Dx, Dy, -1], label="Step")
plt.plot(t[int(0.5*datafiles)+1:], integrated/integrated[-1], "--", label="Integrate of pulse")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Normalized concentration", fontsize=14)
plt.axis([0.6, 1.01, -0.02, 1.05])
plt.legend(loc="best", fontsize=12)
plt.savefig("../powerpoint/figures/step_3.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../powerpoint/figures/step_3.pdf", "../powerpoint/figures/step_3.pdf"))
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