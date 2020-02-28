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

datafiles = 150-50
t = np.linspace(0, 1, datafiles)

C_front = np.zeros((Nx, Ny, datafiles))
C_back  = np.zeros((Nx, Ny, datafiles))

for i in range(datafiles):
  data_back = np.loadtxt("../data/step_C_"+str(i+49)+"_back.txt")
  C_back[:, :, i] = (np.reshape(data_back, (Nx, Ny)))

  data_front = np.loadtxt("../data/step_C_"+str(i+49)+"_front.txt")
  C_front[:, :, i] = (np.reshape(data_front, (Nx, Ny)))

#print(np.sum(C_back[:,:, 0]), np.sum(C_front[:,:, 0]))
C_back  /= np.sum(np.sum(C_back[ :,:, -1]))
C_front /= np.sum(np.sum(C_front[:,:, -1]))

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

plt.figure(1)
plt.plot(t, C_back[Dx, Dy, :]*1e3, label="Return")
plt.plot(t, C_front[Sx, Sy, :]*1e3, "--", label="Original")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Normalized concentration $\times 10^3$", fontsize=14)
#plt.axis([0.03, 1.05, -0.02, 0.4])
plt.legend(loc="best", fontsize=12)
plt.show()

plt.figure(2)
plt.plot(t[:-1], np.diff(C_back[Dx, Dy, :])/max(np.diff(C_back[Dx, Dy, :])), label="Derivative of return")
plt.plot(t, C_front[Sx, Sy, :]/max(C_front[Sx, Sy, :]), "--", label="Original")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Normalized concentration $\times 10^3$", fontsize=14)
#plt.axis([0.03, 1.05, -0.02, 0.4])
plt.legend(loc="best", fontsize=12)
plt.show()

import scipy.integrate as sci 

integrated = sci.cumtrapz(C_front[Sx,Sy,:])
plt.figure(2)
plt.plot(t, C_back[Dx, Dy, :]/C_back[Dx, Dy, -1], label="Return")
plt.plot(t[:-1], integrated/integrated[-1], "--", label="Integrate of initial")
plt.xlabel(r"Time [$T_{max}$]", fontsize=14)
plt.ylabel(r"Normalized concentration $\times 10^3$", fontsize=14)
#plt.axis([0.03, 1.05, -0.02, 0.4])
plt.legend(loc="best", fontsize=12)
plt.show()
#plt.savefig("../figures/reciprocal_symmetry1.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/reciprocal_symmetry1.pdf", "../figures/reciprocal_symmetry1.pdf"))

u = np.loadtxt("../data/step_step_back_u.txt")
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
"""
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import matplotlib
import pandas as pd 
import matplotlib.animation as animation
from scipy import ndimage

plt.style.use("bmh")
sns.color_palette("hls", 1)

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

datafiles = 200
var_in = np.zeros(datafiles)
Nx = 140
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)
C = np.zeros((Nx, Ny, datafiles))

u = np.loadtxt("../data/1902heat_heat_u.txt")
u_x = np.reshape(u[0, :], (Nx, Ny))

def animate_forward(i):
  fig.suptitle(str(i))
  data = np.loadtxt("../data/step_C_"+str(i)+"_front.txt")

  C[:, :, i] = (np.reshape(data, (Nx, Ny)))
  print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), C.flatten()[np.argmax(C[:,:,i])])

  C[np.where( np.abs(C[:,:,i]) < 0.00001) ] = 0
  C[np.where( abs(u_x[:, :]) < 1e-8)] = -1
  ax.clear()
  #ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
  ax.contourf(C[:,:,i], levels=np.linspace(-0.0001, 0.01, 25))
  ax.axis("equal")

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig,ax = plt.subplots()
interval = 1 #in ms     
ani = animation.FuncAnimation(fig,animate_forward,datafiles, interval=interval, blit=False, repeat=False)
#ani.save('im.mp3')

plt.show()

def animate_back(i):
  fig.suptitle(str(i))
  data = np.loadtxt("../data/step_C_"+str(i)+"_back.txt")

  C[:, :, i] = (np.reshape(data, (Nx, Ny)))
  print(np.max(np.max(C[:,:,i])), np.sum(np.sum(C[:,:,i])), np.argmax(C[:,:,i]))

  C[np.where( np.abs(C[:,:,i]) < 0.00001) ] = 0
  C[np.where( abs(u_x[:, :]) < 1e-8)] = -1
  ax.clear()
  #ax.contourf(1-np.exp(-40*C[:,:,i]), cmap='Greys', levels=np.linspace(0, 1, 20))
  ax.contourf(C[:,:,i], levels=np.linspace(-0.001, 1000, 25))
  ax.axis("equal")

fig,ax = plt.subplots()

interval = 1 #in ms     
ani = animation.FuncAnimation(fig,animate_back,datafiles, interval=interval, blit=False, repeat=False)

plt.show()
"""
