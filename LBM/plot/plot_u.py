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

u = np.loadtxt("../data/final_vel.txt")
u_x = u[0, :]
u_y = u[1, :]

Nx = 10
Ny = 80

f = 5*1e-8
mu = (2-0.5)/3
a = int(Ny/2)
U = a*a*f/(3*mu)

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(-Ny, Ny, Ny)

u_x = np.reshape(u_x, (Nx, Ny))
u_y = np.reshape(u_y, (Nx, Ny))

measured = np.mean(u_x, axis=0)
mean_U = np.mean(measured)

print(U, mean_U)
analytic = 3*U*(1 - (y_axis/Ny)**2)/2
measured = measured

plt.plot(y_axis, measured)
plt.plot(y_axis, analytic)
plt.show()

plt.plot((measured-analytic))
plt.show()

plt.quiver(y_axis, x_axis, u_y, u_x)
plt.show()