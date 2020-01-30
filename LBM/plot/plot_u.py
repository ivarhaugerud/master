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

Nx = 256
Ny = 64

#for i in range(len(u[0, :])):
#	print(u[0, i])
print(np.max(np.max(u)))

print(np.where( abs(u[0, :]) < 1e-8))
f = 5*1e-8
mu = (2-0.5)/3
a = int(Ny/2)-1
U = a*a*f/(3*mu)

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(-a, a, Ny)

u_x = np.reshape(u_x, (Nx, Ny))
u_y = np.reshape(u_y, (Nx, Ny))

measured = np.mean(u_x, axis=0)
mean_U = np.mean(measured)

analytic = 3*U*(1 - (y_axis/a)**2)/2
measured = measured

length = np.sqrt(u_y*u_y + u_x*u_x)
print(np.shape(length))
plt.quiver(y_axis, x_axis, np.divide(u_y, 1), np.divide(u_x, 1), length)
plt.axis("equal")
plt.show()
"""

x_streampoints = np.arange(0, 100, 2.5)
stream_points = np.zeros((len(x_streampoints), 2))
stream_points[:, 1] = x_streampoints
stream_points[:, 0] = 30

plt.streamplot(y_axis, x_axis, u_x, u_y, start_points=stream_points, density=30, color=sns.color_palette()[1])

x_streampoints = np.arange(0, 100, 2.5)
stream_points = np.zeros((len(x_streampoints), 2))
stream_points[:, 1] = x_streampoints
stream_points[:, 0] = -30

plt.streamplot(y_axis, x_axis, u_x, u_y, start_points=stream_points, density=30, color=sns.color_palette()[1])

plt.show()
"""