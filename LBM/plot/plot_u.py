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

#u = np.loadtxt("../data/final_vel.txt")
u = np.loadtxt("../data/0203heat_heat_u.txt")
u_x = u[0, :]
u_y = u[1, :]

Nx = 140
Ny = 64

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

u_x = np.reshape(u_x, (Nx, Ny))
u_y = np.reshape(u_y, (Nx, Ny))


length = np.sqrt(u_y*u_y + u_x*u_x)
print(np.shape(length))

"""
plt.quiver(y_axis, x_axis, np.divide(u_y, 1), np.divide(u_x, 1), length)
plt.axis("equal")
plt.show()
"""
"""
#stream_points[:, 0] = x_streampoints
#stream_points[:, 1] = 0

plt.streamplot(y_axis, x_axis, u_y, u_x)#, start_points=stream_points, density=30, color=sns.color_palette()[1])
"""
"""
x_streampoints = np.arange(0, 63, 4)
stream_points = np.zeros((len(x_streampoints), 2))
stream_points[:, 0] = x_streampoints
stream_points[:, 1] = 0
plt.streamplot(y_axis, x_axis, u_y, u_x, start_points=stream_points, density=15, color="k", linewidth=1)

x_streampoints = np.arange(0, 63, 4)
stream_points = np.zeros((len(x_streampoints), 2))
stream_points[:, 0] = x_streampoints
stream_points[:, 1] = 130
plt.streamplot(y_axis, x_axis, u_y, u_x, start_points=stream_points, density=15, color="k", linewidth=1)
"""
x_streampoints = np.arange(0, 63, 0.4)
stream_points = np.zeros((len(x_streampoints), 2))
stream_points[:, 0] = x_streampoints
stream_points[:, 1] = 84
plt.streamplot(y_axis, x_axis, u_y, u_x, start_points=stream_points, density=30, color="k", linewidth=1)

plt.axis("equal")
plt.show()
"""
x_streampoints = np.arange(0, 100, 2.5)
stream_points = np.zeros((len(x_streampoints), 2))
stream_points[:, 1] = x_streampoints
stream_points[:, 0] = -30

plt.streamplot(y_axis, x_axis, u_x, u_y, start_points=stream_points, density=30, color=sns.color_palette()[1])

plt.show()
"""