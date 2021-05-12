import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import seaborn as sns
import os
import matplotlib.ticker as tick


plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

root = "../../../master_latex/results/figures/"
Map = matplotlib.cm.get_cmap('Spectral_r')


u_x = np.loadtxt("../data/benchmark_testing_benchmark_pousielle_ux.txt")
u_y = np.loadtxt("../data/benchmark_testing_benchmark_pousielle_uy.txt")

Nx = 512
Ny = 128

x_axis = np.linspace(0, Nx-1, Nx)
y_axis = np.linspace(0, Ny-1, Ny)

boundary_size = len(np.where(abs(u_x) < 1e-10)[0])
indxs = np.where(abs(u_x) < 1e-10)[0]
indys = np.where(abs(u_x) < 1e-10)[1]

print(np.shape(indxs))
print(indxs)

u_x = np.reshape(u_x, (Nx, Ny))
u_y = np.reshape(u_y, (Nx, Ny))
length = np.sqrt(u_y*u_y + u_x*u_x)
U = abs(np.sum(length)/(Nx*Ny - boundary_size))

average_disc_diameter = 15
visc = (2-0.5)/3
Re = U*average_disc_diameter/visc

print("Reynolds number: ", Re)

p = 10
stream_points    = np.array(list(zip( np.ones(p)*71, np.linspace(1, 127, p))))
stream_points3    = np.array(list(zip( np.ones(p)*232, np.linspace(1, 127, p))))
stream_points2   = np.array(list(zip( np.ones(p)*340, np.linspace(1, 127, p))))
stream_points4   = np.array(list(zip( np.ones(p)*500, np.linspace(1, 127, p))))
p
fig = plt.figure(3, figsize=(8*1.25, 3*0.8))
plt.plot(x_axis[indxs], y_axis[indys], "ko", markersize=4)
X, Y = np.meshgrid(x_axis, y_axis)
import matplotlib
plt.streamplot(x_axis, y_axis, np.transpose(u_x)/U, np.transpose(u_y)/U, start_points=stream_points,  density=1, color='k', arrowsize=1.35)
plt.streamplot(x_axis, y_axis, np.transpose(u_x)/U, np.transpose(u_y)/U, start_points=stream_points2, density=1, color='k', arrowsize=1.35)
plt.streamplot(x_axis, y_axis, np.transpose(u_x)/U, np.transpose(u_y)/U, start_points=stream_points3, density=1, color='k', arrowsize=1.35)
plt.streamplot(x_axis, y_axis, np.transpose(u_x)/U, np.transpose(u_y)/U, start_points=stream_points4, density=1, color='k', arrowsize=1.35)
CS = plt.contourf(X, Y, np.transpose(length)/U, 10, cmap=Map)
cbar = plt.colorbar(CS, fraction=0.02725, pad=0.01)
cbar.ax.set_ylabel('Velocity $[U]$', fontsize=12)
cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
plt.axis("equal")
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='minor', labelsize=12)
plt.ylabel(r"Vertical position", fontsize=12)
plt.xlabel(r"Horizontal position", fontsize=12)
filename = root + "visualize_u.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
