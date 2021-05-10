import numpy as np
import scipy.interpolate as sci
import h5py
from scipy.interpolate import griddata
import matplotlib.pyplot as plt 
import matplotlib
import os
Map = matplotlib.cm.get_cmap('Spectral_r')

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'  
root = "../../../master_latex/results/"
base = "../data_test/vary_geometry/"
root = "../../../master_latex/results/figures/num_streamlines/"

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 5/3

#geometry parameters
epsilon = 0.4
Lx = 6.28
nu = 2.25
D = 1.0
F0 = 1000/nu
eps = 0.4
kappa  = 2*np.pi/Lx
omega = 2*np.pi/tau
dt = 0.004
tau = 3.0


dirr = "data/save_u/Lx10.47_tau3.0_eps0.4_nu2.25_D1.0_fzero0.0_fone1000.0_res150_dt0.004/" #1720
dirr = "data/save_u/Lx6.28_tau3.0_eps0.4_nu2.25_D1.0_fzero0.0_fone1000.0_res125_dt0.005/" #1720
dt = 0.005

tdata = np.loadtxt(dirr+"tdata.dat")


T = np.max(tdata[:,0])
N = int(T/dt)-1

t = tdata[:,0]

name = dirr + "u.h5"
data = h5py.File(name, 'r')
geometry = np.array(list(h5py.File(name, 'r')["Mesh"]["0"]["mesh"]["geometry"]))

Nx = 150
Ny = 100
x = np.linspace(0, Lx, Nx)
y = np.linspace(-1-epsilon, 1+epsilon, Ny)
X, Y = np.meshgrid(x,y)

frames = 100
half_period = int(tau/(dt))
skip = int(half_period/(2*frames))
skip = 4
s = 10


dat = 15
dat2 = 9
y1 = np.array([-1-epsilon, -1-epsilon+0.03, -1-epsilon+0.06, -1-epsilon+0.09])
dat = len(y1)
stream_points     = np.array(list(zip( 0.33*np.ones(dat),  y1)))
stream_points2    = np.array(list(zip( 0.33*np.ones(dat), -y1)))
saves = [6, 16, 43, 57, 58, 59]

fig = plt.figure(4)
for i in range(frames):
    plt.clf()
    u = np.array(list(h5py.File(name, 'r')["VisualisationVector"][str(int(N-skip*i))])) 
    # Interpolate uneven grid onto an even grid

    ux_grid  = np.roll(griddata((geometry[::s,0], geometry[::s,1]), u[::s,0], (X, Y), method='cubic'), 0, axis=1)
    uy_grid  = np.roll(griddata((geometry[::s,0], geometry[::s,1]), u[::s,1], (X, Y), method='cubic'), 0, axis=1)

    x_, y_ = np.meshgrid(kappa, epsilon)
    ax1 = plt.contourf(X,Y, np.sqrt(np.square(ux_grid) + np.square(uy_grid)), cmap=Map, levels=15)
    cbar = fig.colorbar(ax1, format='%1.0f')
    #plt.quiver(x[::s], y[::s], ux_grid[::s, ::s], uy_grid[::s, ::s])
    plt.streamplot(X, Y, ux_grid, uy_grid, color='k',  density=3.25, start_points=stream_points )
    plt.streamplot(X, Y, ux_grid, uy_grid, color='k',  density=3.25, start_points=stream_points2 )
    plt.streamplot(X, Y, ux_grid, uy_grid, color='k',  density=0.75)

    plt.fill_between(x,  1+epsilon*np.ones(len(x)), +1+epsilon*np.cos(kappa*x), color="k")
    plt.fill_between(x, -1-epsilon*np.ones(len(x)), -1-epsilon*np.cos(kappa*x), color="k")
    cbar.ax.set_ylabel(r'Velocity $u$', fontsize=8)
    plt.ylabel(r"Vertical position $y$ [$a$]",    fontsize=8)
    plt.xlabel(r"Horizontal position $x$ [$a$]",    fontsize=8)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.tick_params(axis='both', which='minor', labelsize=8)
    plt.axis([0, 6.18, -1.4, 1.4])
    #time = ( ( (t[N-skip*i]-t[N])*omega ) % 2*np.pi )/(2*np.pi)
    #plt.title(str(int(i)))
    #plt.title(r"$t=$ %2.2f $2\pi$ [$\omega$]" % time, fontsize=8)
    #plt.axis("equal")
    if i in saves:
        filename = root + "frame_"+str(i)+".pdf"
        plt.savefig(filename, bbox_inches="tight")
        os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
    plt.pause(0.01)
plt.show() 