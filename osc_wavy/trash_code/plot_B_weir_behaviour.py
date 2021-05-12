import numpy as np
import scipy.interpolate as sci
import h5py
from scipy.interpolate import griddata
import matplotlib.pyplot as plt 
import matplotlib
Map = matplotlib.cm.get_cmap('Spectral_r')

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'  
root = "../../../master_latex/results/"
base = "../data_test/vary_geometry/"

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 1000

#geometry parameters
epsilon = 0.3
Lx = 7.39
nu = 1.2
D = 1.0
F0 = 12/nu
kappa  = 2*np.pi/Lx
omega = 2*np.pi/tau
dt = 0.004
tau = 3.0

print("Kappa = ", kappa)

dirr = "data_test/write_u/Lx7.39_tau3.0_eps0.3_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt0.004/"
tdata = np.loadtxt(dirr+"tdata.dat")


T = np.max(tdata[:,0])
N = int(T/dt)-1

t = tdata[:,0]
#plt.plot(tdata[:,0], tdata[:,8])
#plt.show()

name = dirr + "B.h5"
data = h5py.File(name, 'r')
B = (data["VisualisationVector"][str(N)])[:,0]
geometry = np.array(list(h5py.File(name, 'r')["Mesh"]["0"]["mesh"]["geometry"]))

Nx = 300
Ny = 100
x = np.linspace(0, Lx, Nx)
y = np.linspace(-1-epsilon, 1+epsilon, Ny)
X, Y = np.meshgrid(x,y)

frames = 10
half_period = int(tau/(dt))
skip = int(half_period/frames)
s = 5
for i in range(frames):
    plt.clf()
    B = (data["VisualisationVector"][str(N-skip*i)])[:,0]
    B_grid  = griddata((geometry[:,0], geometry[:,1]), B, (X, Y), method='cubic')

    grad = np.gradient(B_grid, y, x)
    fig = plt.figure(4)
    x_, y_ = np.meshgrid(kappa, epsilon)
    ax1 = plt.contourf(X,Y, B_grid, cmap=Map, levels=np.linspace(-0.6, 0.65, 75))
    cbar = fig.colorbar(ax1, format='%1.2f')
    #plt.quiver(x[::s], y[::s], grad[0][::s, ::s], grad[1][::s, ::s])
    plt.fill_between(x,  1+epsilon*np.ones(len(x)),  1+epsilon*np.cos(kappa*x), color="k")
    plt.fill_between(x, -1-epsilon*np.ones(len(x)), -1-epsilon*np.cos(kappa*x), color="k")
    cbar.ax.set_ylabel(r'Brenner field $B$', fontsize=8)
    plt.ylabel(r"Vertical position $y$ [$a$]",    fontsize=8)
    plt.xlabel(r"Horizontal position $x$ [$a$]",    fontsize=8)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.tick_params(axis='both', which='minor', labelsize=8)
    time = ( ( (t[N-skip*i]-t[N])*omega ) % 2*np.pi )/(2*np.pi)
    plt.title(r"$t=$ %2.2f $2\pi$ [$\omega$]" % time, fontsize=8)
    plt.axis("equal")
    plt.pause(0.01)
plt.show() 