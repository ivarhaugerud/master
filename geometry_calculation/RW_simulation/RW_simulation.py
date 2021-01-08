import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py
from scipy.interpolate import griddata

dirr = "flow_fields/"

#simulation paramters
dt = 0.01
tau = 5.0 
timesteps = int(tau/dt)
print(timesteps)
periods = 100
datafiles = 5000

#geometry parameters
epsilon = 0.2
kappas = np.array([0.4])
Lx = 2*np.pi/kappas

#flow parameters
omega = 2*np.pi/tau
nu = 16
D = 1
f1 = 10
F0 = f1/nu
Sc = nu
periods_of_flow = 3

exp_u2 = np.zeros(len(kappas))
names = []


#for RW simulation 
N  = int(1e3)    #number of random walkers
prev_pos = np.zeros((2, N))
pos      = np.zeros((2, N))
pos[1,:]      = np.random.uniform(-1+epsilon, 1-epsilon, N)
prev_pos[1,:] = np.copy(pos[1,:])

for j in range(len(kappas)):
    res = int(100)
    filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(epsilon)[:4])+"_nu16.0_D1.0_fzero0.0_fone10.0_res"+str(res)+"_dt0.01/"
    tdat = np.loadtxt(dirr + filename + "tdata.dat")
    names.append(filename)

    time = tdat[:,0]
    u2 = tdat[:,8]
    exp_u2[j] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)

for i in range(len(Lx)):
    kappa = kappas[i]
    l = Lx[i]
    name = names[i] + "u.h5"
    f = h5py.File(dirr + name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))

    Nx = int(600)
    Ny = int(700)
    x = np.linspace(0, l, Nx)
    y = np.linspace(-1-epsilon, 1+epsilon, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    U = np.sqrt(exp_u2[i])
    U_scale = 2*l

    for j in range(timesteps):
        u = np.array(list(f["VisualisationVector"][str(periods_of_flow*timesteps-j)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), U_scale*u[:,0]/U, (X, Y), method='nearest')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), U_scale*u[:,1]/U, (X, Y), method='nearest')

        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)
        print(j, timesteps)

    Pe = 10
    D = U_scale/Pe
    alpha = np.sqrt(2*D*dt)

    for k in range(periods*timesteps):
        pos[0, :] = pos[0, :] + dt * interpolation["x-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False) + alpha*np.random.normal(loc=0, scale=1, size=N)
        pos[1, :] = pos[1, :] + dt * interpolation["y-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False) + alpha*np.random.normal(loc=0, scale=1, size=N)

        pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        if int(k) % int(periods*timesteps/datafiles) == 0:
            np.save('data/run_09_01/RW_positions_' +str(k), pos[:, :])