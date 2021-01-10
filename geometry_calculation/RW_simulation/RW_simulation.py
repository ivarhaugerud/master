import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py
from scipy.interpolate import griddata

dirr = "flow_fields/Lx3.92_tau5.0_eps0.2_nu16.0_D1.0_fzero0.0_fone10.0_res100_dt0.01/"

#simulation paramters
dt = 0.01
tau = 5.0 
timesteps = int(tau/dt)
periods = 2000
datafiles = periods*100

#geometry parameters
epsilon = 0.2
kappas = np.array([1.6])
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
    filename = dirr
    tdat = np.loadtxt(dirr + "tdata.dat")
    names.append(filename)

    time = tdat[:,0]
    u2 = tdat[:,8]
    exp_u2[j] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)
    print("U = ", exp_u2[j])

for i in range(len(Lx)):
    kappa = kappas[i]
    l = Lx[i]
    name = dirr + "u.h5"
    f = h5py.File(name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))

    Nx = int(600)
    Ny = int(700)
    x = np.linspace(0, l, Nx)
    y = np.linspace(-1-epsilon, 1+epsilon, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    U = np.sqrt(exp_u2[i])
    U_scale = 1

    for j in range(timesteps):
        u = np.array(list(f["VisualisationVector"][str(periods_of_flow*timesteps-j)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), U_scale*u[:,0]/U, (X, Y), method='nearest')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), U_scale*u[:,1]/U, (X, Y), method='nearest')

        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)
        print(j, timesteps)

    Pe = 0.5
    D = U_scale/Pe
    alpha = np.sqrt(2*D*dt)

    for k in range(periods*timesteps):
        pos[0, :] = pos[0, :] + dt * interpolation["x-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False) + alpha*np.random.normal(loc=0, scale=1, size=N)
        pos[1, :] = pos[1, :] + dt * interpolation["y-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False) + alpha*np.random.normal(loc=0, scale=1, size=N)

        pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        if int(k) % int(periods*timesteps/datafiles) == 0:
            np.save('data/run_12_01/RW_positions_' +str(k), pos[:, :])