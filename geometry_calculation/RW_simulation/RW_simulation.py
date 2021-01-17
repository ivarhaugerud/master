import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py
from scipy.interpolate import griddata

#dirr = "flow_fields/Lx62.8_tau5.0_eps0.2_nu16.0_D1.0_fzero0.0_fone10.0_res100_dt0.01/"
#Lx1.0_tau5.0_eps0.0_nu0.5_D1.0_fzero0.0_fone3.0_res100_dt0.01
dirr = "flow_fields/zero_eps/nu_0.5/"
#simulation paramters
dt = 0.01
tau = 5.0 
timesteps = int(tau/dt)
periods   = 100
datafiles = periods*25

#geometry parameters
epsilon = 0.0
kappas  = np.array([1e-5])
Lx      = 2*np.pi/kappas

#flow parameters
omega = 2*np.pi/tau
nu = 0.5
D  = 1.0/3.0
F0 = 3
Sc = 1/nu
periods_of_flow = 3

exp_u2 = np.zeros(len(kappas))


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

    Pe = 3.0
    D  = U_scale/Pe
    alpha = np.sqrt(2*D*dt)

    for k in range(2*periods*timesteps):
        if k % 2 == 0:
            #print("even:", str(int((k/2+timesteps)%timesteps)))
            pos[0, :] += alpha*np.random.normal(loc=0, scale=1, size=N) + 0.5* dt * interpolation["x-"+str(int((k/2+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)
            pos[1, :] += alpha*np.random.normal(loc=0, scale=1, size=N) + 0.5* dt * interpolation["y-"+str(int((k/2+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)

        else:
            #print("odd:", str(int(((k-1)/2+timesteps)%timesteps)), str(int(((k+1)/2+timesteps)%timesteps)))
            pos[0, :] += alpha*np.random.normal(loc=0, scale=1, size=N) + 0.5* (dt * interpolation["x-"+str(int(((k+1)/2+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False) + dt * interpolation["x-"+str(int(((k-1)/2+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False))
            pos[1, :] += alpha*np.random.normal(loc=0, scale=1, size=N) + 0.5* (dt * interpolation["y-"+str(int(((k+1)/2+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False) + dt * interpolation["y-"+str(int(((k-1)/2+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False))

        pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        if int(k) % int(periods*timesteps/datafiles) == 0:
            #np.save('data/Lx62_8/RW_positions_' +str(k), pos[:, :])
            np.save(dirr+'pos/RW_positions_' +str(k), pos[:, :])