import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py
from scipy.interpolate import griddata

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 5000
datafiles = periods*20

#geometry parameters
epsilon = 0.25
Lx = np.array([12.56, 15.71, 25.13])
Lx =  np.array([1.05, 2.09, 6.28, 9.42])
kappas  = 2*np.pi/Lx

#flow parameters
periods_of_flow = 8/3
exp_u2 = np.zeros(len(kappas))

dirr = []

for i in range(len(Lx)):
    dirr.append("flow_fields/non_zero_eps/Lx"+str(Lx[i])+"_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

#for RW simulation 
N  = int(1e3)    #number of random walkers

xi = np.linspace(-1, 1, int(1e5))
nu = 3.6
F0 = 12.0/nu
omega = 2*np.pi/tau 
Sc = nu
gamma = np.sqrt(1j*omega/Sc)
ux = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)/2
U_ana = integrate.trapz(ux*np.conj(ux), xi)


for j in range(len(kappas)):
    tdat = np.loadtxt(dirr[j] +"tdata.dat")

    time = tdat[:,0]
    u2   = tdat[:,4]
    exp_u2[j] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)

for i in range(len(Lx)):
    prev_pos = np.zeros((2, N))
    pos      = np.zeros((2, N))
    pos[1,:]      = np.random.uniform(-1+epsilon, 1-epsilon, N)
    prev_pos[1,:] = np.copy(pos[1,:])
    kappa = kappas[i]
    l = Lx[i]
    name = dirr[i] + "u.h5"
    f = h5py.File(name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))

    Nx = int(600)
    Ny = int(700)
    x = np.linspace(0, l, Nx)
    y = np.linspace(-1-epsilon, 1+epsilon, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    for j in range(int(timesteps)):
        u = np.array(list(f["VisualisationVector"][str(int(periods_of_flow*timesteps)-j-1)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='nearest')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='nearest')

        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)
        print(Lx[i], j, str(int(periods_of_flow*timesteps)-j-1))

    U = np.sqrt(exp_u2[i])
    Pe = 6
    D  = U/Pe
    alpha = np.sqrt(2*D*dt)

    for k in range(int(periods*timesteps)):
        pos[0, :] += alpha*np.random.normal(loc=0, scale=1, size=N) + dt * interpolation["x-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)
        pos[1, :] += alpha*np.random.normal(loc=0, scale=1, size=N) + dt * interpolation["y-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)

        pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        if int(k) % int(periods*timesteps/(datafiles)) == 0:
            np.save(dirr[i]+"pos/RW_positions_" +str(k), pos[:, :])

    print("DONE WITH RUN FOR KAPPA: ", kappa)