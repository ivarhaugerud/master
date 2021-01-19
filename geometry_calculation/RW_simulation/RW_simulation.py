import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py
from scipy.interpolate import griddata

dirr = "flow_fields/zero_eps/mu_0.5/"

dt = 0.001
tau = 3.0 
skip = 2
timesteps = int(tau/(dt))
periods   = 5000
datafiles = periods*25

#geometry parameters
epsilon = 0.0
kappas  = np.array([1e-5])
Lx      = 2*np.pi/kappas

#flow parameters
periods_of_flow = 3

exp_u2 = np.zeros(len(kappas))


#for RW simulation 
N  = int(1e3)    #number of random walkers
prev_pos = np.zeros((2, N))
pos      = np.zeros((2, N))
pos[1,:]      = np.random.uniform(-1+epsilon, 1-epsilon, N)
prev_pos[1,:] = np.copy(pos[1,:])

xi = np.linspace(-1, 1, int(1e5))
nu = 0.5
F0 = 3.0/nu
omega = 2*np.pi/tau 
Sc = nu
gamma = np.sqrt(1j*omega/Sc)
ux = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)/2
U_ana = integrate.trapz(ux*np.conj(ux), xi)


for j in range(len(kappas)):
    res = int(100)
    filename = dirr
    tdat = np.loadtxt(dirr + "tdata.dat")

    time = tdat[:,0]
    u2   = tdat[:,4]
    exp_u2[j] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)
    #print(np.mean(u2[-timesteps:]))
    plt.plot(time, u2)
    plt.plot(time[-timesteps:], u2[-timesteps:])
    print("U = ", (exp_u2[j]))
    print("U_ana = ", (U_ana))
    plt.show()
"""
timesteps = int(timesteps/skip)

for i in range(len(Lx)):
    kappa = kappas[i]
    l = Lx[i]
    name = dirr + "u.h5"
    print("before reading")
    f = h5py.File(name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))
    print("after saving geometry")
    Nx = int(60)
    Ny = int(700)
    x = np.linspace(0, l, Nx)
    y = np.linspace(-1-epsilon, 1+epsilon, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    U = np.sqrt(exp_u2[i])
    U_scale = 5

    for j in range(timesteps):
        u = np.array(list(f["VisualisationVector"][str(periods_of_flow*skip*timesteps-skip*j-1)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), U_scale*u[:,0]/U, (X, Y), method='nearest')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), U_scale*u[:,1]/U, (X, Y), method='nearest')

        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)
        print(j, str(periods_of_flow*skip*timesteps-skip*j-1), timesteps, 2*timesteps)

    Pe = 20
    dt = dt*skip
    D  = U_scale/Pe
    alpha = np.sqrt(2*D*dt)

    for k in range(periods*timesteps):
        pos[0, :] += alpha*np.random.normal(loc=0, scale=1, size=N) + dt * interpolation["x-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)
        pos[1, :] += alpha*np.random.normal(loc=0, scale=1, size=N)# + dt * interpolation["y-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)

        pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        if int(k) % int(periods*timesteps/datafiles) == 0:
            #np.save('data/Lx62_8/RW_positions_' +str(k), pos[:, :])
            np.save(dirr+'pos_2/RW_positions_' +str(k), pos[:, :])
"""