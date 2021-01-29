import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py
from scipy.interpolate import griddata

#simulation parameters
dt           = 0.006
tau          = 3.0 
timesteps    = int(tau/(dt))
periods      = 10000
datafiles    = periods*20
N            = int(1e3)
RW_timesteps = 25
Pe           = 1

#geometry parameters
epsilon = 0.0
visc = np.array([1.5, 3.0, 5.0])
Lx = 12.56

#flow parameters
periods_of_flow = 8/3
exp_u2 = np.zeros(len(visc))
exp_D  = np.zeros(len(visc))

dirr = []

for i in range(len(visc)):
    dirr.append("flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

#for RW simulation 
print(np.random.normal(loc=0, scale=1, size=(2, N)))


for j in range(len(visc)):
    tdat = np.loadtxt(dirr[j] +"tdata.dat")

    time = tdat[:,0]
    u2   = tdat[:,4]
    exp_u2[j] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)
    exp_D[j] = integrate.trapz(tdat[-timesteps:, 8], time[-timesteps:])/(tau)

plt.plot(visc, exp_u2)
plt.show()

for i in range(len(visc)):
    prev_pos = np.zeros((2, N))
    pos      = np.zeros((2, N))
    pos[1,:]      = np.random.uniform(-1+epsilon, 1-epsilon, N)
    prev_pos[1,:] = np.copy(pos[1,:])
    l = Lx
    name = dirr[i] + "u.h5"
    f = h5py.File(name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))

    Nx = int(600)
    Ny = int(700)
    x = np.linspace(0,  l, Nx)
    y = np.linspace(-1, 1, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    for j in range(int(timesteps)):
        u = np.array(list(f["VisualisationVector"][str(int(periods_of_flow*timesteps)-j-1)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='nearest')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='nearest')

        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)
        print(visc[i], j, str(int(periods_of_flow*timesteps)-j-1))

    U = np.sqrt(exp_u2[i])
    D  = U/Pe
    alpha = np.sqrt(2*D*dt/RW_timesteps)

    for k in range(int(periods*timesteps)):
        pos[0, :] += dt * interpolation["x-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)
        prev_pos = np.copy(pos)

        for j in range(RW_timesteps):
            pos[:, :] += alpha*np.random.normal(loc=0, scale=1, size=(2, N))
            pos[:, np.where( abs(pos[1, :]) >   1)] = prev_pos[:, np.where( abs(pos[1, :]) >  1)] #checks if y-coordinate outside
            prev_pos = np.copy(pos)
    
        if int(k) % int(periods*timesteps/(datafiles)) == 0:
            np.save(dirr[i]+"pos/RW_positions_" +str(k), pos[:, :])

    print("DONE WITH RUN FOR NU: ", visc[i])