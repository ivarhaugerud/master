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
    result = np.load(dirr[i]+"var_over_t.npy")
    plt.plot(result)
plt.show()

periods      = 1000
N            = int(2000)
RW_timesteps = 30
Pe           = 1
var = np.zeros((len(visc), int(periods*timesteps)))
t   = np.linspace(0, periods*tau, int(periods*timesteps))

for i in range(len(visc)):
    prev_pos = np.zeros((2, N))
    pos      = np.zeros((2, N))
    pos[1,:]      = np.random.uniform(-0.99, 0.99, N)
    prev_pos[1,:] = np.copy(pos[1,:])
    l = Lx
    name = dirr[i] + "u.h5"
    f = h5py.File(name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))

    Nx = int(10)
    Ny = int(700)
    x = np.linspace(0,  l, Nx)
    y = np.linspace(-1, 1, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    for j in range(int(timesteps)):
        u = np.array(list(f["VisualisationVector"][str(int(periods_of_flow*timesteps)-j-1)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='cubic')
        plt.plot(ux_grid)
        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        print(visc[i], j, str(int(periods_of_flow*timesteps)-j-1))

    D  = 0.2
    alpha = np.sqrt(2*D*dt/RW_timesteps)

    for k in range(int(periods*timesteps)):
        pos[0, :] += dt * interpolation["x-"+str((k+timesteps)%timesteps)](pos[1, :], (pos[0, :]+l)%l, grid=False)
        prev_pos = np.copy(pos)

        for j in range(RW_timesteps):
            pos[:, :] += alpha*np.random.normal(loc=0, scale=1, size=(2, N))
            pos[:, np.where( abs(pos[1, :]) > 1)] = prev_pos[:, np.where( abs(pos[1, :]) >  1)]
            prev_pos = np.copy(pos)

        var[i, k] = np.var(pos[0, :])

    np.save(dirr[i]+"var_over_t", var[i, 1:]/(2*D*t[1:]))
    print("DONE WITH RUN FOR NU: ", visc[i], np.mean(var[i, int(periods*timesteps/2):]/(2*D*t[int(periods*timesteps/2):])))