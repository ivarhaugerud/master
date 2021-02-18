import numpy as np
import scipy.interpolate as sci
#import h5py
from scipy.interpolate import griddata

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 1000

#geometry parameters
epsilon = 0.25
Lx = np.array([1.05, 2.09, 6.28, 9.42, 12.56, 15.71, 25.13])
kappas  = 2*np.pi/Lx

#flow parameters
periods_of_flow = 8/3

dirr = []
for i in range(len(Lx)):
    dirr.append("flow_fields/non_zero_eps/Lx"+str(Lx[i])+"_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

#for RW simulation 
N  = int(2000)    #number of random walkers
nu = 3.6
F0 = 12.0/nu
omega = 2*np.pi/tau 
Sc = nu
D  = 1
RW_steps = 30
alpha = np.sqrt(2*D*dt/RW_steps)
var = np.zeros((len(Lx), periods*timesteps))
t = np.linspace(0, periods*tau, periods*timesteps)

for i in range(len(Lx)):
    prev_pos = np.zeros((2, N))
    pos      = np.zeros((2, N))
    pos[1,:]      = np.random.uniform(-1-epsilon+0.01, 1+epsilon-0.01, N)
    prev_pos[1,:] = np.copy(pos[1,:])
    kappa = kappas[i]
    l = Lx[i]



    name = dirr[i] + "u.h5"
    geometry = np.array(list(h5py.File(name, 'r')["Mesh"]["0"]["mesh"]["geometry"]))

    Nx = int(450)
    Ny = int(700)
    x = np.linspace(0, l, Nx)
    y = np.linspace(-1-epsilon, 1+epsilon, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    for j in range(int(timesteps)):
        u = np.array(list(h5py.File(name, 'r')["VisualisationVector"][str(int(periods_of_flow*timesteps)-j-1)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='cubic')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='cubic')

        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)
        print(Lx[i], j, int(periods_of_flow*timesteps*skip)-j*skip-1)
        print(np.mean(np.mean(ux_grid)))

    for k in range(int(periods*timesteps)):
        plt.clf()
        pos[0, :] += dt * interpolation["x-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)
        pos[1, :] += dt * interpolation["y-"+str(int((k+timesteps)%timesteps))](pos[1, :], (pos[0, :]+l)%l, grid=False)

        pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        for j in range(RW_steps):
           pos[0, :] += alpha*np.random.normal(loc=0, scale=1, size=N)
           pos[1, :] += alpha*np.random.normal(loc=0, scale=1, size=N)
           pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
           pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
           prev_pos = np.copy(pos)

        plt.scatter(pos[0, :], pos[1, :])
        plt.pause(0.3)

        var[i, k] = np.var(pos[0, :])

    np.save(dirr[i]+"var_over_2Dmt", var[i, :]/(2*D*t))
    print("DONE WITH RUN FOR KAPPA: ", kappa)