import numpy as np
import scipy.interpolate as sci
import h5py
from scipy.interpolate import griddata
import sys 
import matplotlib.pyplot as plt 

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 3000
#skip = 10

#geometry parameters
#epsilon = float(sys.argv[1])
l   = float(sys.argv[1])
eps = np.array([0.25])
#kappas = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1])
#Lx = 2*np.pi/kappas
#Lx = np.array([10.47, 2.991, 6.283, 4.487, 3.695])
kappa = 2*np.pi/l
#flow parameters
periods_of_flow = 8/3

dirr = []
for i in range(len(eps)):
    dirr.append("RW_simulation/flow_fields/non_zero_eps/Lx"+str(l)+"_tau3.0_eps"+str(eps[i])+"_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/")
    print(l, eps[i])
    a = np.loadtxt(dirr[i] + "tdata.dat")

#for RW simulation 
N  = int(500)    #number of random walkers
nu = 1.2
F0 = 12.0/nu
omega = 2*np.pi/tau 
Sc = nu
D  = 0.01
RW_steps = 12 # -> alpha = 0.01
alpha = np.sqrt(2*D*dt/RW_steps)
var = np.zeros((len(eps), periods*timesteps))
t   = np.linspace(0, periods*tau, periods*timesteps)
xi = np.linspace(-1, 1, 1000)
ETA = np.linspace(-3*l , 3*l, 500)

for i in range(len(eps)):
    epsilon = eps[i]
    prev_pos = np.zeros((2, N))
    pos      = np.zeros((2, N))
    pos[1,:]      = np.random.uniform(-1+epsilon+0.01, 1-epsilon-0.01, N)
    #prev_pos[1,:] = np.copy(pos[1,:])
    #kappa = kappas[i]
    #l = Lx[i]


    name = dirr[i] + "u.h5"
    geometry = np.array(list(h5py.File(name, 'r')["Mesh"]["0"]["mesh"]["geometry"]))

    Nx = int(600)
    Ny = int(700)
    x = np.linspace(0, l, Nx)
    y = np.linspace(-1-epsilon, 1+epsilon, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}

    for j in range(int(timesteps/skip)):
        u = np.array(list(h5py.File(name, 'r')["VisualisationVector"][str(int(periods_of_flow*timesteps)-j*skip-1)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='cubic')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='cubic')

        #print(np.shape(ux_grid))
        """
        fig,ax=plt.subplots(1,1)
        cp = ax.contourf(X, Y, ux_grid)
        fig.colorbar(cp) # Add a colorbar to a plot
        plt.show()

        fig,ax=plt.subplots(1,1)
        cp = ax.contourf(X, Y, uy_grid)
        fig.colorbar(cp) # Add a colorbar to a plot
        plt.show()
        """
        #fig,ax=plt.subplots(1,1)
        #cp = ax.contourf(X, Y, np.abs(ux_grid), levels=20)
        #fig.colorbar(cp) # Add a colorbar to a plot
        #plt.plot(ETA, 1+epsilon*np.cos(kappa*ETA), "k")
        #plt.plot(ETA,-1-epsilon*np.cos(kappa*ETA), "k")
        #plt.show()
        """
        for s in range(len(x)):
            for t in range(len(y)):
                if y[t] > 1+eps*np.cos(kappa*x[s]) or y[t] < -1-eps*np.cos(kappa*x[s]):
                    ux_grid[t, s] = 0
                    uy_grid[t, s] = 0
        """
        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)
    
    print("DONE WITH FINDING U: ", l, epsilon, np.mean(interpolation["x-"+str(j)](xi, 0)) )
        #fig,ax=plt.subplots(1,2)
        #cp = ax.contourf(X, Y, np.abs(ux_grid), levels=20)
        #fig.colorbar(cp) # Add a colorbar to a plot
        #plt.plot(ETA, 1+epsilon*np.cos(kappa*ETA), "k")
        #plt.plot(ETA,-1-epsilon*np.cos(kappa*ETA), "k")
        #plt.show()

    for k in range(int(periods*timesteps)):
        pos[0, :] += skip*dt * interpolation["x-"+str(int((k+timesteps/skip)%timesteps/skip))](pos[1, :], (pos[0, :]+l)%l, grid=False)
        pos[1, :] += skip*dt * interpolation["y-"+str(int((k+timesteps/skip)%timesteps/skip))](pos[1, :], (pos[0, :]+l)%l, grid=False)


        ux_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='cubic')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='cubic')
        #print(np.shape(ux_grid))
        plt.clf()

        plt.scatter(pos[0, :], pos[1, :])
        plt.plot(ETA, 1+epsilon*np.cos(kappa*ETA), "k")
        plt.plot(ETA,-1-epsilon*np.cos(kappa*ETA), "k")
        plt.pause(0.01)

        pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        for j in range(RW_steps):
           pos[0, :] += alpha*np.random.normal(loc=0, scale=1, size=N)
           pos[1, :] += alpha*np.random.normal(loc=0, scale=1, size=N)
           pos[:, np.where( pos[1, :] >   1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
           pos[:, np.where( pos[1, :] <  -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
           prev_pos = np.copy(pos)

        var[i, k] = np.var(pos[0, :])

    np.save(dirr[i]+"var_over_2Dmt_D01", var[i, :]/(2*D*t))
    print("DONE WITH RUN FOR KAPPA: ", kappa)

