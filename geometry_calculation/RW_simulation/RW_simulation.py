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
dt = 1.0
tau = 5.0 

epsilon = 0.2
kappas = np.array([0.4])
Lx = 2*np.pi/kappas

omega = 2*np.pi/tau
nu = 1.2
D = 1
f1 = 3
F0 = f1/nu
Sc = nu
gamma = np.sqrt(1j*omega/Sc)
timesteps = int(tau/dt)
period  = 2*np.pi/omega
print(period, timesteps)
t = np.linspace(0, period, timesteps)

exp_u2 = np.zeros(len(kappas))
periods = 3
names = []


#for RW simulation 
N  = int(1e3)    #number of random walkers
prev_pos = np.zeros((2, N))
pos      = np.zeros((2, N))

pos[1,:]      = np.random.uniform(-1+epsilon, 1-epsilon, N)
prev_pos[1,:] = np.random.uniform(-1+epsilon, 1-epsilon, N)

periods = 100
datafiles = 5000




for j in range(len(kappas)):
    res = int(100)
    try:
        filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(epsilon)[:4])+"_nu16.0_D1.0_fzero0.0_fone10.0_res"+str(res)+"_dt0.01/"
        print(filename)
        tdat = np.loadtxt(dirr + filename + "tdata.dat")
        names.append(filename)

    except:
        filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(epsilon)[:3])+"_nu16.0_D1.0_fzero0.0_fone10.0_res"+str(res)+"_dt0.01/"
        tdat = np.loadtxt(dirr + filename + "tdata.dat")
        names.append(filename)

    time = tdat[:,0]
    u2 = tdat[:,8]
    exp_u2[j] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)
    print(exp_u2[j])

for i in range(len(Lx)):
    kappa = kappas[i]
    l = Lx[i]
    name = names[i] + "u.h5"
    print(dirr + name)

    f = h5py.File(dirr + name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))
    x, y = geometry[:,0], geometry[:,1]
    speed = np.zeros(timesteps)

    Nx = int(800)
    Ny = int(900)
    x = np.linspace(0, l, Nx)
    y = np.linspace(-1-epsilon, 1+epsilon, Ny)
    X, Y = np.meshgrid(x,y)

    interpolation = {}
    U = np.sqrt(exp_u2[i])
    a = 1
    Pe = 5
    D = a*U/Pe
    alpha = np.sqrt(2*D*dt)

    for j in range(timesteps):
        u = np.array(list(f["VisualisationVector"][str(1500-j)])) 
        # Interpolate uneven grid onto an even grid

        ux_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='nearest')
        uy_grid  = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='nearest')

        interpolation["x-"+str(j)]  = sci.RectBivariateSpline(y, x, ux_grid)
        interpolation["y-"+str(j)]  = sci.RectBivariateSpline(y, x, uy_grid)

    print(np.shape(interpolation["x-2"](y, x)))
    plt.imshow(interpolation["x-2"](y, x))
    plt.show()
    plt.plot(y, interpolation["x-2"](y, x[0]  )[:,0])  
    plt.plot(y, interpolation["x-2"](y, x[200])[:,0])    
    plt.plot(y, interpolation["x-2"](y, x[400])[:,0])    
    plt.plot(y, interpolation["x-2"](y, x[600])[:,0])    
    plt.show()

    plt.plot(y, interpolation["y-1"](0, y)[0])  
    plt.plot(y, interpolation["y-1"](0.33*np.pi/kappa, y)[0])    
    plt.plot(y, interpolation["y-1"](0.66*np.pi/kappa, y)[0])    
    plt.plot(y, interpolation["y-1"](0.99*np.pi/kappa, y)[0])    
    plt.show()


    for k in range(periods*timesteps):
        for n in range(N):
            pos[0, n] = pos[0, n] + dt * u_x( [(pos[0, n]+l)%l, pos[1, n], (t[int((k+timesteps)%timesteps)]+period)%period])
            pos[1, n] = pos[1, n] + dt * u_y( [(pos[0, n]+l)%l, pos[1, n], (t[int((k+timesteps)%timesteps)]+period)%period])

        pos[0, :] += alpha*np.random.normal(loc=0, scale=1, size=N)
        pos[1, :] += alpha*np.random.normal(loc=0, scale=1, size=N)
        #print("\n before" , pos[1, (np.where( pos[1, :] < -1-epsilon*np.cos(pos[0,:])))], prev_pos[1, (np.where( pos[1, :] < -1-epsilon*np.cos(pos[0,:])))])
        pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.cos(kappa*pos[0,:]))] #checks if y-coordinate outside
        #update previous position
        #print("\n after" , pos[1, (np.where( pos[1, :] < -1-epsilon*np.cos(pos[0,:])))])
        #print("\n\n\n")
        prev_pos = np.copy(pos)

        if int(k) % int(periods*timesteps/datafiles) == 0:
            np.save('data/run_08_01/RW_positions_' +str(k), pos[:, :])