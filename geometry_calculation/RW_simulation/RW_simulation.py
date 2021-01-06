import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py

dirr = "../results_oscwavychannel/run_12_11/"

#simulation paramters
dt = 0.01
tau = 5.0 

epsilon = 0.2
kappas = np.array([0.1,  0.7, 1.5])
Lx = 2*np.pi/kappas

omega = 2*np.pi/tau
nu = 1.2
D = 1
f1 = 3
F0 = f1/nu
Sc = nu
gamma = np.sqrt(1j*omega/Sc)
timesteps = int(tau/dt)

exp_u2 = np.zeros(len(kappas))
periods = 2
names = []

for j in range(len(kappas)):
    res = int(100*(1+float(epsilon)))
    try:
        filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(epsilon)[:4])+"_nu1.2_D0.3_fzero0.0_fone3.0_res"+str(res)+"_dt0.01/"
        tdat = np.loadtxt(dirr + filename + "tdata.dat")
        names.append(filename)

    except:
        filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(epsilon)[:3])+"_nu1.2_D0.3_fzero0.0_fone3.0_res"+str(res)+"_dt0.01/"
        tdat = np.loadtxt(dirr + filename + "tdata.dat")
        names.append(filename)

    t = tdat[:,0]
    u2 = tdat[:,4]
    exp_u2[j] = integrate.trapz(u2[-timesteps:], t[-timesteps:])/(tau)

for i in range(len(Lx)):
    name = names[i] + "u.h5"
    print(dirr + name)

    f = h5py.File(dirr + name, 'r')

    #index of all lattice-points, 0 could be any timestep as this is indep of time
    geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))
    x, y = geometry[:,0], geometry[:,1]
    speed = np.zeros(timesteps)

    Nx = int(2.5*1e3)
    Ny = int(2.5*1e3)
    x = np.linspace(min(geometry[:,0]), max(geometry[:,0]), Nx)
    y = np.linspace(min(geometry[:,1]), max(geometry[:,1]), Ny)
    X, Y = np.meshgrid(x,y)
    ux = np.zeros((len(x), len(y), timesteps))
    uy = np.zeros((len(x), len(y), timesteps))

    for i in range(timesteps):
        # read data
        u = np.array(list(f["VisualisationVector"][str(-i)]))        
        
        # Interpolate uneven grid onto an even grid
        ux[:,:,i] = griddata((geometry[:,0], geometry[:,1]), u[:,0], (X, Y), method='nearest')
        uy[:,:,i] = griddata((geometry[:,0], geometry[:,1]), u[:,1], (X, Y), method='nearest')

    u_x = RegularGridInterpolator((x,y,t), ux)
    u_y = RegularGridInterpolator((x,y,t), uy)



#to check if flow-field is correct
#plt.quiver(x_axis[:int(b/2)], y_axis[:int(a/2)], frame[:int(a/2),:int(b/2),1], frame[:int(a/2),:int(b/2), 0])
#plt.show()



N  = int(1e3)    #number of random walkers
prev_pos = np.zeros((2, N))
pos = np.zeros((2, N))

pos[1,:] = np.random.uniform(-1+epsilon, 1-epsilon, N)

U = np.sqrt(exp_u2[j])
Pe = 5
D = a*U/Pe
alpha = np.sqrt(2*D*dt)
period  = 2*np.pi/omega
periods = 100
t = np.linspace(0, periods*period, int(periods*timesteps))
l = Lx[i]
datafiles = 1000

for i in range(timesteps):
    pos[0, :] = pos[0, :] + dt * u_x( (pos[0, :]+l)%l, pos[1, :], (t[i]+perio)%period)  + alpha*np.random.normal(loc=0, scale=1, size=N)
    pos[1, :] = pos[1, :] + dt * u_y( (pos[0, :]+l)%l, pos[1, :], (t[i]+perio)%period)  + alpha*np.random.normal(loc=0, scale=1, size=N)

    pos[1, np.where( abs(pos[1, :]) > 1+epsilon*np.cos(pos[0,:]))] = prev_pos[1, np.where( abs(pos[1, :]) > 1+epsilon*np.cos(pos[0,:]))] #checks if y-coordinate outside

    #update previous position
    prev_pos = pos

    if int(i) % int(len(t)/datafilesa) == 0:
        np.save('data/RW_positions_' +str(i), pos[:, :])