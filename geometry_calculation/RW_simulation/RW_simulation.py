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
timesteps = tau/dt

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
"""
x = len(velocity[0,:,0])
y = len(velocity[:,0,0])

a  = int(y/2)
b  = int(x/2)
b2 = int(b/2)
l  = x
print(a, b, b2, l)

frame = np.zeros((y+2, x+2, 2))

#so that we have zero around the edges
frame[1:-1, 1:-1, :] = velocity[:,:,:]
#so that it is periodic in x
frame[1:-1, 0, :]    = velocity[:, -1, :]
frame[1:-1, -1, :]   = velocity[:, 0, :]

x_axis = np.linspace(-b, b, x+2)
y_axis = np.linspace(-a, a, y+2)

#to check if flow-field is correct
#plt.quiver(x_axis[:int(b/2)], y_axis[:int(a/2)], frame[:int(a/2),:int(b/2),1], frame[:int(a/2),:int(b/2), 0])
#plt.show()

#interpoalte function using splines in 2D
u_x = sci.RectBivariateSpline(y_axis, x_axis, frame[:,:,1])
u_y = sci.RectBivariateSpline(y_axis, x_axis, frame[:,:,0])
#argument is (y,x)


N  = int(1e3)    #number of random walkers
Nt = int(1e4)  #number of timesteps
T  = int(1e6)

time = np.linspace(0, T, Nt)
dt = T/Nt

prev_pos = np.zeros((2, N))
pos = np.zeros((2, N))

pos[1,:] = np.random.uniform(-a, a, N)

U = (np.mean(np.mean(frame[:, :b2, 1])) + np.mean(np.mean(frame[:, -b2:, 1])) + np.mean(np.mean(frame[b:-b, b2:-b2, 1])))/3
print(a*U)
Pe = 1
D = a*U/Pe
alpha = np.sqrt(2*D*dt)

for t in range(Nt-1):
    pos[0, :] = pos[0, :] + dt * u_x(pos[1, :], pos[0, :]%l, grid=False)  + alpha*np.random.normal(loc=0, scale=1, size=N)
    pos[1, :] = pos[1, :] + dt * u_y(pos[1, :], pos[0, :]%l, grid=False)  + alpha*np.random.normal(loc=0, scale=1, size=N)

    pos[1, np.where( abs(pos[1, :]) > a)] = prev_pos[1, np.where( abs(pos[1, :]) > a)] #checks if y-coordinate outside

    #check if inside box
    pos[:, np.where(  np.invert(((abs(pos[0, :])%l<b2) |  (abs(pos[0, :])%l>b2+b))) & (np.abs(pos[1, :]) > a - b) )] = prev_pos[:, np.where( np.invert(((abs(pos[0, :])%l<b2) |  (abs(pos[0, :])%l>b2+b))) & (np.abs(pos[1, :]) > a - b) )]

    #update previous position
    prev_pos[:, :] = pos[:, :]

    if int(t+2) % int(1000) == 0:
        np.save('../data/RW_rough_' +str(t+2), pos[:, :])
        print(t+2, "Here")
"""