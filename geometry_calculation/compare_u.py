import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata
import matplotlib.ticker as tick 

name1 = "results_oscwavychannel/Lx0.05_tau5.0_eps0.0_nu100.0_D0.3_fzero0.0_fone1.0_res100_dt0.01/u.h5" #results_oscwavychannel/Lx452.19_tau5.0_eps0.032_nu1.0_D0.3_fzero0.0_fone3.0_res106_dt0.02/u.h5"
name2 = "results_oscwavychannel/Lx0.05_tau2.5_eps0.0_nu100.0_D0.3_fzero0.0_fone1.0_res100_dt0.01/u.h5" #results_oscwavychannel/Lx452.19_tau5.0_eps0.032_nu1.0_D0.3_fzero0.0_fone3.0_res106_dt0.02/u.h5"

f = h5py.File(name1, 'r')
g = h5py.File(name2, 'r')

dt1 = 0.02 
dt2 = 0.01
tau1 = 5
tau2 = 2.5

skipp = 2
timesteps = int((5/0.02)/skipp)

speed = np.zeros((timesteps+1, 2))


max_index = len(list(f["VisualisationVector"]))
timesteps = 5
F0 = 1
omega = 1e-1#2*np.pi/tau2
nu = 100
D = 0.3
Sc = nu/D
gamma = -np.sqrt(1j*omega/Sc)
T = np.linspace(0, 3*2*np.pi/omega, max_index)

for i in range(max_index-1):
	plt.clf()
	u = np.array(list(f["VisualisationVector"][str(int((max_index-skipp*(i*0.5)-1)))]))
	#speed[i, 0] = np.mean(np.sqrt(np.square(u[:,0]) + np.square(u[:,1])))
	plt.plot(u[::11,0])

	u1 = np.array(list(g["VisualisationVector"][str((int(max_index/2)-skipp*i-1))]))
	speed[i, 1] = np.mean(np.sqrt(np.square(u1[:,0]) + np.square(u1[:,1])))
	plt.plot(u1[::11,0])
	x = np.linspace(-1, 1, len(u1[::11,0]))
	#plt.axis([0, len(u1[::11,0]), -0.4, 0.4])

	plt.plot(np.real(F0*((1-np.cosh(gamma*x)/np.cosh(gamma))/(gamma*gamma))*np.exp(1j*omega*T[max_index-skipp*i-1])), "--")
	plt.draw()
	plt.pause(0.01)
#plt.show()

	#plt.plot(u[:,0]/u1[:,0])
	#plt.show()
"""

plt.plot(speed[:,0])
plt.plot(speed[:,1]*4/3)
plt.show()
#np.save("results_oscwavychannel/test_kinetic_2.npy", speed)



speed = np.load("results_oscwavychannel/test_kinetic_2.npy")

Lx = 452.19 
T = 5 
epsilon = 0.017 
nu = 1 
D = 0.3 
F0 = 1 
dt = 0.02 
Sc = nu/D
pi = np.pi 
omega = 2*pi/T
kappa = 2*pi/Lx 
Nt = timesteps + 1

T = np.linspace(0, 2*np.pi/omega, Nt)
eta = np.linspace(0, 2*np.pi/kappa, 100)
xi = np.linspace(-1, 1, 500)

gamma = np.sqrt(1j*omega/Sc)
avg_kinetic = np.zeros(len(T))

for i in range(len(T)):
	avg_kinetic[i] = np.mean(np.mean(np.sqrt(np.square(u_x[i,:,:])+np.square(u_y[i,:,:]))))

import scipy.interpolate as sci 
x = np.linspace(0, len(avg_kinetic)-1, len(avg_kinetic))
interpool_anal = sci.interp1d(x, avg_kinetic, kind='cubic')
interpool_nume = sci.interp1d(x, speed, kind='cubic')
x = np.linspace(0, len(avg_kinetic)-1, int(1e5))
anal = interpool_anal(x)
nume = interpool_nume(x)

cut_an = np.argmin(abs(anal)[:int(0.3*1e5)])
cut_nu = np.argmin(abs(nume)[:int(0.3*1e5)])

plt.plot(x-x[cut_an], anal[:]/max(anal))
plt.plot(x-x[cut_nu], nume[:]/max(nume), "--")
plt.show()
print("mean analytic kinetic:", max(avg_kinetic)/max(speed))
print("mean numerical: ", max(speed))
"""