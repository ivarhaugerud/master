import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata
import matplotlib.ticker as tick 

name1 = "results_oscwavychannel/Lx0.1_tau5.0_eps0.0_nu1.0_D0.3_fzero0.0_fone3.0_res100_dt0.02/u.h5" #results_oscwavychannel/Lx452.19_tau5.0_eps0.032_nu1.0_D0.3_fzero0.0_fone3.0_res106_dt0.02/u.h5"

f = h5py.File(name1, 'r')

dt1 = 0.02
tau1 = 5
timesteps = int((tau1/0.02))
speed = np.zeros((timesteps+1, 2))


max_index = len(list(f["VisualisationVector"]))
geometry = np.array(list(f["Mesh"]["0"]["mesh"]["geometry"]))
eta = geometry[::11, 0]
xi  = geometry[::11, 1]

F0 = 3
omega = 2*np.pi/tau1
nu = 1
D = 0.3
Sc = nu/D
gamma = np.sqrt(1j*omega/Sc)
T = np.linspace(3*2*np.pi/omega, 0, max_index+1)
index_1_p = np.argmin(abs(T-2*2*np.pi/omega))
spatial_avg = np.zeros((index_1_p, 2))


tdat = np.loadtxt("results_oscwavychannel/Lx0.1_tau5.0_eps0.0_nu1.0_D0.3_fzero0.0_fone3.0_res100_dt0.02/tdata.dat")
stop_index = np.argmin(abs(tdat[:,0]-(3*tau1-2.5*tau1)))
end_index  = np.argmin(abs(tdat[:,0]-(3*tau1-0.5*tau1)))
# x = [time, <ux>, <ux/U>, <uy>, <u*u>, ...]

from scipy import integrate
for i in range(max_index-1):
	if i < index_1_p:
		#plt.clf()
		u = np.array(list(f["VisualisationVector"][str(int((max_index-i-1)))]))

		x = np.linspace(-1, 1, int(1e3))#len(u[::11,0]))
		spatial_avg[i, 0] = integrate.trapz(u[::11,0]*u[::11,0], xi)/2
		u_ana = np.real(F0*((1-np.cosh(gamma*x)/np.cosh(gamma))/(gamma*gamma))*np.exp(1j*omega*T[max_index-i-1]))

		spatial_avg[i, 1] = integrate.trapz(u_ana*u_ana, x)/2
		#plt.plot(xi, u[::11,0])
		#plt.plot(x, np.real(F0*((1-np.cosh(gamma*x)/np.cosh(gamma))/(gamma*gamma))*np.exp(1j*omega*T[max_index-i-1])), "--")
		#plt.draw()
		#plt.pause(0.01)
#plt.show()

plt.plot(tdat[stop_index:end_index, 4])
plt.plot(spatial_avg[:,0])
plt.plot(spatial_avg[:,1], "--")
plt.show()	

plt.plot(spatial_avg[:,0]-spatial_avg[:,1])
plt.show()	

plt.plot(abs(spatial_avg[:,0]-spatial_avg[:,1]/spatial_avg[:,1]))
plt.show()	

time_avg_num2 = integrate.trapz(tdat[stop_index:end_index, 4], tdat[stop_index:end_index, 0])/(4*np.pi/omega)
time_avg_num = integrate.trapz(spatial_avg[:,0], np.flip(T[:index_1_p]))/(2*np.pi/omega)
time_avg_ana = integrate.trapz(spatial_avg[:,1], np.flip(T[:index_1_p]))/(2*np.pi/omega)
print(time_avg_ana, time_avg_num, time_avg_num2)
	#plt.plot(u[:,0]/u1[:,0])
	#plt.show()
"""
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