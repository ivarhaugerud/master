import numpy as np
import matplotlib.pyplot as plt

dt = 0.01
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

periods = 2000
datafiles = periods*100
skip = int(periods*timesteps/datafiles)

U_scale = 50
Pe = 10
D = U_scale/Pe

"""	
skip = 100
for i in range(int(datafiles/skip)):
	plt.clf()
	pos = np.load("data/run_08_01/RW_positions_"+str(int(skip*i))+".npy")
	plt.scatter(pos[0, :], pos[1, :])
	plt.pause(0.01)
plt.show()
"""
var = np.zeros(datafiles)

t = np.linspace(0, period*periods, datafiles)
x = np.zeros(datafiles)
y = np.zeros(datafiles)

for i in range(datafiles):
	pos = np.load("data/run_10_01/RW_positions_"+str(int(i*skip))+".npy")
	#plt.scatter(pos[0, :], pos[1, :])
	#plt.pause(0.01)
	x[i] = pos[0, 5]
	y[i] = pos[1, 5]
	var[i] = np.std(pos[0, :])
#plt.show()

plt.plot(np.trim_zeros(x), np.trim_zeros(y), "o")
plt.show()

plt.plot(t[:len(np.trim_zeros(var))]/tau, np.trim_zeros(var))
plt.xlabel("time [periods]")
plt.ylabel("variance")
plt.show()

plt.plot(t[1:len(np.trim_zeros(var))]/tau, np.trim_zeros(var[1:])/t[1:len(np.trim_zeros(var))])
plt.xlabel("time [periods]")
plt.ylabel("variance")
plt.show()

