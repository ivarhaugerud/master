import numpy as np
import matplotlib.pyplot as plt

dt = 0.002
tau = 3.0 

epsilon = 0.0
U_scale = 1
Pe = 3.0
D = U_scale/Pe

kappas = np.array([0.4])
Lx = 2*np.pi/kappas

omega = 2*np.pi/tau
nu = 4.0
F0 = 3/nu
Sc = nu
gamma = np.sqrt(1j*omega/Sc)
timesteps = int(tau/dt)
period    = tau

periods = 5000
datafiles = periods*25
skip = int(periods*timesteps/datafiles)

import scipy.integrate as sci
gamma_c = np.conj(gamma)
a       = np.sqrt(1j*omega)
a_c     = np.conj(a)
xi      = np.linspace(-1, 1, int(1e5))

factor  = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
D_para0 = 1 + factor*0.5 * sci.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)
ux = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)

print("D_parallel: ", D_para0)
print("U_avg: ", 0.5*sci.trapz(ux*np.conj(ux), xi))

var  = np.zeros(datafiles)
mean = np.zeros(datafiles)

t = np.linspace(0, period*periods, datafiles)
x = np.zeros(datafiles)
y = np.zeros(datafiles)

for i in range(datafiles):
	#plt.clf()
	pos = np.load("flow_fields/zero_eps/mu_4.0/pos/RW_positions_"+str(int(i*skip))+".npy")
	#plt.scatter(pos[0, :], pos[1, :])
	#plt.pause(0.5)
	#plt.axis([-3.5, 3.5, -1, 1])
	x[i] = pos[0, 5]
	y[i] = pos[1, 5]
	var[i]  = np.square(np.std(pos[0, :]))/D
	mean[i] = np.mean(pos[0, :])

np.save("flow_fields/zero_eps/nu_0.5/pos/var", var)
#plt.show()

plt.plot(t, mean)
plt.show()

plt.plot(x, y)
plt.show()

plt.plot(t/tau, var)
plt.xlabel("time [periods]")
plt.ylabel("variance")
plt.show()

plt.plot(t[1:len(var)]/tau, var[1:]/t[1:len(var)])
plt.plot(t/tau, np.ones(len(t))*D_para0)
plt.xlabel("time [periods]")
plt.ylabel("variance")
plt.show()