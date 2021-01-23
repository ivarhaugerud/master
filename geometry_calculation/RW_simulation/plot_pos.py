import numpy as np
import matplotlib.pyplot as plt

dirr = "flow_fields/non_zero_eps/Lx15.71_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/"
dt = 0.006
tau = 3.0

epsilon = 0.25
U  = 0.6204847210770329
U = 1.1772393257809046#0.03586352057013804
U = 0.4743218639154592
Pe = 6
D  = U/Pe 

#Pe      = #U_scale/D

kappas = np.array([0.4])
Lx = 2*np.pi/kappas

omega = 2*np.pi/tau
nu = 3.6	
F0 = 12/nu
Sc = nu
gamma = np.sqrt(1j*omega/Sc)
timesteps = int(tau/dt)
period    = tau
periods = 4000
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
	pos = np.load(dirr+"pos/RW_positions_"+str(int(i*skip))+".npy")
	#plt.title(str(t[int(i)]))
	#plt.scatter(pos[0, :], pos[1, :])
	#plt.axis([-4, 4, -1, 1])
	#plt.pause(0.5)
	x[i] = pos[0, 5]
	y[i] = pos[1, 5]
	var[i]  = np.square(np.std(pos[0, :]))
	mean[i] = np.mean(pos[0, :])

print(np.mean(var[int(datafiles/2):]/t[int(datafiles/2):]))
print(np.mean(np.std((var[int(datafiles/2):]/t[int(datafiles/2):]))))

np.save(dirr+"pos/var", var)

plt.plot(t, mean)
plt.show()

plt.plot(x, y)
plt.show()

plt.plot(t/tau, var)
plt.xlabel("time [periods]")
plt.ylabel("variance")
plt.show()

plt.plot(t[1:len(var)]/tau, var[1:]/t[1:len(var)]/D)
plt.plot(t/tau, np.ones(len(t))*D_para0)
plt.xlabel("time [periods]")
plt.ylabel("variance")
plt.show()