import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scp
from numpy import *
"""
F0 = 18
r = np.logspace(-5, 5, 3000)
Ds = np.logspace(-3, 3, 7)
for i in range(len(Ds)):
	D = Ds[i]
	Pe = 1/D
	rho = np.sqrt(1j*r/D)
	gamma = np.sqrt(1j*r/1)
	gamma_c = np.conjugate(gamma)
	rho_c = np.conjugate(rho)
	D0 = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
	amp = ( (1/(rho*tanh(rho)) - 1/(sinh(rho)**2) + 1/(gamma*tanh(gamma)) - 1/(sinh(gamma)**2) - 4*(rho/tanh(rho) - gamma/tanh(gamma))/(rho*rho-gamma*gamma))/2)
	amp *= Pe*Pe*F0*F0*tanh(gamma)*tanh(gamma)/(4*gamma*gamma*(rho*rho-gamma*gamma)**2)

	plt.figure(2)
	plt.plot(r, D0, "--")
	plt.fill_between(r, D0+np.sqrt(4*np.real(amp)**2 + 4*np.imag(amp)**2), D0-np.sqrt(4*np.real(amp)**2 + 4*np.imag(amp)**2), alpha=0.1)
	plt.xscale("log")
	plt.yscale("log")

	#plt.figure(3)
	#plt.fill_between(r, np.sqrt(4*np.real(amp)**2 + 4*np.imag(amp)**2), D0-np.sqrt(4*np.real(amp)**2 + 4*np.imag(amp)**2), alpha=0.1)
	#plt.xscale("log")
	plt.show()
"""
#t = np.linspace(0, 10, int(1e3))
"""
for i in range(len(r)):
	t = np.linspace(2*np.pi/r[i], 4*np.pi/r[i], int(1e3))
	plt.plot(t, np.real(np.exp(1j*r[i]*t)*amp[i]))
plt.xscale("log")
plt.show()
"""
from sympy import * 

#define variables
F0, Pe, r, g     = symbols("F0 Pe r g", real=True)
rho, gamma = symbols("rho gamma")

#rho = sqrt(I*r)
#gamma = sqrt(I*g)
amp = (1/(rho*tanh(rho)) - 1/(sinh(rho)**2) + 1/(gamma*tanh(gamma)) - 1/(sinh(gamma)**2) - 4*(rho/tanh(rho) -gamma/tanh(gamma))/(gamma*gamma-rho*rho))/2
amp *= (Pe*F0*tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma)))**2

#print(limit(amp, rho, oo))
order = 1
approx = series(amp,    r, n=order).removeO()
approx = series(approx, g, n=order).removeO()
print(simplify(approx.subs(F0, sqrt(18))))

my_series = 2/3 - 4*rho*rho/45 + 4*rho**4/315  + 2/3 - 4*gamma*gamma/45 + 4*gamma**4/315 
my_series += (rho*rho/3 - rho**4/45 + 2*rho**6/945 - gamma*gamma/3 + gamma**4/45 -2*gamma**6/945 -rho**8/4725 + gamma**8/4725)*4/(gamma*gamma-rho*rho)
my_series = (simplify(my_series))
my_series = gamma**4 + rho**4 + gamma**4*rho**2 + gamma**2*rho**4
my_series /= (gamma*gamma - rho*rho)**2
my_series = series(my_series, gamma, n=4).removeO()
#my_series = series(my_series, rho, n=4).removeO()
print(simplify(my_series))


