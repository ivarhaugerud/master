import numpy as np 
import matplotlib.pyplot as plt 
from cmath import *
N = int(1e4)
xi = np.linspace(-1, 1, N)
H = np.zeros(N, dtype="complex")

omega = 4
Sc = 1
F_0 = 1
kappa = 1
gamma = 1j*omega/Sc 

kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
P_1 = (gamma*F_0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))

H += P_1*kappa*kappa*( (kappa**4 - gamma**4)*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma)))/(2*gamma*gamma*(gamma**2-kappa**2)**2)
H += P_1*kappa*kappa*( -np.cosh(kappa)*4*gamma*gamma*kappa*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma)))/(2*gamma*gamma*(gamma**2-kappa**2)**2)
H -= P_1*kappa_p*kappa_p*np.sinh(kappa)*(xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p)-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
H -= F_0*(xi*np.sinh(gamma*xi)*(1+kappa*kappa*xi*xi/3) - np.sinh(gamma)*(1+kappa*kappa/3)*np.cosh(gamma*xi)/(np.cosh(gamma)) + gamma*xi*xi*np.cosh(gamma*xi) - gamma*np.cosh(xi*gamma))/(4*gamma*np.cosh(gamma))

T = np.linspace(0, 2*np.pi/omega, 100)

for i in range(len(T)):
	plt.clf()
	#plt.axis([-1, 1, -0.4, 0.4])
	plt.plot(xi, np.real(H))
	plt.plot(xi, np.imag(H))
	#plt.plot(xi, np.real(H)*np.cos(omega*T[i]) - np.imag(H)*np.sin(omega*T[i]))
	plt.draw()
	plt.pause(0.01)
plt.show()