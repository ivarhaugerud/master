import numpy as np 
import matplotlib.pyplot as plt 

Sc = 1
F0 = 1
pi = np.pi 
omega = 1
kappa = 1
Nt = 5
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(0, 6*np.pi/kappa, 300)
xi = np.linspace(-1, 1, 100)

gamma = np.sqrt(1j*omega/Sc)

u_x = np.zeros((len(T), len(xi), len(eta)), dtype="complex")
#psi_1 = F0*gamma*np.tanh(gamma)*np.cosh(kappa)/(kappa*np.cosh(2*kappa))
#phi_1 = - psi_1*np.sinh(kappa)/np.cosh(kappa)

phi_1 = -F0*gamma*np.tanh(gamma)/(np.sinh(kappa)+np.cosh(kappa)*np.cosh(kappa)/np.sinh(kappa))
psi_1 = F0*gamma*np.tanh(gamma)*np.cosh(kappa)/(np.sinh(kappa)*np.sinh(kappa) + np.cosh(kappa)*np.cosh(kappa))

for t in range(len(T)):
	for x in range(len(xi)):
		u_x[t, x, :] -= F0*xi[x]*np.sinh(gamma*xi[x])*np.sin(kappa*eta)/(gamma*np.cosh(gamma))
		u_x[t, x, :] -= Sc*kappa*np.cosh(kappa*xi[x])*(-psi_1*np.sin(kappa*eta) + phi_1*np.cos(kappa*eta))/(1j*omega)
		u_x[t, x, :] *= np.exp(1j*omega*T[t])

for i in range(len(T)):
	plt.imshow(np.real(u_x[i, :,:]))
	plt.show()

