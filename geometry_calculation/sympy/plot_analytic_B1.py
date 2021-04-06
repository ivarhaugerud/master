import numpy as np 
import matplotlib.pyplot as plt 

xi = np.linspace(-1, 1, int(1e4))
kappa = 0.2
Pe    = 3
omega = 0.2
F0   = 6
nu = 100
D = 20

for i in range(8):
	omega *= 2
	rho = np.sqrt(1j*omega/D) 
	gamma = np.sqrt(1j*omega/nu)

	B0_grad = F0*Pe*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))
	B0_grad_approx = F0*Pe*(np.sinh(rho*xi)/np.sinh(rho)-1)/(2*rho*rho)
	B0_grad_approx = F0*Pe*(np.cosh(rho*xi)/np.sinh(rho))/(2*rho)
	ux0 = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
	RHS_full = Pe*ux0 - 2*np.gradient(B0_grad, xi)

	plt.plot(xi, RHS_full)
	plt.plot(xi, Pe*F0*( (1-xi*xi)/4 - np.cosh(rho*xi)/(rho*np.sinh(rho)) ), "--")
	plt.legend(loc="best")
	plt.show()


