import numpy as np 
import matplotlib.pyplot as plt 

xi = np.linspace(-1, 1, int(1e4))
kappa = 0.2
Pe    = 3
omega = 0.01
F0   = 6
rho = np.sqrt(1j*omega)

for i in range(8):
	omega *= 2
	rho = np.sqrt(omega) 

	beta0  = F0*F0*Pe*Pe*kappa*(xi**4 - 2*xi**2 + 7/15 + 8/(kappa*kappa))/(288*rho*rho) + kappa*(1/6 - xi*xi/2 -1/(kappa*kappa))
	beta0p = Pe*Pe*F0*F0*(1-xi*xi-2/(kappa*kappa))/(24*rho*rho) + np.cosh(kappa*xi)*(F0*F0*Pe*Pe/(12*rho*rho*kappa*kappa) - 1)/np.sinh(kappa)
	beta1 = -F0*Pe*np.ones(len(xi))/(6*rho*rho)
	beta2 = F0*F0*Pe*Pe*kappa*(  (15*xi**4 - 30*xi**2 + 7)/(4320) + 1/(72*rho*rho) )/(rho*rho)

	beta1_p = F0*Pe*(-540/(rho**4) - 270*xi*xi/(rho*rho) + 30/(rho*rho) - 7 + 30*xi*xi -15*xi**4 )/360 + 3*F0*Pe*np.cosh(rho*xi)/(2*rho*rho*rho*np.sinh(rho))
	#beta2 = F0*F0*Pe*Pe*(1/(rho*rho*kappa*360) - 1/(180*kappa*kappa*kappa) + rho*rho/(90*kappa**5))*np.ones(len(xi))
	#plt.plot(xi, beta2-beta2[0], "--")
	#beta2 = F0**2*Pe**2*(2*kappa**6*rho**2*(21*xi**6 - 105*xi**4 + 147*xi**2 - 31) + kappa**4*(kappa**2*(21*kappa**2*xi**6 - 105*kappa**2*xi**4 + 147*kappa**2*xi**2 - 31*kappa**2 + 630*xi**4 - 1260*xi**2 + 294) + 5040) - 10080*kappa**2*rho**2 + 20160*rho**4)/(181440*kappa**5*rho**2)
	#plt.plot(xi, beta0-beta0[0], label=r"$\beta=0$")
	#plt.plot(xi, beta0p-beta0p[0], "--", label=r"$\beta=0$")
	plt.plot(xi, beta1-beta1[0], label=r"$\beta=1$")
	plt.plot(xi, beta1_p-beta1_p[0], "--", label=r"$\beta=1$")
	#plt.plot(xi, beta2-beta2[0], label=r"$\beta=2$")
	plt.legend(loc="best")
	plt.show()



F0**2*Pe**2*(2*kappa**6*rho**2*(21*xi**6 - 105*xi**4 + 147*xi**2 - 31) + kappa**4*(kappa**2*(21*kappa**2*xi**6 - 105*kappa**2*xi**4 + 147*kappa**2*xi**2 - 31*kappa**2 + 630*xi**4 - 1260*xi**2 + 294) + 5040) - 10080*kappa**2*rho**2 + 20160*rho**4)/(181440*kappa**5*rho**2)
