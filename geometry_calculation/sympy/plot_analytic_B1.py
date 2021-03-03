import numpy as np 
import matplotlib.pyplot as plt 

xi = np.linspace(-1, 1, int(1e4))
kappa = 0.2
Pe    = 3
omega = 0.2
F0   = 6
rho = np.sqrt(1j*omega)
kappa = 0

for i in range(8):
	omega *= 2
	rho = np.sqrt(omega) 
	rho_p = np.sqrt(kappa*kappa+rho*rho)
	gamma = np.sqrt(omega/100)
	#beta1 = F0*Pe*(-7 + 30*xi*xi - 15*xi**4 + (30-270*xi*xi)/(rho*rho) -540 /(rho**4))/360
	#beta1 += 3*F0*Pe*np.cosh(rho*xi)/(2*rho*rho*rho*np.sinh(rho))
	#plt.plot(xi, beta1, label=r"$\beta=1$")
	#beta1 = F0*Pe*(rho**2*(rho**2*(21*xi**6 - 105*xi**4 + 147*xi**2 - 31) + 210*xi**4 - 420*xi**2 + 98) - 1680)/(10080*rho**2)
	#beta1_p = F0*Pe*( - 1680)/(10080*rho**2)*np.ones(len(xi))	
	#beta0  = F0*F0*Pe*Pe*kappa*(xi**4 - 2*xi**2 + 7/15 + 8/(kappa*kappa))/(288*rho*rho) + kappa*(1/6 - xi*xi/2 -1/(kappa*kappa))
	#beta0p = Pe*Pe*F0*F0*(1-xi*xi-2/(kappa*kappa))/(24*rho*rho) + np.cosh(kappa*xi)*(F0*F0*Pe*Pe/(12*rho*rho*kappa*kappa) - 1)/np.sinh(kappa)
	#beta1 = -F0*Pe*np.ones(len(xi))/(6*rho*rho)
	#beta2 = F0*F0*Pe*Pe*kappa*(  (15*xi**4 - 30*xi**2 + 7)/(4320) + 1/(72*rho*rho) )/(rho*rho)
	#B0 = F0*Pe*(gamma**2*rho**2*( 630*xi**4 - 1260*xi**2 + 294) + 5040*gamma**2 - 15120)/(30240*gamma**2*rho**2)
	#B0_p = 0.5*Pe*F0*np.tanh(gamma)/(gamma*(rho*rho-gamma*gamma))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma)))
	beta1_p = F0*Pe*(-540/(rho**4) - 270*xi*xi/(rho*rho) + 30/(rho*rho) - 7 + 30*xi*xi -15*xi**4 )/360 + 3*F0*Pe*np.cosh(rho*xi)/(2*rho*rho*rho*np.sinh(rho))
	B1 = -rho*rho*xi**4/(24*rho_p*rho_p) + xi*xi*(rho*rho - 3 - 6*rho*rho/(rho_p*rho_p))/(12*rho_p*rho_p)  - (12*rho*rho/(rho_p**4) - 2*rho*rho/(rho_p*rho_p) + 6/(rho_p*rho_p) + 1 + 7*rho*rho/30)/(12*rho_p*rho_p)
	B1 *= Pe*F0
	#beta2 = F0*F0*Pe*Pe*(1/(rho*rho*kappa*360) - 1/(180*kappa*kappa*kappa) + rho*rho/(90*kappa**5))*np.ones(len(xi))
	#plt.plot(xi, beta2-beta2[0], "--")
	#beta2 = F0**2*Pe**2*(2*kappa**6*rho**2*(21*xi**6 - 105*xi**4 + 147*xi**2 - 31) + kappa**4*(kappa**2*(21*kappa**2*xi**6 - 105*kappa**2*xi**4 + 147*kappa**2*xi**2 - 31*kappa**2 + 630*xi**4 - 1260*xi**2 + 294) + 5040) - 10080*kappa**2*rho**2 + 20160*rho**4)/(181440*kappa**5*rho**2)
	plt.plot(xi, B1-B1[0], label=r"$\beta=0$")
	plt.plot(xi, beta1_p-beta1_p[0], "--", label=r"$\beta=0$")
	#plt.plot(xi, beta0p-beta0p[0], "--", label=r"$\beta=0$")
	#plt.plot(xi, beta1, label=r"$\beta=1$")
	#plt.plot(xi, beta1_p, "--", label=r"$\beta=1$")
	print(omega)
	#plt.plot(xi, beta2-beta2[0], label=r"$\beta=2$")
	plt.legend(loc="best")
	plt.show()



F0**2*Pe**2*(2*kappa**6*rho**2*(21*xi**6 - 105*xi**4 + 147*xi**2 - 31) + kappa**4*(kappa**2*(21*kappa**2*xi**6 - 105*kappa**2*xi**4 + 147*kappa**2*xi**2 - 31*kappa**2 + 630*xi**4 - 1260*xi**2 + 294) + 5040) - 10080*kappa**2*rho**2 + 20160*rho**4)/(181440*kappa**5*rho**2)
