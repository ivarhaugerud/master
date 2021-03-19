import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci


base = "data_test/taylor_approximation/"
tau  = np.logspace(-2, 2, 10)
T    = 500
dt   = tau/T
D  = np.zeros((len(tau)))
U  = np.zeros(np.shape(D))
difference = np.zeros(np.shape(D))
D_ana = np.zeros(len(tau))

nu = 8000
epsilon = 0.4
Dm = 0.1
F0 = 96000/nu 
omega = 2*np.pi/tau 
RHO = np.sqrt(1j*omega/Dm)
Pe = 1/Dm
from numpy import *
D0 = np.zeros(len(tau))

for i in range(len(tau)):
	rho = RHO[i]
	gamma = np.sqrt(1j*omega[i]/nu)
	D_ana[i] = 0.5*F0**2*Pe**2*(0.583333333333333*rho**8*conjugate(rho)**2 + 0.5*rho**8 + 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho)) + 1.75*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 1.75*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 0.5*rho**7*conjugate(rho)**4/tanh(rho) + 0.5*rho**7*conjugate(rho)**4/tanh(rho)**3 - 2.75*rho**6*conjugate(rho)**4 - 1.0*rho**6*conjugate(rho)**2 - 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4.0*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4.0*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 2.25*rho**6*conjugate(rho)**4/tanh(rho)**2 + 1.0*rho**5*conjugate(rho)**6/tanh(rho) - 0.25*rho**5*conjugate(rho)**4/tanh(rho) - 1.0*rho**5*conjugate(rho)**6/tanh(rho)**3 + 2.75*rho**4*conjugate(rho)**6 + 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 0.25*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 2.25*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4.0*rho**4*conjugate(rho)**6/tanh(rho)**2 - 0.5*rho**3*conjugate(rho)**8/tanh(rho) + 4.0*rho**3*conjugate(rho)**6/tanh(rho) + 0.5*rho**3*conjugate(rho)**8/tanh(rho)**3 - 0.583333333333333*rho**2*conjugate(rho)**8 + 1.0*rho**2*conjugate(rho)**6 + 1.75*rho**2*conjugate(rho)**8/tanh(rho)**2 - 1.75*rho*conjugate(rho)**8/tanh(rho) - 0.5*conjugate(rho)**8)/(rho**4*(1.0*rho**6 - 3.0*rho**4*conjugate(rho)**2 + 3.0*rho**2*conjugate(rho)**4 - 1.0*conjugate(rho)**6)*conjugate(rho)**4)
	D0[i]    = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))
	try:
		data = np.loadtxt(base+"Lx31.41_tau" + str(tau[i])[:4] + "_eps0.5_nu8000.0_D0.1_fzero0.0_fone96000.0_res150_dt"+str(dt[i])[:6] +"/tdata.dat")
		print(np.shape(data), np.shape(D))
		D[i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau[i]
		U[i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau[i]
	#plt.scatter(nu[i], data[-1, 0]/tau)
		difference[i] = abs(D[i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau[i])/D[i]
		plt.plot(np.trim_zeros(data[:, 0])/tau[i], np.trim_zeros(data[:, 8]))
		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[i], np.trim_zeros(data[:, 8])[-T:])
		plt.show()
	except:
		print(tau[i])
	



plt.plot(omega, D_ana)
plt.xscale("log")
plt.show()
#plt.xscale("log")
#plt.show()
plt.figure(1)
plt.plot(tau, difference, "o")
plt.ylabel(r"check for convergence")
plt.xlabel(r"viscosity $\nu$")
plt.xscale("log")
plt.legend(loc="best")

plt.figure(2)
plt.plot(tau, D, "o")
#plt.plot(tau, D0)
#plt.plot(tau, epsilon*epsilon*D_ana)
plt.plot(tau, D0 + epsilon*epsilon*D_ana)
plt.ylabel(r"$D_{eff}$")
plt.xlabel(r"viscosity $\nu$")

plt.xscale("log")
plt.legend(loc="best")
plt.show()
"""
plt.figure(3)
plt.plot(tau, U, "o")
plt.xscale("log")
plt.ylabel(r"$\langle u^2 \rangle $")
plt.xlabel(r"viscosity $\nu$")
plt.legend(loc="best")


plt.figure(4)
plt.plot(tau, abs((D)/(U*U)), "o")
plt.ylabel(r"$D/U^2$")
plt.xlabel(r"viscosity $\nu$")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="best")

plt.show()
"""
