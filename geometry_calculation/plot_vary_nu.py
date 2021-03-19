import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

Lx = np.array([4.48])
kappa = 2*np.pi/Lx
nu = np.logspace(-2, 2, 10)[:9]

base = "data_test/vary_nu/"
D = np.zeros((len(nu)))
difference = np.zeros(np.shape(D))
U  = np.zeros(np.shape(D))

T = 750
tau = 3.0

print(nu)

for i in range(len(nu)):
	data = np.loadtxt(base+"Lx4.48_tau3.0_eps0.5_nu" + str(nu[i])[:4]+"_D1.0_fzero0.0_fone" + str(12*nu[i])[:5] + "_res150_dt0.004/tdata.dat")
	print(np.shape(data), np.shape(D))
	D[i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
	U[i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
	#plt.scatter(nu[i], data[-1, 0]/tau)
	difference[i] = abs(D[i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[i]
	#plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
	#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
	#plt.show()
#plt.xscale("log")
#plt.show()
plt.figure(1)
plt.plot(nu, difference, "o")
plt.ylabel(r"check for convergence")
plt.xlabel(r"viscosity $\nu$")
plt.xscale("log")
plt.legend(loc="best")

plt.figure(2)
plt.plot(nu, D, "o")
plt.ylabel(r"$D_{eff}$")
plt.xlabel(r"viscosity $\nu$")

plt.xscale("log")
plt.legend(loc="best")


plt.figure(3)
plt.plot(nu, U, "o")
plt.xscale("log")
plt.ylabel(r"$\langle u^2 \rangle $")
plt.xlabel(r"viscosity $\nu$")
plt.legend(loc="best")


plt.figure(4)
plt.plot(nu, abs((D)/(U*U)), "o")
plt.ylabel(r"$D/U^2$")
plt.xlabel(r"viscosity $\nu$")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="best")

plt.show()
