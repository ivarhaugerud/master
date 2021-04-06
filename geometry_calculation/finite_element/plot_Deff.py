import numpy as np 
import matplotlib.pyplot as plt 
import os
import scipy.integrate as sci
from scipy import integrate

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'

D_parallels = np.load("data/D_parallels_kappa_D1.npy")[:]

dt  = 0.04
nu  = 1.2
D   = 1
F0  = 12/nu 
tau = 3.0
omega = 2*np.pi/tau
Sc = nu 
Pe = 1/D
gamma = np.sqrt(1j*omega/Sc)
gamma_c = np.conj(gamma)
rho = np.sqrt(1j*omega/D)
rho_c = np.conj(rho)
kappas   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])

D_para0 = np.real(1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c))) )

filename = "figures/data_para_vs_kappa.pdf"
plt.scatter(kappas, D_parallels, s=4)
plt.plot(kappas, D_parallels, "-", linewidth=1)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"2nd order Parallel Diffusion $D_\parallel^{(2)}$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
#plt.show()

"""
numeric = np.load("../data_test/tdata_03_03_D01_.npy")
epsilon = np.array([0.1, 0.2, 0.3, 0.5]) #np.arange(0.05, 0.51, 0.05)
kappa_num   = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1]) #0.2
D_num = np.zeros((len(epsilon), len(kappa_num)))
T = int(3/0.004)

for i in range(len(epsilon)):
	for j in range(len(kappa_num)):
		#plt.plot(numeric[i, j, :, 0]/tau, numeric[i, j, :, 8], label=r"$\kappa=$"+str(kappa_num[j]))
		#plt.plot(numeric[i, j, -T:, 0]/tau, numeric[i, j, -T:, 8])
		D_num[i, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)

	#plt.show()

for i in range(len(epsilon)-1):
	eps = epsilon[i]
	plt.plot(kappas, D_para0+epsilon[i]*epsilon[i]*D_parallels,  color="C"+str(i))
	plt.plot(kappa_num, D_num[i,:], "o", markersize=2.5, color="C"+str(i), label=r"$\epsilon=$"+str(epsilon[i]))

plt.plot(kappa_num, D_para0*np.ones(len(kappa_num)), "k")
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.savefig("figures/approximate_result.pdf")
plt.show()

for i in range(len(kappas)):
	plt.plot(epsilon, D_para0+epsilon*epsilon*D_parallels[i],  color="C"+str(i))
	plt.plot(epsilon, D_num[:,i], "o", color="C"+str(i), label=r"$\epsilon=$"+str(kappas[i]))

plt.plot(epsilon, D_para0*np.ones(len(epsilon)), "k")
plt.show()


for i in range(len(kappas)):
	plt.plot(epsilon, (D_num[:,i]-D_para0*np.ones(len(epsilon))), color="C"+str(i), label=r"$\epsilon=$"+str(kappas[i]))
plt.show()
"""


numeric = np.load("../data_test/tdata_04_03_D1_.npy")
epsilon = np.array([0.0, 0.1, 0.2, 0.3])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
T = 3
dt = 0.004
T = int(tau/dt)
D_num[0, :] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
plt.show()

for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		#plt.plot(numeric[i, j,   :, 0], numeric[i, j,   :, 8])
		#plt.plot(numeric[i, j, -T:, 0], numeric[i, j, -T:, 8])
		if i == 2:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -2*T:-T, 8], numeric[i, j, -2*T:-T, 0])/(tau)
		else:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)
plt.figure(3)
for i in range(len(epsilon)):
	plt.plot(kappa, D_num[i,:], color="C"+str(i), label=r"$\kappa=$"+str(kappa[i]))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.show()


for i in range(len(epsilon)):
	plt.plot(kappas, D_para0+epsilon[i]*epsilon[i]*D_parallels,  color="C"+str(i))
	plt.plot(kappas, D_num[i,:], "o", markersize=2.5, color="C"+str(i), label=r"$\epsilon=$"+str(epsilon[i]))

plt.plot(kappas, D_para0*np.ones(len(kappas)), "k")
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.savefig("figures/approximate_result.pdf")
plt.show()

for i in range(len(kappas)):
	plt.plot(epsilon, D_para0+epsilon*epsilon*D_parallels[i],  color="C"+str(i))
	plt.plot(epsilon, D_num[:,i], "o", color="C"+str(i), label=r"$\epsilon=$"+str(kappas[i]))

plt.plot(epsilon, D_para0*np.ones(len(epsilon)), "k")
plt.show()


for i in range(len(kappas)):
	plt.plot(epsilon, (D_num[:,i]-D_para0*np.ones(len(epsilon))), color="C"+str(i), label=r"$\epsilon=$"+str(kappas[i]))
plt.show()
