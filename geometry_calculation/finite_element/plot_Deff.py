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

D_parallels = np.load("data/D_parallels_kappa.npy")[:]
kappas =  np.load("data/D_parallels_kappa.npy")[:]

#print(D_parallels)
#Lx = np.linspace(0.7, 27, 20)
#kappas = 2*np.pi/Lx

#num_kappas = np.load("../RW_simulation/data/D_eff_vs_kappa.npy")
nu = 1.2
tau = 3.0
omega = 2*np.pi/tau
dt = 0.006
timesteps = int(tau/dt)
F0 = 12/nu
Sc = nu
D = 0.1
Pe = 1/D
kappas = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1]) #0.2

#kappas = 2*np.pi/Lx
#kappas =  np.linspace(0.5, 3, 15) #9.42,


gamma   = np.sqrt(1j*omega/Sc)
rho     = np.sqrt(1j*omega/D)
rho_c   = np.conj(rho)
gamma_c = np.conj(gamma)

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

eps = 0.25
D_para = eps*eps*np.real((D_parallels))
D_para += np.real(D_para0)

filename = "figures/data_para_tot_vs_kappa.pdf"
plt.scatter(kappas, D_para, s=4)
plt.plot(kappas, D_para, "-", linewidth=1)
plt.plot(kappas, np.ones(len(kappas))*D_para0, "-")

#Lx =  np.array([1.05, 2.09, 6.28, 12.56, 15.71, 25.13]) #9.42,
#kappas = 2*np.pi/Lx
#D_RW = np.load("../RW_simulation/data/D_eff_vs_kappa.npy")
#plt.errorbar(kappas, D_RW[:,0], yerr=D_RW[:,1], fmt="o")
#plt.fill_between(kappas, D_RW[:,0]-eps*eps*eps*eps, D_RW[:,0]+eps*eps*eps*eps, alpha=0.7, color="green")

plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Total Parallel Diffusion $D_\parallel + O(\epsilon^4)$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()



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
	plt.plot(kappa_num, D_num[i,:], "o", color="C"+str(i), label=r"$\epsilon=$"+str(epsilon[i]))

plt.plot(kappa_num, D_para0*np.ones(len(kappa_num)), "k")
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.xlabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.show()

for i in range(len(kappas)):
	plt.plot(epsilon, D_para0+epsilon*epsilon*D_parallels[i],  color="C"+str(i))
	plt.plot(epsilon, D_num[:,i], "o", color="C"+str(i), label=r"$\epsilon=$"+str(kappas[i]))

plt.plot(epsilon, D_para0*np.ones(len(epsilon)), "k")
plt.show()


for i in range(len(kappas)):
	plt.plot(epsilon, (D_num[:,i]-D_para0*np.ones(len(epsilon)))/(epsilon*epsilon), color="C"+str(i), label=r"$\epsilon=$"+str(kappas[i]))
plt.show()