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

D_parallels = np.load("data/D_parallels_kappa.npy")[0, :]
kappas =  np.load("data/D_parallels_kappa.npy")[1, :]

#print(D_parallels)
#Lx = np.linspace(0.7, 27, 20)
#kappas = 2*np.pi/Lx

#num_kappas = np.load("../RW_simulation/data/D_eff_vs_kappa.npy")
nu = 1.2
tau = 3.0
omega = tau/(2*np.pi)
dt = 0.006
timesteps = int(tau/dt)
F0 = 12/nu
Sc = nu
Pe = 1
Lx = np.linspace(0.7, 27, 20)
Lx =  np.array([1.05, 2.09, 6.28, 12.56, 15.71, 25.13]) #9.42,
#kappas = 2*np.pi/Lx
#kappas =  np.linspace(0.5, 3, 15) #9.42,


gamma   = np.sqrt(1j*omega/Sc)
gamma_r = np.real(gamma)
print("resonance at", np.sqrt(2)*gamma_r)
analytic = np.load("../sympy/data/analytic_Dpara.npy")

gamma_c = np.conj(gamma)
a = np.sqrt(1j*omega)
a_c = np.sqrt(-1j*omega)
xi = np.linspace(-1, 1, int(1e5))
factor = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
D_para0 = 1+factor*0.5*sci.trapz(np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)

filename = "figures/data_para_vs_kappa.pdf"
plt.scatter(kappas, D_parallels, s=4)
plt.plot(kappas, D_parallels, "-", linewidth=1)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"2nd order Parallel Diffusion $D_\parallel^{(2)}$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

eps = 0.25
D_para = eps*eps*np.real((D_parallels))
D_para += np.real(D_para0)

filename = "figures/data_para_tot_vs_kappa.pdf"
plt.scatter(kappas, D_para, s=4)
plt.plot(kappas, D_para, "-", linewidth=1)
plt.plot(kappas, np.ones(len(kappas))*D_para0, "-")

Lx =  np.array([1.05, 2.09, 6.28, 12.56, 15.71, 25.13]) #9.42,
kappas = 2*np.pi/Lx
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

Gaute_deff = np.load("data/tdatas.npy")
kappa  = np.array([0.5, 0.9, 1.1, 1.3, 1.7, 2.1])

gate_D = np.zeros(len(kappa))

for i in range(len(Gaute_deff[:,0,0])):
	plt.plot(Gaute_deff[i, :, 0], Gaute_deff[i, :, 8])
	plt.plot(Gaute_deff[i, -timesteps:, 0], Gaute_deff[i, -timesteps:, 8], "--")
	gate_D[i] = integrate.trapz(Gaute_deff[i, -timesteps:, 8], Gaute_deff[i, -timesteps:, 0])/tau
plt.show()

plt.plot(kappa, gate_D)
plt.show()