import numpy as np 
import matplotlib.pyplot as plt 
import os
import scipy.integrate as sci

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'

D_parallels = np.load("data/D_parallels_kappa.npy")
print(D_parallels)

num_kappas = np.load("../RW_simulation/data/D_eff_vs_kappa.npy")
nu = 16
omega = 5/(2*np.pi)
F0 = 10
D = 0.3
Sc = nu#/D
Pe = 1/D
kappas = np.arange(0.1, 2.5, 0.4)

print(len(kappas), len(D_parallels))
gamma   = np.sqrt(1j*omega/Sc)
gamma_r = np.real(gamma)
print("resonance at", np.sqrt(2)*gamma_r)

gamma_c = np.conj(gamma)
a = np.sqrt(1j*omega)
a_c = np.sqrt(-1j*omega)
xi = np.linspace(-1, 1, int(1e5))
factor = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
D_para0 = 1+factor*0.5*sci.trapz(np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)

filename = "figures/data_para_vs_kappa.pdf"
plt.scatter(kappas, D_parallels, s=4)
plt.plot(kappas, D_parallels, "-", linewidth=1)
#plt.yscale("log")
#plt.xscale("log")
tol = 0.2*min(kappas)
#plt.axis([min(kappas)-tol, max(kappas)*1.1, min(D_parallels)*0.8, max(D_parallels)*1.2])
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"2nd order Parallel Diffusion $D_\parallel^{(2)}$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

eps = 0.2
D_para = eps*eps*np.real(D_parallels)
D_para += np.real(D_para0)

filename = "figures/data_para_tot_vs_kappa.pdf"
plt.scatter(kappas, D_para, s=4)
plt.plot(kappas, D_para, "-", linewidth=1)
plt.plot(kappas, np.ones(len(kappas))*D_para0, "-")
#plt.yscale("log")
#plt.xscale("log")
kappas = np.array([0.1, 0.4, 0.7, 1.0, 1.3, 1.6])
plt.plot(kappas, num_kappas[:, 0])
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Total Parallel Diffusion $D_\parallel + O(\epsilon^4)$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()