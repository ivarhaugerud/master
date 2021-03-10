import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
from scipy.signal import savgol_filter
import os 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

numeric = np.load("data_test/vary_omega.npy")
tau     = np.logspace(-2, 2, 10)
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.2])
D0_ana  = np.zeros(len(tau))

epsilon = 0.3
nu  = 1.3
D   = 0.7
F0  = 12/nu 
Sc = nu 
Pe = 1/D


for i in range(len(tau)):
	T = tau[i]
	omega = 2*np.pi/T
	gamma = np.sqrt(1j*omega/Sc)
	gamma_c = np.conj(gamma)
	rho = np.sqrt(1j*omega/D)
	rho_c = np.conj(rho)

	D0_ana[i] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))


D_eff  = np.zeros((len(tau), len(kappa)))
datapoints = 750

for i in range(len(tau)):
	for j in range(len(kappa)):
		t = np.trim_zeros(numeric[i, j, :, 0])
		D = np.trim_zeros(numeric[i, j, :, 8])

		plt.plot(t, D)
		plt.plot(t[-datapoints:], D[-datapoints:], "--")

		D_eff[i, j]   = sci.trapz(D[-datapoints:], t[-datapoints:])/tau[i]
plt.show()


plt.figure(1)
for i in range(len(kappa)):
	plt.errorbar(2*np.pi/tau, D_eff[:, i], label=r"$\kappa=%3.2f$" % kappa[i])

plt.plot(2*np.pi/tau, D0_ana, "k", label=r"$\epsilon=0$")
plt.xscale("log")
plt.legend(loc="best", ncol=1, fontsize=8)
plt.xlabel(r"Frequency $\rho$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = "figures/D_eff_vs_rho.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()