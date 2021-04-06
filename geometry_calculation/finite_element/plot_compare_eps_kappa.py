import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
from scipy.signal import savgol_filter
import os 

root = "../../../master_latex/results/"
plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


numeric = np.load("../data_test/tdata_04_03_D01.npy")

eps_num = np.array([0.1, 0.2, 0.3])#np.arange(0.05, 0.51, 0.05)
epsilon = eps_num

kappa       = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1]) #0.2
kappa_num   = kappa#np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.5]) #0.2
print(np.shape(numeric))
tau = 3
dt  = 0.004
nu  = 1.2
D   = 0.1
F0  = 12/nu 
omega = 2*np.pi/tau
Sc = nu 
Pe = 1/D
gamma = np.sqrt(1j*omega/Sc)
gamma_c = np.conj(gamma)
rho = np.sqrt(1j*omega/D)
rho_c = np.conj(rho)
T = int(tau/dt)

D0_ana = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
D_num = np.zeros((len(epsilon), len(kappa)))
differene = np.zeros(np.shape(D_num))

print(np.shape(numeric))
for i in range(len(epsilon)):
	for j in range(len(kappa)):
		D_num[i, j]        = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)
		differene[i, j]   = abs(D_num[i,j]-sci.trapz(numeric[i, j, -2*T:-T, 8], numeric[i, j, -2*T:-T, 0])/(tau))/D_num[i,j]
		#plt.plot(numeric[i, j, :, 0], numeric[i, j, :, 8])
		#plt.plot(numeric[i, j, -T:, 0], numeric[i, j, -T:, 8])
	#plt.show()
"""
for i in range(len(kappa)):
	plt.plot(epsilon, differene[:, i])
plt.yscale("log")
plt.show()
"""
numeric = np.load("../data_test/tdata_03_03_D01_.npy")
print(np.shape(numeric))
epsilon = np.array([0.0, 0.1, 0.2, 0.3, 0.5])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
D_num[0, :] = D0_ana

for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		D_num[i+1, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)

plt.figure(2)
ana_kappa     = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
ana_deff = np.load("../finite_element/data/vary_kappa_single.npy")[:-1]
for i in range(len(kappa)):
	plt.plot(epsilon, D_num[:,i], color="C"+str(i), label=r"$\kappa=$"+str(kappa_num[i]))
	plt.plot(epsilon, epsilon*epsilon*ana_deff[i,0]+D0_ana, "o", markersize=3, color="C"+str(i))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)


plt.figure(5)
ana_kappa     = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
ana_deff = np.load("../finite_element/data/vary_kappa_single.npy")[:-1]
for i in range(len(kappa)):
	plt.plot(epsilon, 1e-5+abs(D_num[:,i]-(D0_ana+epsilon*epsilon*ana_deff[i,0]))/(D_num[:,i]), color="C"+str(i), label=r"$\kappa=$"+str(kappa_num[i]))
	#plt.plot(epsilon, epsilon*epsilon*ana_deff[i,0], "--", color="C"+str(i))
E = np.linspace(0.05, max(epsilon), int(1e4))
plt.plot(E, E**4, "ok", label=r"$\epsilon^4$", markersize=1.5)
plt.xscale("log")
plt.legend(loc="best", fontsize=8)
plt.yscale("log")
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"$(D_\parallel-D_\parallel^{(2)})$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#plt.show()
#filename = root+"figures/D_eff_vs_eps_D01.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))

numeric = np.load("../data_test/tdata_04_03_D1_.npy")
epsilon = np.array([0.0, 0.1, 0.2, 0.3])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
D = 1.0
Pe = 1/D
rho = np.sqrt(1j*omega/D)
rho_c = np.conj(rho)
D_num[0, :] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))

for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		if i == 2:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -2*T:-T, 8], numeric[i, j, -2*T:-T, 0])/(tau)
		else:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)

ana_kappa     = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
ana_deff = np.load("../finite_element/data/vary_kappa_single_D1.npy")[:-1]

plt.figure(3)
for i in range(len(epsilon)):
	plt.plot(kappa, D_num[i,:], color="C"+str(i), label=r"$\epsilon=$"+str(epsilon[i]))
	plt.plot(kappa, D_num[i,:], "o", markersize=3, color="C"+str(i))
	plt.plot(ana_kappa, epsilon[i]*epsilon[i]*ana_deff+D_num[0,0], "--", color="C"+str(i))


plt.legend(loc="best", fontsize=8, ncol=2)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)

plt.figure(4)
print(np.shape(D_num))
print(np.shape(epsilon))
print(np.shape(kappa))

for i in range(len(kappa)):
	plt.plot(epsilon, abs((1e-5+D_num[:,i]-(epsilon*epsilon*ana_deff[i]+D_num[0,0]))), color="C"+str(i), label=r"$\kappa=$"+str(kappa[i]))
	plt.plot(epsilon, abs((1e-5+D_num[:,i]-(epsilon*epsilon*ana_deff[i]+D_num[0,0]))), "o", color="C"+str(i), markersize=3)
	#plt.plot(ana_kappa, epsilon[i]*epsilon[i]*ana_deff+D_num[0,0], "--", color="C"+str(i))

E = np.linspace(0.05, 0.3, int(1e4))
plt.plot(E, E**4, "ok", label=r"$\epsilon^4$", markersize=1.5)
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="best", fontsize=8, ncol=2)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"$(D_\parallel-D_\parallel^{(2)})$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)


plt.show()
"""
numeric = np.load("../data_test/tdata_04_03_D10_.npy")
epsilon = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
tau = 30
dt  = 0.04
nu  = 12
D   = 10
F0  = 12/nu 
omega = 2*np.pi/tau
Sc = nu 
Pe = 1/D
gamma = np.sqrt(1j*omega/Sc)
gamma_c = np.conj(gamma)
rho = np.sqrt(1j*omega/D)
rho_c = np.conj(rho)
T = int(tau/dt)
D_num[0, :] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))


for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		#plt.plot(np.trim_zeros(numeric[i, j,   :, 0]), np.trim_zeros(numeric[i, j,   :, 8]))
		#plt.plot(numeric[i, j, -T:, 0], numeric[i, j, -T:, 8])
		D_num[i+1, j]   = sci.trapz(numeric[i, j, -2*T:-T, 8], numeric[i, j, -2*T:-T, 0])/(tau)
	#plt.show()

plt.figure(4)
for i in range(len(kappa)):
	plt.plot(epsilon, D_num[:,i], color="C"+str(i), label=r"$\kappa=$"+str(kappa[i]))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root+"figures/D_eff_vs_eps_D10.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()




"""