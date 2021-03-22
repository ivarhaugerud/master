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

RW_sim_old  = np.load("../data_test/final_run_RW_D0.1.npy")
RW_sim  = np.load("../data_test/RW_pos_03_03__D01.npy")

numeric = np.load("../data_test/tdata_04_03_D01.npy")
t = np.linspace(0, 300*3, len(RW_sim[0,0,:]))

EPS = np.array([0.1, 0.2, 0.3])
epsilon = np.array([0.1, 0.2, 0.3])
eps_num = epsilon #np.arange(0.05, 0.51, 0.05)

kappa       = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1]) #0.2
kappa_num   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.5]) #0.2

print(len(kappa_num), np.shape(numeric))
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

cutoff = int(len(t)*0.5)
D_RW  = np.zeros((len(epsilon), len(kappa), 2))
RW_sim2  = np.zeros((len(EPS), len(kappa), 2))

D_num = np.zeros((len(eps_num), len(kappa_num)))

for i in range(len(eps_num)):
	for j in range(len(kappa_num)):
		D_num[i, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		D_RW[i, j, 0] = np.mean(RW_sim[i, j, cutoff:]/t[cutoff:])
		D_RW[i, j, 1] = np.std( RW_sim[i, j, cutoff:]/t[cutoff:])

for i in range(len(EPS)):
	for j in range(len(kappa)):
		RW_sim2[i, j, 0] = (np.mean(RW_sim_old[i, j, cutoff:]) + D_RW[i, j, 0])/2
		RW_sim2[i, j, 1] = np.sqrt( np.std( RW_sim_old[i, j, cutoff:])**2 + D_RW[i, j, 1]**2 )


plt.figure(1)
for i in range(len(epsilon)):
	plt.errorbar(kappa, RW_sim2[i, :, 0], yerr=RW_sim2[i, :, 1], markersize=2, fmt="o", color="C"+str(i), label=r"$\epsilon = $"+str(epsilon[i]))
	plt.plot(kappa_num, D_num[i,:], color="C"+str(i))
plt.legend(loc="best", ncol=3, fontsize=8)
plt.plot(kappa_num+np.linspace(-2, 2, len(kappa_num)), np.ones(len(kappa_num))*D0_ana, "k", label=r"$\epsilon = 0$")
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.axis([0.05, 2.55, 1.65, 2.5])
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root+"figures/comparison_RW_brenner.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))



numeric = np.load("../data_test/tdata_03_03_D01_.npy")
epsilon = np.array([0.0, 0.1, 0.2, 0.3, 0.5])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
D_num[0, :] = D0_ana

for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		D_num[i+1, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)

plt.figure(2)
for i in range(len(kappa)):
	plt.plot(epsilon, D_num[:,i], color="C"+str(i), label=r"$\kappa=$"+str(kappa_num[i]))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
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
		#plt.plot(numeric[i, j,   :, 0], numeric[i, j,   :, 8])
		#plt.plot(numeric[i, j, -T:, 0], numeric[i, j, -T:, 8])
		if i == 2:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -2*T:-T, 8], numeric[i, j, -2*T:-T, 0])/(tau)
		else:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)
	#plt.show()

ana_kappa = np.arange(0.2, 2.201, 0.4)
ana_deff = np.load("../finite_element/data/D_parallels_kappa_D1.npy")

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
#filename = root+"figures/D_eff_vs_eps_D1.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))


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




