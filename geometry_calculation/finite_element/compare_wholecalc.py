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







ana_kappa = np.arange(0.2, 2.201, 0.4)
ana_deff = np.load("../finite_element/data/D_parallels_kappa_D1.npy")

numeric = np.load("../data_test/tdata_04_03_D1_.npy")
epsilon = np.array([0.0, 0.1, 0.2, 0.3])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
difference = np.zeros(np.shape(D_num))

D = 1.0
Pe = 1/D
rho = np.sqrt(1j*omega/D)
rho_c = np.conj(rho)
D_num[0, :] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))

for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		#plt.plot(numeric[i, j,   :, 0], numeric[i, j,   :, 8	])
		#plt.plot(numeric[i, j, -T:, 0], numeric[i, j, -T:, 8	])
		if i == 2:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -2*T:-T, 8	], numeric[i, j, -2*T:-T, 0])/(tau)
			difference[i+1, j]   = sci.trapz(numeric[i, j, -3*T:-2*T, 8	], numeric[i, j, -3*T:-2*T, 0])/(tau)
		else:
			D_num[i+1, j]   = sci.trapz(numeric[i, j, -T:, 8	], numeric[i, j, -T:, 0])/(tau)
			difference[i+1, j]   = sci.trapz(numeric[i, j, -2*T:-1*T, 8	], numeric[i, j, -2*T:-1*T, 0])/(tau)
	#plt.show()

for i in range(len(epsilon)):
	plt.plot(kappa, abs(D_num[i,:]-difference[i,:])/D_num[i,:])
plt.show()
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
plt.show()
#filename = root+"figures/D_eff_vs_eps_D1.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
