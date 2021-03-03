import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
from scipy.signal import savgol_filter

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

RW_sim_old  = np.load("data_test/final_run_RW_D0.1.npy")
RW_sim  = np.load("data_test/RW_pos_03_03__D01.npy")

numeric = np.load("data_test/final_run_tdata_fixed_D01_res150_dt0.004.npy")
numeric = np.load("data_test/tdata_03_03_D01_.npy")

print(np.shape(RW_sim))
#numeric = np.load("data_test/tdatas_large_run_2.npy")

print(numeric)
t = np.linspace(0, 300*3, len(RW_sim[0,0,:]))

EPS = np.array([0.1, 0.2, 0.3])
epsilon = np.array([0.1, 0.2, 0.3, 0.5])
eps_num = epsilon #np.arange(0.05, 0.51, 0.05)

kappa       = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
kappa_num   = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1]) #0.2

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

D_num = np.zeros((len(eps_num), len(kappa)))
"""
#yhat = savgol_filter(y, 51, 3) # window size 51, polynomial order 3
for i in range(len(epsilon)):
	for j in range(len(kappa)):
		#plt.plot(t[::10], ((RW_sim[i, j, :]))[::10])
		#plt.plot(t[::10], (savgol_filter(RW_sim[i, j, :], 2501, 3))[::10])
		plt.plot(t[::10], np.gradient(savgol_filter(RW_sim[i, j, :], 2501, 3), t)[::10])
	plt.show()
"""
for i in range(len(eps_num)):
	for j in range(len(kappa_num)):
		plt.plot(numeric[i, j, :-T, 0]/tau, numeric[i, j, :-T, 8], label=r"$\kappa=$"+str(kappa_num[j]))
		plt.plot(numeric[i, j, -2*T:-T, 0]/tau, numeric[i, j, -2*T:-T, 8])
		D_num[i, j]   = sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau)

	plt.title(r"$\epsilon$ = " + str(eps_num[i]))
	plt.legend(loc="best")
	plt.xlabel("Time [periods]")
	plt.ylabel("Effective diffusion coefficient")

plt.show()

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		plt.plot(t[::10], RW_sim[i, j, ::10])
		D_RW[i, j, 0] = np.mean(RW_sim[i, j, cutoff:]/t[cutoff:])
		D_RW[i, j, 1] = np.std( RW_sim[i, j, cutoff:]/t[cutoff:])

for i in range(len(EPS)):
	for j in range(len(kappa)):
		RW_sim2[i, j, 0] = (np.mean(RW_sim_old[i, j, cutoff:]) + D_RW[i, j, 0])/2
		RW_sim2[i, j, 1] = (np.std( RW_sim_old[i, j, cutoff:]) + D_RW[i, j, 1])/np.sqrt(2)
		#print(D_num[i,j])
	#print("\n")
plt.show()

print(D0_ana)

for i in range(len(epsilon)-1):
	plt.errorbar(kappa, RW_sim2[i, :, 0], yerr=RW_sim2[i, :, 1], fmt="o", color="C"+str(i), label=r"$\epsilon = $"+str(epsilon[i]))
	#plt.errorbar(kappa, RW_sim2[i, :, 0], yerr=RW_sim2[i, :, 1], fmt="o", color="k")
	plt.plot(kappa_num, D_num[i,:], color="C"+str(i), label=r"$\epsilon=$"+str(eps_num[i]))
plt.plot(kappa, np.ones(len(kappa))*D0_ana, "k", label=r"$\epsilon = 0$")
plt.legend(loc="best", ncol=2)
plt.xlabel(r"Wave number $\kappa$")
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$")
plt.show()



for i in range(len(kappa_num)):
	plt.plot(eps_num, D_num[:,i], color="C"+str(i), label=r"$\kappa=$"+str(kappa_num[i]))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.savefig("figures/D_eff_vs_eps.pdf")
plt.show()

"""
for i in range(len(kappa_num)):
	plt.plot(eps_num, D_num[:,i], color="C"+str(i), label=r"$\kappa=$"+str(kappa_num[i]))
	plt.errorbar(epsilon, D_RW[:, i, 0], yerr=D_RW[:, i, 1], fmt="o", color="C"+str(i), label=r"$\epsilon = $"+str(kappa_num[i]))
plt.legend(loc="best")
plt.xlabel(r"Boundary amplitude $\epsilon$")
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$")
plt.show()
"""