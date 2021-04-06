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


periods = 400
dt = 0.004
t = np.linspace(0, 400*3, int(400*3/dt))
epsilon = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
kappa   = np.array([0.20, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2

tau = 3
dt  = 0.004
nu  = 12
D   = 0.1
F0  = 12/nu 
omega = 2*np.pi/tau
Sc = nu 
Pe = 1/D
gamma = np.sqrt(1j*omega/Sc)
gamma_c = np.conj(gamma)
rho = np.sqrt(1j*omega/D)
rho_c = np.conj(rho)
T = int(0.5*tau/dt)

D0_ana = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))

cutoff = int(len(t)*0.5)
D_RW = np.zeros((len(epsilon), len(kappa), 2))
D_B  = np.zeros((len(epsilon), len(kappa)))

"""
#yhat = savgol_filter(y, 51, 3) # window size 51, polynomial order 3
for i in range(len(epsilon)):
	for j in range(len(kappa)):
		#plt.plot(t[::10], ((RW_sim[i, j, :]))[::10])
		#plt.plot(t[::10], (savgol_filter(RW_sim[i, j, :], 2501, 3))[::10])
		plt.plot(t[::10], np.gradient(savgol_filter(RW_sim[i, j, :], 2501, 3), t)[::10])
	plt.show()
"""
dirr = "data_test/RW_variance/"
B = np.load("data_test/B_solver_12_03_nu12_fone12_D01_dt0.004_res150.npy")
print(np.shape(B))

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		var = np.load(dirr + 'var_over_2Dm_D01_kappa{:.2f}_eps'.format(kappa[j])+str(epsilon[i])[:3] + "_periods400.npy")
		D_RW[i, j, 0]   = np.mean(var[cutoff:]/t[cutoff:])
		D_RW[i, j, 1]   = np.std( var[cutoff:]/t[cutoff:])
		#print(sci.trapz(B[i, j, -T:, 8], B[i, j, -T:, 0])/tau)

		D_B[i, j] = sci.trapz(B[i, j, -T:, 8], B[i, j, -T:, 0])/(tau/2)
		plt.plot(B[i, j, :, 0], B[i, j, :, 8], label="%3.2f" % kappa[j])
		plt.plot(B[i, j, -T:, 0], B[i, j, -T:, 8])


		#plt.plot(t, var)
		#plt.plot(t[cutoff:], var[cutoff:])
	#plt.legend(loc="best")
	#plt.title("%3.2f" % epsilon[i])
plt.show()

plt.plot(kappa, np.ones(len(kappa))*D0_ana 	)
for i in range(len(epsilon)):
	plt.title("%3.2f" % (epsilon[i]))
	plt.errorbar(kappa, D_RW[i, :, 0], yerr=D_RW[i,:,1], fmt="o")
	plt.plot(kappa, D_B[i, :])
plt.show()
"""
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
filename = "figures/comparison_RW_brenner.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))




numeric = np.load("data_test/tdata_03_03_D01_.npy")
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
filename = "figures/D_eff_vs_eps_D01.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
#plt.show()


numeric = np.load("data_test/tdata_04_03_D1_.npy")
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

plt.figure(3)
for i in range(len(epsilon)):
	plt.plot(kappa, D_num[i,:], color="C"+str(i), label=r"$\epsilon=$"+str(epsilon[i]))
plt.legend(loc="best", fontsize=8, ncol=2)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = "figures/D_eff_vs_eps_D1.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


numeric = np.load("data_test/tdata_04_03_D10_.npy")
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
		plt.plot(np.trim_zeros(numeric[i, j,   :, 0]), np.trim_zeros(numeric[i, j,   :, 8]))
		#plt.plot(numeric[i, j, -T:, 0], numeric[i, j, -T:, 8])
		D_num[i+1, j]   = sci.trapz(numeric[i, j, -2*T:-T, 8], numeric[i, j, -2*T:-T, 0])/(tau)
	plt.show()

plt.figure(3)
for i in range(len(epsilon)):
	plt.plot(kappa, D_num[i,:], color="C"+str(i), label=r"$\epsilon=$"+str(epsilon[i]))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = "figures/D_eff_vs_eps_D10.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()




"""