import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

RW_sim  = np.load("data_test/final_run_tdata_analytic_D.npy")
numeric = np.load("data_test/final_run_tdata_D01_res200_dt0.003.npy")
print(np.shape(RW_sim))
#numeric = np.load("data_test/tdatas_large_run_2.npy")

print(numeric)
t = np.linspace(0, 1000, len(RW_sim[0,0,:]))

epsilon = np.array([0.1, 0.2])
eps_num = np.array([0.1, 0.2, 0.3, 0.5])

kappa       = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
kappa_num   = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1]) #0.2

tau = 3
dt  = 0.003
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
#D_RW2  = np.zeros((len(epsilon), len(kappa), 2))

D_num = np.zeros((len(eps_num), len(kappa)))

for i in range(len(eps_num)):
	for j in range(len(kappa_num)):
		plt.plot(numeric[i, j, :-T, 0]/tau, numeric[i, j, :-T, 8], label=r"$\kappa=$"+str(kappa_num[j]))
		plt.plot(numeric[i, j, -2*T:-T, 0]/tau, numeric[i, j, -2*T:-T, 8])
		D_num[i, j]   = sci.trapz(numeric[i, j, -2*T:-T, 8], numeric[i, j, -2*T:-T, 0])/(tau)

	plt.title(r"$\epsilon$ = " + str(eps_num[i]))
	plt.legend(loc="best")
	plt.xlabel("Time [periods]")
	plt.ylabel("Effective diffusion coefficient")

	plt.show()

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		#plt.plot(t[::100], RW_sim[i, j, ::100])
		D_RW[i, j, 0] = np.mean(RW_sim[i, j, cutoff:])
		D_RW[i, j, 1] = np.std( RW_sim[i, j, cutoff:])
		#D_RW2[i, j, 0] = np.mean(RW_sim2[i, j, cutoff:])
		#D_RW2[i, j, 1] = np.std( RW_sim2[i, j, cutoff:])
		#print(D_num[i,j])
	#print("\n")

print(D0_ana)

plt.plot(kappa, np.ones(len(kappa))*D0_ana, label=r"$\epsilon = 0$")

for i in range(len(epsilon)):
	plt.errorbar(kappa, D_RW[i, :, 0], yerr=D_RW[i, :, 1], fmt="o", color="C"+str(i+1), label=r"$\epsilon = $"+str(epsilon[i]))
for i in range(len(eps_num)):
	plt.plot(kappa_num, D_num[i,:], color="C"+str(i+1), label=r"$\epsilon=$"+str(eps_num[i]))
plt.legend(loc="best")
plt.show()