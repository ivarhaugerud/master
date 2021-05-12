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

t = np.linspace(0, 300*3, len(RW_sim[0,0,:]))
epsilon = np.array([0.1, 0.2, 0.3])

kappa       = np.array([0.2, 0.6, 1.0, 1.4, 1.7, 2.1]) #0.2
kappa_num   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.5]) #0.2

eps2 = np.arange(0.4, 0.501, 0.1)
additional_RW = np.zeros((len(eps2), len(kappa_num), 2))

for i in range(len(eps2)):
	for j in range(len(kappa)):
		data = np.load("../data_test/RW_eps_04_05/var_over_2Dm_D01_kappa%3.2f" %(kappa_num[j])+"_eps"+str(eps2[i])+"_periods400.npy")
		t2 = np.linspace(0, 600*3, len(data))
		cutoff = int(len(t2)/2)
		additional_RW[i, j, 0] = np.mean(data[cutoff:]/t2[cutoff:])
		additional_RW[i, j, 1] = np.std( data[cutoff:]/t2[cutoff:])

eps3 = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])
T3 = int(3.0/0.004)
Lx3 = np.array([31.41, 10.47, 6.283, 4.487, 3.49, 2.855])
kappa3 = 2*np.pi/Lx3
bren_D01 = np.zeros((len(eps3), len(kappa3)))
bren_D01_amp = np.zeros((len(eps3), len(kappa3)))
difference = np.zeros(np.shape(bren_D01))
U_bren_D01 = np.zeros(np.shape(bren_D01))

for i in range(len(eps3)):
	for j in range(len(kappa3)):
		data = np.loadtxt("../data_test/tdatas_tau3.0_nu1.2_D0.1_fone12.0/Lx"+str(Lx3[j])+"_tau3.0_eps"+str(eps3[i])+"_nu1.2_D0.1_fzero0.0_fone12.0_res150_dt0.004/tdata.dat")
		bren_D01[i,j] = sci.trapz(data[-T3:, 8], data[-T3:, 0])/3.0
		difference[i, j] = abs(bren_D01[i,j]-sci.trapz(data[-2*T3:-T3, 8], data[-2*T3:-T3, 0])/3.0)/(bren_D01[i,j])
		U_bren_D01[i,j]   = sci.trapz(data[-T3:, 4], data[-T3:, 0])/3.0
		bren_D01_amp[i,j] = (np.max( abs(data[-T3:, 8] - sci.trapz(data[-T3:, 8], data[-T3:, 0])/3.0 )))/bren_D01[i,j]
		#plt.plot(data[:, 0], data[:, 8])
		#plt.plot(data[-T3:, 0], data[-T3:, 8])
	#plt.show()

plt.figure(4)
for i in range(len(eps3)):
	plt.plot(kappa3, U_bren_D01[i,:])
#plt.yscale("log")
plt.show()

plt.figure(4)
for i in range(len(eps3)):
	plt.plot(kappa3, difference[i,:])
plt.yscale("log")
#plt.show()

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
RW_sim2  = np.zeros((len(epsilon), len(kappa), 2))

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		D_RW[i, j, 0] = np.mean(RW_sim[i, j, cutoff:]/t[cutoff:])
		D_RW[i, j, 1] = np.std( RW_sim[i, j, cutoff:]/t[cutoff:])

plt.figure(1)
for i in range(len(epsilon)):
	for j in range(len(kappa)):
		RW_sim2[i, j, 0] = (np.mean(RW_sim_old[i, j, cutoff:]) + D_RW[i, j, 0])/2
		RW_sim2[i, j, 1] = np.sqrt( np.std( RW_sim_old[i, j, cutoff:])**2 + D_RW[i, j, 1]**2 )/np.sqrt(2)
		t2 = np.linspace(0, 1, len(RW_sim_old[i, j,:]))
		if j == 3:
			plt.plot(t2, RW_sim_old[i, j,:], label=r"$\epsilon=%2.1f$" % epsilon[i])

for i in range(len(eps2)):
	for j in range(len(kappa)):
		data = np.load("../data_test/RW_eps_04_05/var_over_2Dm_D01_kappa%3.2f" %(kappa_num[j])+"_eps"+str(eps2[i])+"_periods400.npy")
		t2 = np.linspace(0, 600*3, len(data))
		if j == 3:
			plt.plot(t2/(600*3), data/t2, label=r"$\epsilon=%2.1f$" % eps2[i])
plt.legend(loc="best", ncol=3, fontsize=8)
plt.xlabel(r"Time $[T_{max}]$", fontsize=8)
plt.axis([-0.02, 1.02, 0.8, 2.5])
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$ [$D_m$]", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/RW_vs_time_oscflow.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
plt.figure(1)
for i in range(len(epsilon)):
	plt.errorbar(kappa, RW_sim2[i, :, 0], yerr=RW_sim2[i, :, 1], markersize=2, fmt="o", color="C"+str(i), label=r"$\epsilon = $"+str(epsilon[i]))
	#plt.plot(kappa_num, D_num[i,:], color="C"+str(i))
for i in range(len(eps2)):
	plt.errorbar(kappa_num, additional_RW[i, :, 0], yerr=additional_RW[i, :, 1], markersize=2, fmt="o", color="C"+str(i+3), label=r"$\epsilon = $"+str(eps2[i]))

counter = 0
for i in range(len(eps3)):
	if (i + 1)%2 == 0:
		plt.plot(kappa3, bren_D01[i,:], color="C"+str(counter))
		counter += 1
plt.plot(kappa_num+np.linspace(-2, 2, len(kappa_num)), np.ones(len(kappa_num))*D0_ana, "k", label=r"$\epsilon = 0$")
plt.legend(loc="upper center", ncol=3, fontsize=8)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.axis([0.05, 2.35, 1.15, 2.76])
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$ [$D_m$]", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/comparison_RW_brenner.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))

plt.show()
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


numeric = np.load("../data_test/tdata_03_03_D01_.npy")
epsilon = np.array([0, 0.1, 0.2, 0.3, 0.5])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
D_nuA = np.zeros((len(epsilon), len(kappa)))
U     = np.zeros((len(epsilon), len(kappa)))
T =  750
tau = 3
omega = 2*np.pi/tau
F0 = 12

D_num[0, :] = np.real(D0_ana)

for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		D_num[i+1, j]   = sci.trapz(np.trim_zeros(numeric[i, j, :, 8])[-T:], np.trim_zeros(numeric[i, j, :, 0])[-T:])/(tau)
		U[i+1, j]       = np.sqrt(sci.trapz(numeric[i, j, -T:, 8], numeric[i, j, -T:, 0])/(tau))
		D_nuA[i+1, j]   = np.max(abs(numeric[i,j, -T:, 8] - D_num[i+1, j])/D_num[i+1, j])


plt.figure(1)
for i in range(len(kappa)):
	plt.plot(epsilon[1:], D_num[1:,i]/(1+2/105*(U[1:,i]/0.1)**2), color="C"+str(i), label=r"$\kappa=$"+str(kappa_num[i]))
	print("Pe = ", U[1:,i]/0.1)

plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Relative Effective Dispersion $D_\parallel/D_\parallel^{aris}$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.axis([0.08, 0.518, 0.36, 0.42999])
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/rel_D_eff_vs_eps_D01.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))

plt.figure(2)
for i in range(len(kappa)):
	plt.plot(epsilon[1:], D_num[1:,i], color="C"+str(i), label=r"$\kappa=$"+str(kappa_num[i]))

plt.legend(loc="lower left", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Effective Dispersion $D_\parallel$ [$D_m$]", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/D_eff_vs_eps_D01.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

numeric = np.load("../data_test/tdata_04_03_D1_.npy")
epsilon = np.array([0.0, 0.1, 0.2, 0.3])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2
D_num = np.zeros((len(epsilon), len(kappa)))
D_nuA = np.zeros((len(epsilon), len(kappa)))
U     = np.zeros((len(epsilon), len(kappa)))

D = 1.0
Pe = 1/D
rho = np.sqrt(1j*omega/D)
rho_c = np.conj(rho)
xi = np.linspace(-1, 1, int(1e4))
D0_ana = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
u0 = 0.5*sci.trapz( (F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma))*np.conjugate(F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)) , xi)
U[0, :] = np.sqrt(u0)
D_num[0, :] = np.real(D0_ana)
for i in range(len(epsilon)-1):
	for j in range(len(kappa)):
		D_num[i+1, j]   = sci.trapz( np.trim_zeros(numeric[i, j, :, 8])[-T:], np.trim_zeros(numeric[i, j, :, 0])[-T:])/(tau)
		U[i+1, j]       = np.sqrt(sci.trapz( np.trim_zeros(numeric[i, j, :, 4])[-T:], np.trim_zeros(numeric[i, j, :, 0])[-T:])/(tau))
		print("Pe = ", U[i+1,j]/1)


plt.figure(3)
for i in range(len(epsilon)-1):
	plt.plot(kappa, D_num[i+1,:]/(1+2/105*(U[i+1,:]/1)**2), color="C"+str(i+1), label=r"$\epsilon=$"+str(epsilon[i+1]))
	plt.plot(kappa, D_num[i+1,:]/(1+2/105*(U[i+1,:]/1)**2), "o", markersize=3, color="C"+str(i+1))

plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Relative Effective Dispersion $D_\parallel/D_\parallel^{aris}$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/rel_D_eff_vs_eps_D1.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))


plt.figure(30)
for i in range(len(epsilon)):
	plt.plot(kappa, D_num[i,:], color="C"+str(i), label=r"$\epsilon=$"+str(epsilon[i]))
	plt.plot(kappa, D_num[i,:], "o", markersize=3, color="C"+str(i))


plt.legend(loc="best", fontsize=8, ncol=2)
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Effective Dispersion $D_\parallel$ [$D_m$]", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/D_eff_vs_kappa_D1.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


numeric = np.load("../data_test/tdata_04_03_D10_.npy")
epsilon = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
kappa   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])
D_num = np.zeros((len(epsilon), len(kappa)))
D_nuA = np.zeros((len(epsilon), len(kappa)))
U     = np.zeros((len(epsilon), len(kappa)))

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
		D_num[i+1, j]   = sci.trapz( np.trim_zeros(numeric[i, j, :, 8])[-T:], np.trim_zeros(numeric[i, j, :, 0])[-T:])/(tau)
		U[i+1, j]       = np.sqrt(sci.trapz( np.trim_zeros(numeric[i, j, :, 4])[-T:], np.trim_zeros(numeric[i, j, :, 0])[-T:])/(tau))
		plt.plot(np.trim_zeros(numeric[i, j, :, 0]), np.trim_zeros(numeric[i, j, :, 4]))
	plt.show()


plt.figure(4)
for i in range(len(kappa)):
	plt.plot(epsilon, D_num[:,i]/(1+2/105*(U[i,:]/10)**2), color="C"+str(i), label=r"$\kappa=$"+str(kappa[i]))
	print((U[i,:]/10))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Relative Effective Dispersion $D_\parallel/D_\parallel^{aris}$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/rel_D_eff_vs_eps_D10.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))




plt.figure(5)
for i in range(len(kappa)):
	plt.plot(epsilon, D_num[:,i], color="C"+str(i), label=r"$\kappa=$"+str(kappa[i]))
	print("PE=", (U[i,:]/10))
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Effective Dispersion $D_\parallel$ [$D_m$]", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root+"figures/D_eff_vs_eps_D10.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()