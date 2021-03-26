import os
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import scipy.interpolate as scint

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	
root = "../../../master_latex/results/"


Lx = np.array([125.66, 62.83])#, 41.89])
kappa = 2*np.pi/Lx
epsilon = np.array([0.0178, 0.0316, 0.056, 0.1, 0.177])
base = "data_test/benchmark/"

D          = np.zeros((len(epsilon), len(Lx)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = int(3.0/(0.5*0.003))
tau = 3.0*2
omega = 2*np.pi/tau

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		res = int(100*(1+4*epsilon[i]))
		data = np.loadtxt(base+"Lx"+str(Lx[j])+"_tau3.0_eps"+str(epsilon[i])+"_nu1.2_D1.0_fzero0.0_fone12.0_res"+str(res)+"_dt0.003/tdata.dat")
		D[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[i, j] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
		difference[i, j] = abs(D[i,j] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[i,j]
		plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
	plt.show()


tau = 3
nu  = 1.2
Dm   = 1.0
F0  = 12/nu 
omega = 2*np.pi/tau
Pe = 1/Dm
gamma = np.sqrt(1j*omega/nu)
gamma_c = np.conj(gamma)
rho = np.sqrt(1j*omega/Dm)
rho_c = np.conj(rho)

D0 = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))


#epsilon[0] += 1e-3
D_ana = np.load("finite_element/data/vary_kappa_small.npy")

for i in range(len(kappa)):
	plt.plot(epsilon, difference[:,i])
plt.yscale("log")
plt.xscale("log")
plt.show()
print(np.shape(D0), np.shape(D_ana))

for j in range(len(kappa)):
	plt.plot(epsilon, D0 + epsilon*epsilon*D_ana[j])
	plt.plot(epsilon, D[:,j], "o", markersize=3)
plt.xlabel(r" Womersley number $\frac{\omega a^2}{\nu}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/D_eff_vs_nu.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()



for j in range(len(kappa)):
	plt.plot(epsilon, abs(D[:,j]-(D0 + epsilon*epsilon*D_ana[j]))/D[:,j], "o")
	#plt.plot(epsilon, D[:,j], "o", markersize=3)

epsilon = np.linspace(0.01, 0.177, 1000)
plt.plot(epsilon, epsilon**4, "k")
plt.plot(epsilon, epsilon**2, "r")
plt.xlabel(r" Womersley number $\frac{\omega a^2}{\nu}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.yscale("log")
plt.xscale("log")
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/D_eff_vs_nu.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()