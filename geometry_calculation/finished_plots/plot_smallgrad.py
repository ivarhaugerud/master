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


Lx = np.array([41.89])#, 62.83])#, 41.89])
kappa = 2*np.pi/Lx
epsilon = np.array([0.0, 0.0178, 0.0316, 0.056, 0.1, 0.177, 0.245, 0.31])
base = "../data_test/benchmark_3/"

D          = np.zeros((len(epsilon), len(Lx)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = int(6/(0.003))
tau = 6
omega = 2*np.pi/tau

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		res = int(170)
		data = np.loadtxt(base+"Lx"+str(Lx[j])+"_tau6.0_eps"+str(epsilon[i])+"_nu100.0_D25.0_fzero0.0_fone5000.0_res"+str(res)+"_dt0.003/tdata.dat")
		D[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[i, j] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
		difference[i, j] = abs(D[i,j] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[i,j]
		#plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
		#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
		#plt.show()


tau = 6.0
nu  = 100
Dm   = 25
F0  = 5000/nu 
omega = 2*np.pi/tau
Pe = 1/Dm
gamma = np.sqrt(1j*omega/nu)
gamma_c = np.conj(gamma)
rho = np.sqrt(1j*omega/Dm)
rho_c = np.conj(rho)

D0 = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
from numpy import *
D2  = 1/2*F0**2*Pe**2*(7/12*rho**8*conjugate(rho)**2 + rho**8/2 + rho**8*conjugate(rho)**3/(2*tanh(conjugate(rho))) + 7/4*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 7/4*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1/2*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 1/2*rho**7*conjugate(rho)**4/tanh(rho) + 1/2*rho**7*conjugate(rho)**4/tanh(rho)**3 - 11/4*rho**6*conjugate(rho)**4 - rho**6*conjugate(rho)**2 - rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 9/4*rho**6*conjugate(rho)**4/tanh(rho)**2 + rho**5*conjugate(rho)**6/tanh(rho) - 1/4*rho**5*conjugate(rho)**4/tanh(rho) - rho**5*conjugate(rho)**6/tanh(rho)**3 + 11/4*rho**4*conjugate(rho)**6 + 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 1/4*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 9/4*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4*rho**4*conjugate(rho)**6/tanh(rho)**2 - 1/2*rho**3*conjugate(rho)**8/tanh(rho) + 4*rho**3*conjugate(rho)**6/tanh(rho) + 1/2*rho**3*conjugate(rho)**8/tanh(rho)**3 - 7/12*rho**2*conjugate(rho)**8 + rho**2*conjugate(rho)**6 + 7/4*rho**2*conjugate(rho)**8/tanh(rho)**2 - 7/4*rho*conjugate(rho)**8/tanh(rho) - 1/2*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4)
D2 += -1/(kappa*kappa) + 1/(kappa*tanh(kappa)) + 1 #approx 8/3 + O(kappa^2)
DA = D[0,0] + epsilon*epsilon*D2
D_ana = np.load("../finite_element/data/vary_kappa_small.npy")
"""
for i in range(len(kappa)):
	plt.plot(epsilon, difference[:,i])
plt.yscale("log")
plt.xscale("log")
plt.show()
"""
plt.figure(1)
for j in range(len(kappa)):
	plt.plot(epsilon, D[:,j]-1, "-", label="Numeric")
	plt.plot(epsilon, D[0,0]-1 + epsilon*epsilon*D_ana[j], "o", markersize=3, label="Semi-analytic")
	plt.plot(epsilon, DA-1, "o", markersize=3, label="Analytic")
plt.xscale("log")
plt.yscale("log")

plt.xlabel(r" Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel-1$",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_eff_ana_and_num.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
#plt.show()


plt.figure(2)
for j in range(len(kappa)):
	plt.plot(epsilon, abs(D[:,j]-(D[0,0] + epsilon*epsilon*D_ana[j]))/D[:,j], "o", markersize=3, label="Semi-analytic")
	#plt.plot(epsilon, abs(D[:,j]-DA), "o", markersize=3, label="Analytic")

epsilon = np.linspace(0.01, 0.35, 1000)
gamma = np.real(gamma)
plt.plot(epsilon, epsilon**4/(1-epsilon), label=r"$\epsilon^4$")
#plt.plot(epsilon, epsilon*epsilon*(gamma*gamma+kappa*kappa), label=r"$\epsilon^2(\gamma^2+\kappa^2)$")
plt.xlabel(r" Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r" Relative difference $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.yscale("log")
plt.xscale("log")
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/rel_diff_semi_analytic.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
#plt.show()



plt.figure(3)
epsilon = np.array([0.0, 0.0178, 0.0316, 0.056, 0.1, 0.177, 0.245, 0.31])
for j in range(len(kappa)):
	#plt.plot(epsilon, abs(D[:,j]-(D[0,0] + epsilon*epsilon*D_ana[j]))/D[:,j], "o", markersize=3, label="Semi-analytic")
	plt.plot(epsilon, abs(D[0,0]+epsilon*epsilon*D_ana[j]-DA), "o", markersize=3, label="Analytic")
epsilon = np.linspace(0.01, 0.35, 1000)
gamma = np.real(gamma)
#plt.plot(epsilon, epsilon**4/(1-epsilon), label=r"$\epsilon^4$")
plt.plot(epsilon, epsilon*epsilon*(gamma*gamma+kappa*kappa), label=r"$\epsilon^2(\gamma^2+\kappa^2)$")
plt.plot(epsilon, 5*epsilon*epsilon*(gamma*gamma+kappa*kappa), label=r"$5\epsilon^2(\gamma^2+\kappa^2)$")

plt.xlabel(r" Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r" Difference $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.yscale("log")
plt.xscale("log")
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.axis([0.015, 0.40, 5*1e-6, 0.03])
filename = root + "figures/rel_diff_analytic.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()