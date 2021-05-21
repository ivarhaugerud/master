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


Lx = np.array([17.95, 31.41, 41.89])
Lx2 = np.array([17.95, 41.88])
kappa  = 2*np.pi/Lx
kappa2 = 2*np.pi/Lx2
epsilon = np.array([0.0, 0.0178, 0.0316, 0.056, 0.1, 0.177, 0.245, 0.31])
base = "../data_test/benchmark_3/"

D          = np.zeros((len(epsilon), len(Lx)))
D2         = np.zeros((len(epsilon), len(Lx2)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = int(6/(0.003))
tau = 6
omega = 2*np.pi/tau
rho  = np.sqrt(2*np.pi/(tau*25))
rho2 = np.sqrt(2*np.pi/(tau*4))

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		res = int(170)
		data = np.loadtxt(base+"Lx"+str(Lx[j])+"_tau6.0_eps"+str(epsilon[i])+"_nu5.0_D25.0_fzero0.0_fone250.0_res"+str(res)+"_dt0.003/tdata.dat")
		D[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[i, j] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
		difference[i, j] = abs(D[i,j] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[i,j]

for i in range(len(epsilon)):
	for j in range(len(kappa2)):
		data = np.loadtxt(base+"Lx"+str(Lx2[j])+"_tau6.0_eps"+str(epsilon[i])+"_nu5.0_D4.0_fzero0.0_fone250.0_res"+str(res)+"_dt0.003/tdata.dat")
		D2[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau


tau = 6.0
nu  = 5
Dm   = 4
F0  = 250/nu 
omega = 2*np.pi/tau
Pe = 1/Dm
gamma = np.sqrt(1j*omega/nu)
gamma_c = np.conj(gamma)
rho = np.sqrt(1j*omega/Dm)
rho_c = np.conj(rho)

from numpy import *
D2_b  = 1/2*F0**2*Pe**2*(7/12*rho**8*conjugate(rho)**2 + rho**8/2 + rho**8*conjugate(rho)**3/(2*tanh(conjugate(rho))) + 7/4*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 7/4*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1/2*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 1/2*rho**7*conjugate(rho)**4/tanh(rho) + 1/2*rho**7*conjugate(rho)**4/tanh(rho)**3 - 11/4*rho**6*conjugate(rho)**4 - rho**6*conjugate(rho)**2 - rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 9/4*rho**6*conjugate(rho)**4/tanh(rho)**2 + rho**5*conjugate(rho)**6/tanh(rho) - 1/4*rho**5*conjugate(rho)**4/tanh(rho) - rho**5*conjugate(rho)**6/tanh(rho)**3 + 11/4*rho**4*conjugate(rho)**6 + 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 1/4*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 9/4*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4*rho**4*conjugate(rho)**6/tanh(rho)**2 - 1/2*rho**3*conjugate(rho)**8/tanh(rho) + 4*rho**3*conjugate(rho)**6/tanh(rho) + 1/2*rho**3*conjugate(rho)**8/tanh(rho)**3 - 7/12*rho**2*conjugate(rho)**8 + rho**2*conjugate(rho)**6 + 7/4*rho**2*conjugate(rho)**8/tanh(rho)**2 - 7/4*rho*conjugate(rho)**8/tanh(rho) - 1/2*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4) + (-1/(kappa2*kappa2) + 1/(kappa2*tanh(kappa2)) + 1) #approx 8/3 + O(kappa^2)

Dm   = 25
Pe = 1/Dm
rho = np.sqrt(1j*omega/Dm)
rho_c = np.conj(rho)

from numpy import *
D2_a  = 1/2*F0**2*Pe**2*(7/12*rho**8*conjugate(rho)**2 + rho**8/2 + rho**8*conjugate(rho)**3/(2*tanh(conjugate(rho))) + 7/4*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 7/4*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1/2*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 1/2*rho**7*conjugate(rho)**4/tanh(rho) + 1/2*rho**7*conjugate(rho)**4/tanh(rho)**3 - 11/4*rho**6*conjugate(rho)**4 - rho**6*conjugate(rho)**2 - rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 9/4*rho**6*conjugate(rho)**4/tanh(rho)**2 + rho**5*conjugate(rho)**6/tanh(rho) - 1/4*rho**5*conjugate(rho)**4/tanh(rho) - rho**5*conjugate(rho)**6/tanh(rho)**3 + 11/4*rho**4*conjugate(rho)**6 + 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 1/4*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 9/4*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4*rho**4*conjugate(rho)**6/tanh(rho)**2 - 1/2*rho**3*conjugate(rho)**8/tanh(rho) + 4*rho**3*conjugate(rho)**6/tanh(rho) + 1/2*rho**3*conjugate(rho)**8/tanh(rho)**3 - 7/12*rho**2*conjugate(rho)**8 + rho**2*conjugate(rho)**6 + 7/4*rho**2*conjugate(rho)**8/tanh(rho)**2 - 7/4*rho*conjugate(rho)**8/tanh(rho) - 1/2*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4) + (-1/(kappa*kappa) + 1/(kappa*tanh(kappa)) + 1 )#approx 8/3 + O(kappa^2)



plt.figure(1)
for j in range(len(kappa)):
	plt.plot(epsilon, abs(D[:,j]-D[0,j]-epsilon*epsilon*D2_a[j]), "o", color="C"+str(j), markersize=4, label=r"$\kappa=%3.2f$" % (kappa[j]))
	plt.plot(epsilon, epsilon**2*(kappa[j]*kappa[j]+np.real(gamma*np.conjugate(gamma))),  color="C"+str(j))

for j in range(len(kappa2)):
	plt.plot(epsilon, abs(D2[:,j]-D2[0,j]-epsilon*epsilon*D2_b[j]), "+", color="C"+str(j+len(kappa)), markersize=4, label=r"$\kappa=%3.2f$" % (kappa[j]))
	plt.plot(epsilon, epsilon**2*(kappa2[j]*kappa2[j]+np.real(gamma*np.conjugate(gamma))),  color="C"+str(j+len(kappa)))

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r" Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r" Second order effective Diffusion Coefficient $ D_\parallel^{(2)}$",  fontsize=8)
plt.legend(loc="best", fontsize=8, ncol=1)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_eff_analytic_and_num.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


"""
plt.figure(1)
for j in range(len(kappa)):
	plt.plot(epsilon, D[:,j]-D[0,j], "-", color="C"+str(j))
	plt.plot(epsilon, epsilon*epsilon*D2_a[j], "o", color="C"+str(j), markersize=4, label=r"$\kappa=%3.2f$" % (kappa[j]))
plt.xscale("log")
plt.yscale("log")
plt.show()
"""

D_ana = np.load("../finite_element/data/D_eff_vary_D_nu5_D25_tau6_F0250.npy")
D_ana2 = np.load("../finite_element/data/D_eff_vary_D_nu5_D4_tau6_F0250.npy")
D_ana = np.array([D_ana[0], D_ana[2], D_ana[3]])



plt.figure(1)
for j in range(len(kappa)):
	plt.plot(epsilon, D[:,j]-D[0,j], "-", color="C"+str(j))
	plt.plot(epsilon, epsilon*epsilon*D_ana[j], "o", color="C"+str(j), markersize=4, label=r"$\kappa=%3.2f$" % (kappa[j]))

for j in range(len(kappa2)):
	plt.plot(epsilon, D2[:,j]-D2[0,j], "-", color="C"+str(j+len(kappa)))
	plt.plot(epsilon,  epsilon*epsilon*D_ana2[j], "x", color="C"+str(j+len(kappa)), markersize=4, label=r"$\kappa=%3.2f$" % (kappa[j]))
plt.xscale("log")
plt.yscale("log")

plt.xlabel(r" Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r" Second order effective dispersion $ D_\parallel^{(2)}$ [$D_m$]",  fontsize=8)
plt.legend(loc="best", fontsize=8, ncol=1)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_eff_ana_and_num.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


plt.figure(2)
for j in range(len(kappa)):
	plt.plot(epsilon, abs(D[:,j]-(D[0,j] + epsilon*epsilon*D_ana[j])), "o", markersize=4, color="C"+str(j), label=r"$\kappa=%3.2f$" % (kappa[j]))
	plt.plot(epsilon, abs(D[:,j]-(D[0,j] + epsilon*epsilon*D_ana[j])), color="C"+str(j), linewidth=0.5)


for j in range(len(kappa2)):
	plt.plot(epsilon, abs(D2[:,j]-(D2[0,j] + epsilon*epsilon*D_ana2[j])), "x", markersize=4, color="C"+str(j+len(kappa)), label=r"$\kappa=%3.2f$" % (kappa[j]))
	plt.plot(epsilon, abs(D2[:,j]-(D2[0,j] + epsilon*epsilon*D_ana2[j])), color="C"+str(j+len(kappa)), linewidth=0.5)

epsilon = np.linspace(0.01, 0.45, 1000)
plt.plot(epsilon, epsilon**4/(1-epsilon), "k", label=r"$\epsilon^4$")

#plt.plot(epsilon, epsilon*epsilon*(gamma*gamma+kappa*kappa), label=r"$\epsilon^2(\gamma^2+\kappa^2)$")
plt.xlabel(r" Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r" Difference in $D_\parallel$ [$D_m$]",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.yscale("log")
plt.xscale("log")
plt.axis([0.01, 0.413, 5*1e-8, 0.05])
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/rel_diff_semi_analytic.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
