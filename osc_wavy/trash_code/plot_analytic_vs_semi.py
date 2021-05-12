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

kappas = np.logspace(-3, -0.5, 6)
D_ana = np.zeros(len(kappas))

from numpy import *

for i, kappa in enumerate(kappas):
	print(i, kappa)
	D_ana[i]  = 1/2*F0**2*Pe**2*(7/12*rho**8*conjugate(rho)**2 + rho**8/2 + rho**8*conjugate(rho)**3/(2*tanh(conjugate(rho))) + 7/4*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 7/4*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1/2*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 1/2*rho**7*conjugate(rho)**4/tanh(rho) + 1/2*rho**7*conjugate(rho)**4/tanh(rho)**3 - 11/4*rho**6*conjugate(rho)**4 - rho**6*conjugate(rho)**2 - rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 9/4*rho**6*conjugate(rho)**4/tanh(rho)**2 + rho**5*conjugate(rho)**6/tanh(rho) - 1/4*rho**5*conjugate(rho)**4/tanh(rho) - rho**5*conjugate(rho)**6/tanh(rho)**3 + 11/4*rho**4*conjugate(rho)**6 + 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 1/4*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 9/4*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4*rho**4*conjugate(rho)**6/tanh(rho)**2 - 1/2*rho**3*conjugate(rho)**8/tanh(rho) + 4*rho**3*conjugate(rho)**6/tanh(rho) + 1/2*rho**3*conjugate(rho)**8/tanh(rho)**3 - 7/12*rho**2*conjugate(rho)**8 + rho**2*conjugate(rho)**6 + 7/4*rho**2*conjugate(rho)**8/tanh(rho)**2 - 7/4*rho*conjugate(rho)**8/tanh(rho) - 1/2*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4)
	D_ana[i] += -1/(kappa*kappa) + 1/(kappa*tanh(kappa)) + 1 #approx 8/3 + O(kappa^2)

D_num = np.load("finite_element/data/benchmark_kappa.npy")


plt.plot(kappas, D_ana, "o", markersize=3, label="Analytic")
plt.plot(kappas, D_num, "o", markersize=3, label="Numerical")
plt.xscale("log")
plt.legend(loc="best", fontsize=8)
plt.show()


plt.plot(kappas, D_ana-D_num, "o", markersize=3)
kappas = np.logspace(-3, -0.5, 1000)
plt.plot(kappas, kappas*kappas, label=r"$\kappa^2$")
plt.xscale("log")
plt.show()



"""

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
"""