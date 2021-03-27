import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd
import matplotlib
import os

root = "../../../master_latex/results/"
matplotlib.rc('xtick', labelsize=14)
plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

F0 = 3
Pe = 1.0
D   = np.logspace(-3, 3, 300)
n   = np.logspace(-3.001, 3.001, 300)
Gamma = np.sqrt(1j/n)
rho   = np.sqrt(1j/D)
Gamma_c = np.conjugate(Gamma)
rho_c   = np.conjugate(rho)
from numpy import *

D_eff  = np.zeros((len(D), len(n), 2))
D_eff2 = np.zeros((len(D), len(n), 2))

for i in range(len(n)):
    gamma = Gamma[i]
    gamma_c = Gamma_c[i]
    D_eff[:, i, 0] = np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
    D_eff2[:, i, 0] = 1+ F0*F0*Pe*Pe*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))    
    amp = F0*F0*Pe*Pe*tanh(gamma)*tanh(gamma)/(4*gamma*gamma*(rho*rho-gamma*gamma)**2) * ( (1/(rho*tanh(rho)) - 1/(sinh(rho)**2) + 1/(gamma*tanh(gamma)) - 1/(sinh(gamma)**2) - 4*(rho/tanh(rho) - gamma/tanh(gamma))/(rho*rho-gamma*gamma))/2)
    D_eff[:, i, 1] = 2*np.sqrt(np.real(amp)**2 + np.imag(amp)**2)
D_eff *= 18*105/2

"""
fig = plt.figure(1)
x_, y_ = np.meshgrid(D, n)

ax1 = plt.contourf(x_,y_, np.transpose((D_eff[:,:,0])), levels=np.linspace(0, 1.001, 11))
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Geometric factor $g$', fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"Molecular Womserley number $\frac{\omega a^2}{D}$", fontsize=8)
plt.ylabel(r"Womserley number $\frac{\omega a^2}{\nu}$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_0_eff.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
"""

Schmidt = np.logspace(-3.6, 2.6, int(1e4))
D = np.logspace(-3, 2, 6)
Pe = 1/D
rhos = np.sqrt(1j/D)
rhos_c = np.conjugate(rhos)
D_eff2 = np.zeros((len(Schmidt), len(rhos), 2))


for i in range(len(rhos)):
	rho   = rhos[i]
	rho_c = rhos_c[i]
	gamma = rho*np.sqrt(Schmidt)
	gamma_c = rho_c*np.sqrt(Schmidt)

	D_eff2[:, i, 0] = 1+ F0*F0*Pe[i]*Pe[i]*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))    
	amp = F0*F0*Pe[i]*Pe[i]*tanh(gamma)*tanh(gamma)/(4*gamma*gamma*(rho*rho-gamma*gamma)**2) * ( (1/(rho*tanh(rho)) - 1/(sinh(rho)**2) + 1/(gamma*tanh(gamma)) - 1/(sinh(gamma)**2) - 4*(rho/tanh(rho) - gamma/tanh(gamma))/(rho*rho-gamma*gamma))/2)
	D_eff2[:, i, 1] = np.sqrt(np.real(amp)**2 + np.imag(amp)**2)


fig = plt.figure(1)
for i in range(len(rhos)):
	#plt.plot(Schmidt, D_eff2[:,0], color="C"+str(i))
	plt.fill_between(Schmidt, D_eff2[:,i,0]+D_eff2[:,i, 1]-1, D_eff2[:,i,0]-D_eff2[:,i, 1]-1, alpha=0.77, color="C"+str(i), label=r"$\rho^2=10^{%1.0f}$" % np.log10(1/D[i]))

plt.xscale('log')
plt.yscale('log')
plt.legend(loc="lower center", fontsize=8, ncol=3, handlelength=1.4)
plt.xlabel(r"Schmidt number $\frac{\nu}{D}$", fontsize=8)
plt.ylabel(r"Effective Diffusion coefficient $D_\parallel -1 $", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D0_vs_Schmidt.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()