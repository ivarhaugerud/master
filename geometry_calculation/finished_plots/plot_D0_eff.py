import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd
import matplotlib
import os

matplotlib.rc('xtick', labelsize=14)
plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

D   = np.logspace(-3, 3, 300)
n   = np.logspace(-3.001, 3.001, 300)
Gamma = np.sqrt(1j/n)
rho   = np.sqrt(1j/D)
Gamma_c = np.conjugate(Gamma)
rho_c   = np.conjugate(rho)
from numpy import *

D_eff = np.zeros((len(D), len(n), 2))
for i in range(len(n)):
    gamma = Gamma[i]
    gamma_c = Gamma_c[i]
    D_eff[:, i, 0] = np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
    amp = tanh(gamma)*tanh(gamma)/(4*gamma*gamma*(rho*rho-gamma*gamma)**2) * ( (1/(rho*tanh(rho)) - 1/(sinh(rho)**2) + 1/(gamma*tanh(gamma)) - 1/(sinh(gamma)**2) - 4*(rho/tanh(rho) - gamma/tanh(gamma))/(rho*rho-gamma*gamma))/2)
    D_eff[:, i, 1] = 2*np.sqrt(np.real(amp)**2 + np.imag(amp)**2)
D_eff *= 18*105/2

fig = plt.figure(1)
x_, y_ = np.meshgrid(D, n)

ax1 = plt.contourf(x_,y_, np.transpose((D_eff[:,:,0])), levels=np.linspace(0, 1.001, 11))
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Geometric factor $g$', fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"Diffusion coefficient $D$", fontsize=8)
plt.ylabel(r"Kinematic viscosity $\nu$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = "figures/D_0_eff.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
