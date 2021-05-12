import os
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
import scipy.integrate as scpi 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	
root = "../../../master_latex/results/"
Map = matplotlib.cm.get_cmap('Spectral_r')

tau = 3.0
omega    = 2*np.pi/tau

D_D  = np.load("../finite_element/data/res_vs_D_kappa_2.npy")
Ds       = np.arange(0.6, 2.7, 0.1)
kappas   = np.arange(0.2, 1.75, 0.1)
fig = plt.figure(3)
x_, y_ = np.meshgrid(kappas, np.sqrt(omega/Ds))
#x_, y_ = np.meshgrid(kappas, Ds)

ax1 = plt.contourf(x_,y_, np.transpose(D_D), cmap=Map, levels=18)
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Second order effective dispersion $D_\parallel^{(2)}$ [$D_m$]', fontsize=8)
plt.ylabel(r"Diffusive Womserley number $\rho$",    fontsize=8)
plt.xlabel(r"Wave number $\kappa$",         fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/semi_ana_kappa_vs_rho.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))

D_nu = np.load("../finite_element/data/res_vs_nu_kappa.npy")
kappas   = np.arange(0.8, 1.6, 0.1)
nus      = np.arange(0.8, 1.6, 0.1)

fig = plt.figure(4)
x_, y_ = np.meshgrid(kappas, np.sqrt(omega/nus))
ax1 = plt.contourf(x_,y_, D_nu, cmap=Map, levels=18)
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Second order effective dispersion $D_\parallel^{(2)}$ [$D_m$]', fontsize=8)
plt.ylabel(r"Womserley number $\gamma$",    fontsize=8)
plt.xlabel(r"Wave number $\kappa$",       fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/semi_ana_kappa_vs_gamma.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()