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

tau = 3.0
omega    = 2*np.pi/tau
Ds       = np.arange(0.8, 1.6, 0.1)
kappas   = np.arange(0.8, 1.6, 0.1)
nus      = np.arange(0.8, 1.6, 0.1)


D_D  = np.load("finite_element/data/res_vs_D_kappa.npy")
D_nu = np.load("finite_element/data/res_vs_nu_kappa.npy")

for i in range(len(nus)):
	plt.figure(1)
	plt.plot(kappas, D_nu[i, :])
	plt.figure(2)
	plt.plot(kappas, D_D[i, :])

plt.xlabel(r" Wave number number $\kappa$", fontsize=8)
plt.ylabel(r" Second order Effective Diffusion coefficient $ D_\parallel^{(2)}$",  fontsize=8)
#handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1], loc="best", fontsize=8, ncol=1, markerscale=0.1, title=r"$\rho^2$", columnspacing=1.0, labelspacing=0.8)
#plt.axis([0.45, 1.8, 0.42, 1.3])
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/analytic_resonance.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
