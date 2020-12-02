import numpy as np 
import matplotlib.pyplot as plt 
import os

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'

D_parallels = np.load("data/D_parallels_kappa.npy")


kappas = np.arange(0.25, 2.25+1e-3, 0.25)
Sc = 1.2
omega = 1
F0 = 3
Pe = F0*Sc

filename = "figures/data_para_vs_kappa.pdf"
plt.scatter(kappas, D_parallels, s=4)
plt.plot(kappas, D_parallels, "-", linewidth=1)
plt.plot([np.sqrt(2*omega/Sc), np.sqrt(2*omega/Sc)], [-30, 30])
tol = 0.2*min(kappas)
plt.axis([min(kappas)-tol, max(kappas)+tol, min(D_parallels)*0.8, max(D_parallels)*1.2])
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Parallel Diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()