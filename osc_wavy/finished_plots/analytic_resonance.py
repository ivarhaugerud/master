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

nu = 1.2
D = 1.0
F0 = 12/nu 
Pe = 1/D
tau = np.arange(2.0, 4.0, 0.25)
kappas = np.arange(0.5, 1.5, 0.05)
omegas = 2*np.pi/tau 
D_parallels = np.zeros((len(kappas), len(tau)))


D_k_o = np.load("../finite_element/data/vary_kappa_omega.npy")
contin_kappa = np.linspace(min(kappas), max(kappas), int(1e4))
D0 = np.zeros(len(omegas))
colors = plt.cm.viridis(np.linspace(0,1-0.00001,int(len(omegas))))#jet is also good

from numpy import *
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
for i in range(len(omegas)):
	omega = omegas[i]
	gamma   = np.sqrt(1j*omega/nu)
	rho     = np.sqrt(1j*omega/D)

	D0[i]    = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))
	#plt.plot(kappas, D_k_o[:, i], label=r"$\omega=%3.2f$" % omegas[i])
	res = np.sqrt(omegas[i]/(2*D))
	contin_D = interp1d(kappas, D_k_o[:, i], kind="cubic")(contin_kappa)
	if i == 0:
		plt.plot(res, contin_D[np.argmin(abs(contin_kappa-res))], "ko", markersize=3)#, label="Analytic resonance wavelength")
	else:
		plt.plot(res, contin_D[np.argmin(abs(contin_kappa-res))], "ko", markersize=3)
	plt.plot(contin_kappa, contin_D, label=r"$%3.2f$" % omegas[i], color=colors[i])

plt.xlabel(r" Wave number number $\kappa$", fontsize=8)
plt.ylabel(r" Second order Effective Diffusion coefficient $ D_\parallel^{(2)}$",  fontsize=8)
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc="best", fontsize=8, ncol=1, markerscale=0.1, title=r"$\rho^2$", columnspacing=1.0, labelspacing=0.8)
#plt.axis([0.45, 1.8, 0.42, 1.3])
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/analytic_resonance.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
