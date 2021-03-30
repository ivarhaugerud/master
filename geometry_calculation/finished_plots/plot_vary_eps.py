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


Lx = 4.48
kappa = 2*np.pi/Lx
epsilon = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
base = "../data_test/vary_eps/"
D          = np.zeros(len(epsilon))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = int(3.0/0.004)
tau = 3.0
omega = 2*np.pi/tau

for i in range(len(epsilon)-2):
	res = int(100*(1+2*epsilon[i]))
	if epsilon[i] < 0.65:
		data = np.loadtxt(base+"Lx4.48_tau3.0_eps"+str(epsilon[i])+"_nu1.2_D0.5_fzero0.0_fone12.0_res"+str(res)+"_dt0.004/tdata.dat")

	else:
		data = np.loadtxt(base+"Lx4.48_tau3.0_eps"+str(epsilon[i])+"_nu1.2_D0.5_fzero0.0_fone12.0_res130_dt0.003/tdata.dat")
	try:
		D[i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
		difference[i] = abs(D[i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[i]
		plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
	except:
		print("empty file", epsilon[i])
plt.show()

plt.plot(epsilon, difference)
plt.show()

plt.plot(epsilon, D, "o", markersize=3)
plt.xlabel(r" Womersley number $\frac{\omega a^2}{\nu}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/D_eff_vs_nu.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

"""
plt.figure(3)
for i in range(len(Lx)):
	plt.plot(nu, U[i, :], "o")

plt.xscale("log")
plt.ylabel(r"$\langle u^2 \rangle $")
plt.xlabel(r"viscosity $\nu$")
plt.legend(loc="best")


plt.figure(4)
for i in range(len(Lx)):
	plt.plot(nu, abs((D[i,:])/(U[i,:])))
plt.ylabel(r"$D/U^2$")
plt.xlabel(r"viscosity $\nu$")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="best")

plt.show()
"""