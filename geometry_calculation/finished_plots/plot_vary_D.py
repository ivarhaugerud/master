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


Lx = np.array([1.96, 2.41, 3.14, 4.48, 7.85])
kappa = 2*np.pi/Lx
kappa[0] = 3.20
kappa[1] = 2.60
Dm       = np.logspace(-2, 2, 10)

base = "../data_test/vary_D/"
D          = np.zeros((len(kappa), len(Dm)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = int(3/0.004)
tau = 3.0
omega = 2*np.pi/tau

for j in range(len(kappa)):
	for i in range(len(Dm)):
		try:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + "_tau3.0_eps0.5_nu1.2_D"+str(Dm[i])[:4]+"_fzero0.0_fone12.0_res125_dt0.004/tdata.dat")
			D[j, i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
			U[j, i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
			difference[j, i] = abs(D[j,i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[j, i]
			#plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
			#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])

		except:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + "_tau3.0_eps0.5_nu1.2_D"+str(Dm[i])[:5]+"_fzero0.0_fone12.0_res125_dt0.004/tdata.dat")
			D[j, i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
			U[j, i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
			difference[j, i] = abs(D[j,i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[j, i]
			#plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
			#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
#plt.show()

plt.figure(1)
for i in range(len(Lx)):
	ind = len(Lx)-1-i
	plt.plot(U[ind,:]/Dm, D[ind, :], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[ind], markersize=3)
	plt.plot(U[ind,:]/Dm, D[ind, :], color="C"+str(i), linewidth=1)
	#plt.plot(interpool_nu, interpool)

plt.xlabel(r" Peclet number $\frac{aU}{D}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xscale("log")
#plt.yscale("log")
filename = root + "figures/D_eff_vs_Pe.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


plt.figure(2)
for i in range(len(Lx)):
	ind = len(Lx)-1-i
	plt.plot(np.sqrt(omega/Dm), D[ind, :], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[ind], markersize=3)
	plt.plot(np.sqrt(omega/Dm), D[ind, :], color="C"+str(i), linewidth=1)

plt.xlabel(r" Diffusive Womersley number $\sqrt{\frac{\omega a^2}{D}}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel$",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xscale("log")
filename = root + "figures/D_eff_vs_rho.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))

plt.show()
"""
plt.figure(2)
for i in range(len(Lx)):
	ind = len(Lx)-1-i
	plt.plot(Dm, difference[ind, :], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[ind], markersize=3)
	#plt.plot((omega/nu), D[ind, :], color="C"+str(i), linewidth=1)
plt.legend(loc="best")
plt.yscale("log")
plt.xscale("log")
plt.show()
"""