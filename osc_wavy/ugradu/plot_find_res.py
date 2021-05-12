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
root = "../../../master_latex/results/figures/ugradu/"

Lx = np.flip(np.array([3.49, 4.488, 6.28, 8.5, 10.47, 15.0]))# 15.7, 31.41 
Dm  = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0])
kappa = 2*np.pi/Lx

base = "data/vary_geometry_res/"
D          = np.zeros((len(Dm), len(kappa)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))
amp        = np.zeros(np.shape(D))

T = 750
tau = 3.0
omega = 2*np.pi/tau

for i in range(len(Dm)):
	for j in range(len(kappa)):
		try:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + '_tau3.0_eps0.3_nu2.25_D'+str(Dm[i]) + '_fzero0.0_fone1.0_res150_dt0.004/tdata.dat')
		except:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + '_tau3.0_eps0.3_nu2.25_D'+str(Dm[i]) + '_fzero0.0_fone1.0_res150_dt0.004/tdata.dat')
		D[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[i, j] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
		amp[i, j] = np.max( abs(-D[i,j] +  data[-T:, 8]))/D[i,j]
		difference[i, j] = abs(D[i,j] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[i, j]
		print(data[-1, 0])
		plt.plot(data[:, 0]/tau,   data[:, 8])
		plt.plot(data[-T:, 0]/tau, data[-T:, 8])
		print(amp[i,j])
plt.show()
#plt.show()
plt.figure(1)
for i in range(len(Dm)):
	plt.plot(kappa, difference[i, :])
plt.xscale("log")
plt.yscale("log")
plt.show()


plt.figure(1)
for i in range(len(Dm)):
	plt.plot(kappa, D[i, :], "o", color="C"+str(i), label=r"$D_m = %3.2f$" % Dm[i], markersize=3)
	plt.plot(kappa, D[i, :], color="C"+str(i), linewidth=1)
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", ncol=3, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "D_eff_vs_Re.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))


plt.figure(2)
indx = 0
for i in range(len(kappa)):
	if abs(kappa[i]-0.74) > 0.02: 
		plt.plot(np.sqrt(omega/Dm), D[:, i], "o", color="C"+str(indx), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
		plt.plot(np.sqrt(omega/Dm), D[:, i], color="C"+str(indx), linewidth=1)
		print(U[:, i]/Dm)
		indx += 1
plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r" Diffusive Womersley number $\sqrt{\omega a^2/D_m}$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", ncol=2, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "D_eff_vs_rho.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))


plt.figure(17)
indx = 0
for i in range(len(kappa)):
	if abs(kappa[i]-0.74) > 0.02: 
		plt.plot(np.sqrt(omega/Dm), amp[:, i], "o", color="C"+str(indx), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
		plt.plot(np.sqrt(omega/Dm), amp[:, i], color="C"+str(indx), linewidth=1)
		print(U[:, i]/Dm)
		indx += 1
plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r" Diffusive Womersley number $\sqrt{\omega a^2/D_m}$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", ncol=2, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "D_amp_vs_rho.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))



plt.figure(3)
indx = 0
for i in range(len(kappa)):
	if abs(kappa[i]-0.74) > 0.02: 
		plt.plot(U[:, i]/Dm, D[:, i], "o", color="C"+str(indx), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
		plt.plot(U[:, i]/Dm, D[:, i], color="C"+str(indx), linewidth=1)
		indx += 1

plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r" Peclet number $aU/D_m$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", ncol=2, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "D_eff_vs_Pe_and_rho.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()