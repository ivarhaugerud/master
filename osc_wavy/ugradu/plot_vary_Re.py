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

Lx = np.array([3.49, 4.488, 6.28, 7.85, 10.47, 15.7])# 
F  = np.logspace(0, 3, 7)
kappa = 2*np.pi/Lx

base = "data/vary_geometry_vary_F/"
D          = np.zeros((len(F), len(kappa)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = 750
tau = 3.0
omega = 2*np.pi/tau

for i in range(len(F)):
	for j in range(len(kappa)):
		print("%4.1f" % F[i])
		try:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + '_tau3.0_eps0.4_nu2.25_D1.0_fzero0.0_fone%0.1f'% F[i] + "_res150_dt0.004/tdata.dat")
		except:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + '_tau3.0_eps0.4_nu2.25_D1.0_fzero0.0_fone%0.2f'% F[i] + "_res150_dt0.004/tdata.dat")
		D[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[i, j] = np.sqrt(sci.trapz(  data[-T:, 4],  data[-T:, 0] ))/tau
		difference[i, j] = abs(D[i,j] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[i, j]
		#plt.plot(data[:, 0]/tau,   data[:, 8])
		#plt.plot(data[-T:, 0]/tau, data[-T:, 8])
	#plt.show()
plt.figure(1)
for i in range(len(kappa)):
	plt.plot(U[:,i]/1.0, difference[:, i], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
	plt.plot(U[:,i]/1.0, difference[:, i], color="C"+str(i), linewidth=1)
plt.xscale("log")
plt.yscale("log")
plt.show()

plt.figure(1)
for i in range(len(kappa)):
	plt.plot(U[:, i]/2.25, D[:, i]/(1+2/105*(U[:,i]/1)**2), "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
	plt.plot(U[:, i]/2.25, D[:, i]/(1+2/105*(U[:,i]/1)**2), color="C"+str(i), linewidth=1)
plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r" Reynolds number $aU/\nu$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", ncol=2, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "D_eff_vs_Re.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))


plt.figure(2)
for i in range(len(kappa)):
	plt.plot(U[:, i]/1, D[:, i]/(1+2*U[:, i]*U[:, i]/105), "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
	plt.plot(U[:, i]/1, D[:, i]/(1+2*U[:, i]*U[:, i]/105), color="C"+str(i), linewidth=1)

plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r" Peclet number $aU/D_m$", fontsize=8)
plt.ylabel(r" Relative effecitve dispersion $ D_\parallel/D_\parallel^{aris} $",  fontsize=8)
plt.legend(loc="best", ncol=2, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "D_eff_vs_Pe.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))





plt.figure(3)
for i in range((5)):
	plt.plot(kappa, D[i, :]/(1+2/105*(U[i,:]/1)**2), "o", color="C"+str(i), label=r"Pe= %3.2f" % np.mean(U[i, :]/1), markersize=3)
	plt.plot(kappa, D[i, :]/(1+2/105*(U[i,:]/1)**2), color="C"+str(i), linewidth=1)
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Relative effecitve dispersion $ D_\parallel/D_\parallel^{aris} $",  fontsize=8)
plt.axis([0.35, 1.85, 0.82, 1.75])
plt.legend(loc="best", ncol=2, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "D_eff_vs_kappa_and_Pe.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()