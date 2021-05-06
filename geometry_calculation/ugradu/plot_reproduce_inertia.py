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

epsilon = np.arange(0.1, 0.501, 0.1)
Lx = np.array([2.86, 3.49, 4.488, 6.28, 10.47])#, 31.41])   #
kappa = 2*np.pi/Lx

base = "data/vary_geometry/"
D          = np.zeros((len(epsilon), len(kappa)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = 750
tau = 3.0
omega = 2*np.pi/tau

for i in range(len(epsilon)):
	for j in range(len(kappa)):
		data = np.loadtxt(base+"Lx"+str(Lx[j]) + "_tau3.0_eps%1.1f" % (epsilon[i]) +"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt0.004/tdata.dat")
		D[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[i, j] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
		print(data[-1, 0])
		difference[i, j] = abs(D[j,i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[j, i]
		#plt.plot(data[:, 0]/tau,   data[:, 8])
		#plt.plot(data[-T:, 0]/tau, data[-T:, 8])
#plt.show()
plt.figure(1)
for i in range(len(epsilon)):
	plt.plot(kappa, difference[i, :], "o", color="C"+str(i), label=r"$\epsilon = %3.2f$" % epsilon[i], markersize=3)
	plt.plot(kappa, difference[i, :], color="C"+str(i), linewidth=1)
plt.yscale("log")
plt.legend(loc="best")
plt.show()

plt.figure(1)
for i in range(len(epsilon)):
	plt.plot(kappa, D[i, :], "o", color="C"+str(i), label=r"$\epsilon = %3.2f$" % epsilon[i], markersize=3)
	plt.plot(kappa, D[i, :], color="C"+str(i), linewidth=1)

plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", ncol=3, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "D_eff_vs_kappa.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


"""
for i in range(len(epsilon)):
	plt.plot(kappa, U[i, :]/1.2, "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
	plt.plot(kappa, U[i, :]/1.2, color="C"+str(i), linewidth=1)
	#plt.plot(interpool_nu, interpool)

plt.xlabel(r" Peclet number $\frac{aU}{D}$", fontsize=8)
plt.ylabel(r" Reynolds number Re",  fontsize=8)
plt.legend(loc="best", ncol=3, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "Re_vs_geo.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()



"""
plt.figure(1)
for i in range(len(epsilon)):
	plt.plot(kappa, D[i, :], "-", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)

D_re = np.copy(D)
kappa_original = np.copy(kappa)
kappa_cont = np.linspace(min(kappa_original), max(kappa_original), int(1e4))

epsilon = np.arange(0.1, 0.501, 0.10)
Lx = np.array([12.0, 9.0, 6.0, 3.0])
kappa = 2*np.pi/Lx 
tau = 3.0

D = np.zeros((len(kappa), len(epsilon)))
difference = np.zeros(np.shape(D))
T = 750 
base = "../data_test/vary_geometry/"
plt.figure(1)

for i in range(len(kappa)):
	for j in range(len(epsilon)):
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau3.0_eps"+str(epsilon[j])[:3]+"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt0.004/tdata.dat")
		D[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau
		U[i, j] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau

for i in range(len(epsilon)):
	plt.plot(kappa, D[:, i], "o", color="C"+str(i), markersize=3)


plt.xlabel(r" Peclet number $\frac{aU}{D}$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", ncol=3, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "compare_ugradu.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()




for i in range(len(epsilon)):
	print(kappa_original, kappa_cont)
	interpoo = scint.interp1d(kappa_original, D_re[i, :], kind="cubic")(kappa_cont)
	print(kappa_original, kappa)
	for j in range(len(kappa)-1):
		j += 1
		if j == 1:
			plt.plot(U[j,i]/1.2, 100*abs(D[j, i]-interpoo[np.argmin(abs(kappa_cont-kappa[j]))])/D[j, i], "o", color="C"+str(i), label=r"$\epsilon=%1.1f$" % epsilon[i], markersize=3)
		else:
			plt.plot(U[j,i]/1.2, 100*abs(D[j, i]-interpoo[np.argmin(abs(kappa_cont-kappa[j]))])/D[j, i], "o", color="C"+str(i), markersize=3)
plt.xlabel(r" Reynolds number $\frac{aU}{\nu}$", fontsize=8)
plt.ylabel(r" Rel change in $D_\parallel$ with inertia [%]",  fontsize=8)
plt.legend(loc="best", ncol=2, fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "ugradu_change.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
