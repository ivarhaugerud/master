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
nu = np.logspace(-2, 2, 10)[:9]

base = "../data_test/vary_nu_try2/"
D          = np.zeros((len(kappa), len(nu)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))

T = 750
tau = 3.0
omega = 2*np.pi/tau
print(nu)

for j in range(len(kappa)):
	for i in range(len(nu)):
		try:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + "_tau3.0_eps0.5_nu" + str(nu[i])[:4]+"_D1.0_fzero0.0_fone" + str(12*nu[i])[:5] + "_res150_dt0.004/tdata.dat")
			print(np.shape(data))
			T = 750
			D[j, i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
			U[j, i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
			#plt.scatter(nu[i], data[-1, 0]/tau)
			difference[j, i] = abs(D[j,i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[j, i]
			#plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
			#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
			#plt.show()
		except:
			print("Did not work for ", nu[i], kappa[j])

interpool_nu = np.logspace(np.log10(min(nu)), np.log10(max(nu)), int(1e4))
plt.figure(1)
for i in range(len(Lx)):
	ind = len(Lx)-1-i
	plt.plot(np.sqrt(omega/nu), D[ind, :], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[ind], markersize=3)
	plt.plot(np.sqrt(omega/nu), D[ind, :], color="C"+str(i), linewidth=1)
	#plt.plot(interpool_nu, interpool)

plt.xlabel(r" Womersley number $\sqrt{\frac{\omega a^2}{\nu}}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xscale("log")
filename = root + "figures/D_eff_vs_nu.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
#plt.show()

plt.figure(2)
for i in range(len(Lx)):
	ind = len(Lx)-1-i
	plt.plot(np.sqrt(omega/nu), difference[ind, :], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[ind], markersize=3)
	#plt.plot((omega/nu), D[ind, :], color="C"+str(i), linewidth=1)
plt.yscale("log")
#plt.show()


plt.figure(3)
for i in range(len(Lx)):
	plt.plot(omega/nu, U[i, :], label=r"$\kappa = %3.2f$" % kappa[i])

plt.xscale("log")
plt.ylabel(r"$\langle u^2 \rangle $")
plt.xlabel(r"viscosity $\nu$")
plt.legend(loc="best")


plt.figure(4)
g = (D - 1)*(105/2)
for i in range(len(Lx)):
	plt.plot(nu, abs(g[i, :]*nu/U[i,:]))
plt.ylabel(r"$D/U^2$")
plt.xlabel(r"viscosity $\nu$")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="best")




plt.figure(5)
for i in range(len(Lx)):
	ind = len(Lx)-1-i
	plt.plot(U[ind,:], D[ind, :], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[ind], markersize=3)
	plt.plot(U[ind,:], D[ind, :], color="C"+str(i), linewidth=1)
	#plt.plot(interpool_nu, interpool)

plt.xlabel(r" Womersley number $\frac{\omega a^2}{\nu}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xscale("log")
#filename = root + "figures/D_eff_vs_nu.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))

plt.show()