import os
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	
root = "../../../master_latex/results/"

#set 0
tau = np.array([2, 2.6, 3.69])
epsilon = "0.5"
dt = tau/500

#set 1
"""
tau = np.array([2.25, 3.0, 4.5])
epsilon = "0.3"
dt = tau/750
"""
#set 2 
"""
tau = np.array([1.0, 7.38, 12.56])
epsilon = "0.5"
dt = tau/750
"""

Lx = np.array([4.487, 4.654, 4.833, 5.026, 5.235, 5.463, 5.711, 5.983, 6.283, 6.613, 6.981, 7.391, 7.853, 8.377, 8.975])
kappa = 2*np.pi/Lx 

base = "../data_test/find_resonance_try4/"
D = np.zeros((len(kappa), len(tau)))
difference = np.zeros(np.shape(D))
#dt = tau/500
Ts = np.ones(len(tau), dtype="int")*int(500)

for i in range(len(kappa)):
	for j in range(len(tau)):
		T = Ts[j]
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau"+ str(round(tau[j], 3)) +"_eps"+epsilon+"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt" + str(round(dt[j], 6)) + "/tdata.dat")
		print(kappa[i], tau[j], np.shape(data), np.shape(D))
		D[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau[j]
		difference[i, j] = abs(D[i, j] - sci.trapz(  np.trim_zeros(data[:, 8])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau[j])/D[i,j]

		#if D[i, j] > 1.0:
		plt.plot(np.trim_zeros(data[:, 0])/tau[j], np.trim_zeros(data[:, 8]))
		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[j], np.trim_zeros(data[:, 8])[-T:])
		plt.title(str(kappa[i]) + ","+ str(tau[j]))
		plt.xlabel(r" Time [periods]", fontsize=8)
		plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
		plt.show()

gamma = np.sqrt(2*np.pi/(tau*1.2))
rho   = np.sqrt(2*np.pi/(tau*1))


for i in range(len(tau)):
	plt.plot(kappa, difference[:, i], "o", markersize=3, label=r"$\gamma=%3.2f$" % rho[i])
plt.yscale("log")
plt.show()

for i in range(len(tau)):
	plt.plot(kappa, D[:, i], "o", markersize=3, label=r"$\gamma=%3.2f$" % rho[i])

#plt.axis([0.65, 1.45, 0.84, 0.98])
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8, ncol=len(tau))
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/semi_ana_resonance.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()



for i in range(len(tau)):
	plt.plot(kappa, np.gradient(D[:, i], kappa), "o", markersize=3, label=r"$\gamma=%3.2f$" % rho[i])

#plt.axis([0.65, 1.45, 0.84, 0.98])
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8, ncol=len(tau))
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/semi_ana_resonance.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
