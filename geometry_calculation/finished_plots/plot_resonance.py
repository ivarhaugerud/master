import os
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import scipy.interpolate as scp 

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

#Lx = np.array([4.487, 4.654, 4.833, 5.026, 5.235, 5.463, 5.711, 5.983, 6.283, 6.613, 6.981, 7.391, 7.853, 8.377, 8.975])
Lx = np.array([2.855, 3.49, 4.487, 5.235, 6.283, 10.47]) #, 7.853
tau = np.array([2.25, 2.6, 3.0, 3.6, 4.5])
kappa = 2*np.pi/Lx 
epsilon = "0.3"
base = "../data_test/find_resonance_try5/"
D = np.zeros((len(kappa), len(tau)))
difference = np.zeros(np.shape(D))
dt = tau/750
Ts = np.ones(len(tau), dtype="int")*int(500)
T = 750 

plt.figure(1)
for i in range(len(kappa)):
	for j in range(len(tau)):
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau"+ str(round(tau[j], 3)) +"_eps"+epsilon+"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt" + str(round(dt[j], 6)) + "/tdata.dat")
		print(kappa[i], tau[j], np.shape(data), np.shape(D))
		D[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau[j]
		difference[i, j] = abs(D[i, j] - sci.trapz(  np.trim_zeros(data[:, 8])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau[j])/D[i,j]

		#if D[i, j] > 1.0:
		#plt.plot(np.trim_zeros(data[:, 0])/tau[j], np.trim_zeros(data[:, 8]))
		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[j], np.trim_zeros(data[:, 8])[-T:])
		plt.title(str(kappa[i]) + ","+ str(tau[j]))
		plt.xlabel(r" Time [periods]", fontsize=8)
		plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
#plt.show()

gamma = np.sqrt(2*np.pi/(tau*1.2))
rho   = np.sqrt(2*np.pi/(tau*1))

plt.figure(2)
for i in range(len(tau)):
	plt.plot(kappa, difference[:, i], "o", markersize=3, label=r"$\rho=%3.2f$" % rho[i])
plt.yscale("log")
#plt.show()

plt.figure(3)
for i in range(len(tau)):
	plt.plot(kappa, D[:, i], "o", color="C"+str(i), markersize=3, label=r"$\rho=%3.2f$" % rho[i])
	plt.plot(kappa, D[:, i], color="C"+str(i))

#plt.axis([0.65, 1.45, 0.84, 0.98])
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/semi_ana_resonance.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


plt.figure(4)
kappa_cont = np.linspace(min(kappa), max(kappa), int(1e4))
for i in range(len(tau)):
	grad = np.gradient(D[:, i], kappa)
	interpoo = scp.interp1d(kappa, grad, "cubic")(kappa_cont)
	plt.plot(kappa, np.gradient(D[:, i], kappa), "o", color="C"+str(i), markersize=3, label=r"$\rho=%3.2f$" % rho[i])
	plt.plot(kappa_cont, interpoo,      color="C"+str(i), markersize=3)


plt.plot(kappa, np.zeros(len(kappa)), "k")
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/interpolated_resonance_deriv.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()




kappa_res = np.zeros(len(tau))

for i in range(len(tau)):
	grad = np.gradient(D[:, i], kappa)
	interpoo = scp.interp1d(kappa, grad, "cubic")(kappa_cont)
	kappa_res[i] = kappa_cont[np.argmin(abs(interpoo))]

#plt.plot(np.sqrt(2*np.pi/tau[i]), kappa_cont[np.argmin(abs(interpoo))], "o", color="C"+str(i), markersize=3, label=r"$\rho=%3.2f$" % rho[i])

plt.plot( np.sqrt(2*np.pi/tau), kappa_res, "o", markersize=3, label="Numeric")
plt.plot( np.sqrt(2*np.pi/tau),  np.sqrt(2*np.pi/(2*tau)), label="Analytic" )
plt.xlabel(r" Womersley number $\$", fontsize=8)
plt.ylabel(r" Resonance wave number $ \kappa_{res} $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/test_analytic_resonance.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
