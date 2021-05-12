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

tau = np.array([2, 2.6, 3.69])
epsilon = "0.5"
dt = tau/500

Lx = np.array([2.855, 3.49, 4.487, 4.833, 5.235, 5.711, 6.283, 6.613, 8.377, 10.47, 13.96]) #, 7.853
tau = np.array([2.25, 2.6, 3.0, 3.6, 4.5])
kappa = 2*np.pi/Lx 
epsilon = "0.3"
base = "../data_test/find_resonance_try5/"
D = np.zeros((len(kappa), len(tau)))
U = np.zeros(np.shape(D))
D0 = np.zeros(np.shape(D))
difference = np.zeros(np.shape(D))
dt = tau/750
Ts = np.ones(len(tau), dtype="int")*int(500)
T = 750 

nu = 1.2
F0 = 12/nu
Pe = 1
Dm  = 1
for i in range(len(tau)):
	omega = 2*np.pi/tau[i]
	rho = np.sqrt(1j*omega/Dm)
	rho_c = np.conjugate(rho)
	gamma = np.sqrt(1j*omega/nu)
	gamma_c = np.conjugate(gamma)
	D0[:,i] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))


plt.figure(1)
kappa_cont = np.linspace(min(kappa), max(kappa), int(1e4))

for i in range(len(kappa)):
	for j in range(len(tau)):
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau"+ str(round(tau[j], 3)) +"_eps"+epsilon+"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt" + str(round(dt[j], 6)) + "/tdata.dat")
		print(kappa[i], tau[j], np.shape(data), np.shape(D))
		D[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau[j]
		difference[i, j] = abs(D[i, j] - sci.trapz(  np.trim_zeros(data[:, 8])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau[j])/D[i,j]
		U[i, j] = sci.trapz(  data[:, 4][-T:],  data[:, 0][-T:] )/tau[j]
		#if D[i, j] > 1.0:	
		#plt.plot(np.trim_zeros(data[:, 0])/tau[j], np.trim_zeros(data[:, 8]))
		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[j], np.trim_zeros(data[:, 8])[-T:])
		plt.title(str(kappa[i]) + ","+ str(tau[j]))
		plt.xlabel(r" Time [periods]", fontsize=8)
		plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.show()

plt.figure(2)
for i in range(len(tau)):
	plt.plot(kappa, U[:, i]/1.2, "o", markersize=3, label=r"$\tau=%3.2f$" % tau[i])
#plt.yscale("log")
plt.show()

G = np.sqrt(2*np.pi/(tau*1.2))
R   = np.sqrt(2*np.pi/(tau*1))


tau2 = np.array([7.0])
Lx2 = np.array([8.377, 8.975, 9.666, 10.47])
kappa2 = 2*np.pi/Lx2
dt = tau/500

omega = 2*np.pi/tau2[0]
rho = np.sqrt(1j*omega/Dm)
rho_c = np.conjugate(rho)
gamma = np.sqrt(1j*omega/nu)
gamma_c = np.conjugate(gamma)
D2_0 = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))

D2 = np.zeros((len(kappa2), len(tau2)))
difference2 = np.zeros(np.shape(D2))
dt2 = tau2/750
T = 750 
rho2 = np.sqrt(2*np.pi/tau2[0])

for i in range(len(Lx2)):
	for j in range(len(tau2)):
		data = np.loadtxt(base+"Lx" +  str(Lx2[i]) +"_tau"+ str(round(tau2[j], 3)) +"_eps"+epsilon+"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt" + str(round(dt2[j], 6)) + "/tdata.dat")
		D2[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau2[j]
		difference2[i, j] = abs(D2[i, j] - sci.trapz(  np.trim_zeros(data[:, 8])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau2[j])/D2[i,j]



plt.figure(3)
for i in range(len(tau)):
	plt.plot(kappa, D[:, i], "o", color="C"+str(i), markersize=3, label=r"Wo$_D=%3.2f$" % R[i])
	plt.plot(kappa, D[:, i], color="C"+str(i))
plt.plot(kappa2, D2[:, 0], color="C"+str(len(tau)), markersize=3)
plt.plot(kappa2, D2[:, 0], "o", markersize=3, color="C"+str(len(tau)), label=r"Wo$_D=%3.2f$" % rho2)

#plt.axis([0.65, 1.45, 0.84, 0.98])
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/semi_ana_resonance.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

"""
plt.figure(5)
for i in range(len(tau)):
	plt.plot(kappa, D0[:,i], "o", color="C"+str(i), markersize=3, label=r"$\rho=%3.2f$" % R[i])
	plt.plot(kappa, D0[:,i], color="C"+str(i))
plt.plot(kappa2, D2_0*np.ones(len(kappa2)), color="C"+str(len(tau)), markersize=3)
plt.plot(kappa2, D2_0*np.ones(len(kappa2)), "o", markersize=3, color="C"+str(len(tau)), label=r"$\rho=%3.2f$" % rho2)

#plt.axis([0.65, 1.45, 0.84, 0.98])
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Change in dispersion Coefficient $ D_\parallel-D_\parallel^{(0)} $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/semi_ana_resonance_change.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
"""
plt.figure(5)
for i in range(len(tau)):
	plt.plot(kappa, D[:, i]-D0[:,i], "o", color="C"+str(i), markersize=3, label=r"Wo$_D=%3.2f$" % R[i])
	plt.plot(kappa, D[:, i]-D0[:,i], color="C"+str(i))
plt.plot(kappa2, D2[:, 0]-D2_0, color="C"+str(len(tau)), markersize=3)
plt.plot(kappa2, D2[:, 0]-D2_0, "o", markersize=3, color="C"+str(len(tau)), label=r"Wo$_D=%3.2f$" % rho2)

#plt.axis([0.65, 1.45, 0.84, 0.98])
plt.xlabel(r"Wave number $\kappa$", fontsize=8)
plt.ylabel(r"Change in Dispersion $ D_\parallel-D_\parallel^{(0)} $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/semi_ana_resonance_change.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


plt.figure(4)
kappa_cont = np.linspace(min(kappa), max(kappa), int(1e4))
for i in range(len(tau)):
	grad = np.gradient(D[:, i], kappa)
	interpoo = scp.interp1d(kappa, grad, "linear")(kappa_cont)
	plt.plot(kappa, np.gradient(D[:, i], kappa), "o", color="C"+str(i), markersize=3, label=r"$\rho=%3.2f$" % R[i])
	plt.plot(kappa_cont, interpoo,      color="C"+str(i), markersize=3)

grad = np.gradient(D2[:, 0], kappa2)
kappa_cont2 = np.linspace(min(kappa2), max(kappa2), int(1e4))
interpoo = scp.interp1d(kappa2, grad, "linear")(kappa_cont2)
plt.plot(kappa2, np.gradient(D2[:, 0], kappa2), "o", color="C"+str(i), markersize=3, label=r"$\rho=%3.2f$" % R[i])
plt.plot(kappa_cont2, interpoo,      color="C"+str(i), markersize=3)

plt.plot(kappa, np.zeros(len(kappa)), "k")
plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Dispersion $ D_\parallel $ [$D_m$]",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/interpolated_resonance_deriv.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()




kappa_res = np.zeros(len(tau)+1)

for i in range(len(tau)):
	grad = np.gradient(D[:, i], kappa)
	interpoo = scp.interp1d(kappa, grad, "cubic")(kappa_cont)
	kappa_res[i] = kappa_cont[np.argmin(abs(interpoo))]

grad = np.gradient(D2[:, 0], kappa2)
interpoo = scp.interp1d(kappa2, grad, "cubic")(kappa_cont2)
kappa_res[-1] = kappa_cont2[np.argmin(abs(interpoo))]
tau = np.array([2.25, 2.6, 3.0, 3.6, 4.5,  7.0])
#tau_long = np.linspace(0.5, 15, int(1e5))
plt.plot( np.sqrt(2*np.pi/tau), kappa_res, "o", markersize=3, label="Numeric")
#plt.plot( tau_long + 0*np.sqrt(2*np.pi/tau_long),  np.sqrt(2*np.pi/(2*tau_long)),    label="Analytic")
plt.plot( np.sqrt(2*np.pi/tau),  np.sqrt(2*np.pi/(2*tau)),    label=r"$\sqrt{\omega a^2/(2D_m)}$")
plt.xlabel(r" Diffusive Womersley number $\sqrt{\omega a^2/D}$", fontsize=8)
plt.ylabel(r" Resonance wave number $ \kappa_{res}$",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/test_analytic_resonance.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()