import os
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import scipy.interpolate as scint
import scipy.interpolate as scp 

def linear_regresion(x, y):
    n = float(len(x))
    D = float(np.sum(np.square(x)) - (np.sum(x)**2)/n)
    E = float(np.sum(x*y) - np.sum(x)*np.sum(y)/n)
    F = float(np.sum(np.square(y)) - (np.sum(y)**2)/n)

    delta_m = np.sqrt((1/(n-2))*(D*F-E**2)/(D**2))
    delta_c = np.sqrt(1/(n-2)*(D/n+np.mean(x)**2)*(D*F-E**2)/(D**2))
    m = E/D
    c = np.mean(y)-m*np.mean(x)

    return m, c, delta_m, delta_c
    #using linear regression from Squires, with uncertainty to find slope and constant term


plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	
root = "../../../master_latex/results/"

"""
Lx = np.array([3.92, 4.48, 5.23, 5.71, 6.28, 6.98, 7.85, 10.47])
F0 = np.array([12.0, 18.0, 24.0, 30.0, 36.0])
kappa = 2*np.pi/Lx 

base = "../data_test/vary_F/"
D          = np.zeros((len(kappa), len(F0)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))
res        = np.zeros(len(F0))
res_cont   = np.zeros(len(F0))

T = 750
tau = 3.0
omega = 2*np.pi/tau

for j in range(len(kappa)):
	for i in range(len(F0)):
		try:
			data = np.loadtxt(base+"Lx"+str(Lx[j]) + "_tau3.0_eps0.3_nu1.2_D1.0_fzero0.0_fone" + str(F0[i]) + "_res150_dt0.004/tdata.dat")
			T = 750
			D[j, i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
			U[j, i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
			#plt.scatter(nu[i], data[-1, 0]/tau)
			difference[j, i] = abs(D[j,i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[j, i]
			#plt.plot(np.trim_zeros(data[:, 0])/tau, np.trim_zeros(data[:, 8]))
			#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
			#plt.show()
		except:
			print("Did not work for ", F0[i], kappa[j])

plt.figure(1)
kappa_cont = np.linspace(min(kappa), max(kappa), int(1e4))

for i in range(len(F0)):
	plt.plot(kappa, D[:, i], "o", color="C"+str(i), label=r"$F_0 = %3.2f$" % F0[i], markersize=3)
	#plt.plot(kappa, D[:, i], color="C"+str(i), linewidth=1)
	#res[i] = kappa[np.argmax(D[:, i])]
	res[i] = (kappa[np.argmax(D[:,i])])
	interpoo = scp.interp1d(kappa, D[:, i], "cubic")(kappa_cont)
	res_cont[i] = kappa_cont[np.argmax(interpoo)]
	plt.plot(kappa_cont, interpoo, color="C"+str(i), linewidth=1)

	#plt.plot(kappa, np.gradient(D[:, i], kappa), "o", color="C"+str(i), markersize=3, label=r"$\rho=%3.2f$" % R[i])
	#plt.plot(kappa_cont, interpoo,      color="C"+str(i), markersize=3)
	#plt.plot(interpool_nu, interpool)

plt.xlabel(r" Womersley number $\sqrt{\frac{\omega a^2}{\nu}}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_eff_vs_F0.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()



plt.plot(F0, res, "o")
plt.plot(F0, res_cont)
plt.show()


plt.figure(1)
kappa_ind = [0, 1, 2, 4, 6, 7]
for i in range(len(kappa_ind)):
	ind = kappa_ind[i]
	plt.plot(U[ind,:], D[ind, :], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[ind], markersize=3)
	plt.plot(U[ind,:], D[ind, :], color="C"+str(i), linewidth=1)
	m, c, delta_m, delta_c = linear_regresion(U[ind,:], D[ind, :])
	print(m, delta_m)
	#plt.plot(interpool_nu, interpool)

#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r" Womersley number $\sqrt{\frac{\omega a^2}{\nu}}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_eff_vs_Pe_varyF0.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


plt.figure(2)
for i in range(len(F0)):
	plt.plot(kappa, difference[:, i], "o", color="C"+str(i), label=r"$F_0 = %3.2f$" % F0[i], markersize=3)
plt.yscale("log")
plt.legend(loc="best")
plt.show()


plt.figure(3)
for i in range(len(kappa)):
	plt.plot(F0, U[i, :], label=r"$\kappa = %3.2f$" % kappa[i])

plt.ylabel(r"$\langle u^2 \rangle $")
plt.xlabel(r"viscosity $\nu$")
plt.legend(loc="best")
plt.show()
"""



Lx = np.array([3.14, 3.92, 5.23, 7.85]) #4.48, 5.23, 5.71, 6.28, 6.98, 7.85, 10.47])
F0 = np.array([3.981, 7.437, 13.89, 25.95, 48.49, 90.6, 169.2, 316.2])
kappa = 2*np.pi/Lx 
Dm    = 3 

base = "../data_test/vary_F0_2/"
D          = np.zeros((len(kappa), len(F0)))
amp        = np.zeros((len(kappa), len(F0)))
difference = np.zeros(np.shape(D))
U          = np.zeros(np.shape(D))
res        = np.zeros(len(F0))
res_cont   = np.zeros(len(F0))

T = 750
tau = 3.0
omega = 2*np.pi/tau

for j in range(len(kappa)):
	for i in range(len(F0)):
		#try:
		data = np.loadtxt(base+"Lx"+str(Lx[j]) + "_tau3.0_eps0.3_nu2.0_D3.0_fzero0.0_fone" + str(F0[i]) + "_res150_dt0.004/tdata.dat")
		T = 750
		D[j, i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau
		U[j, i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau
		amp[j, i] = (max(data[-T:, 8]) - np.mean(data[-T:, 8]))
		#plt.plot(data[:, 0], data[:, 8])
		#plt.plot(data[-T:, 0], data[-T:, 8])
		#plt.plot(data[-T:, 0], np.ones(T)*D[j,i], "k")
		#plt.plot(data[-1, 0], D[j, i]+amp[j,i], "ko")
		#plt.show()
		difference[j, i] = abs(D[j,i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau)/D[j, i]
		#difference[j, i] = abs(D[j,i] - np.mean(  data[-T:, 8]))/D[j, i]
		#except:
		#	print("Did not work for ", F0[i], kappa[j])


plt.figure(1)
for i in range(len(kappa)):
	plt.plot(U[i,:]/Dm, amp[i, :]/D[i,:], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
	plt.plot(U[i,:]/Dm, amp[i, :]/D[i,:], color="C"+str(i), linewidth=1)

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r" Peclet number $\frac{aU}{D}$", fontsize=8)
plt.ylabel(r" Relative dispersion amplitude",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_eff_vs_Pe_varyF0.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()



plt.figure(1)
for i in range(len(kappa)):
	plt.plot(U[i,:]/Dm, D[i,:], "o", color="C"+str(i), label=r"$\kappa = %3.2f$" % kappa[i], markersize=3)
	plt.plot(U[i,:]/Dm, D[i,:], color="C"+str(i), linewidth=1)

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r" Peclet number $\frac{aU}{D}$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/D_eff_vs_Pe_varyF0.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

plt.figure(2)
for i in range(len(F0)):
	plt.plot(kappa, difference[:, i], "o", color="C"+str(i), label=r"$F_0 = %3.2f$" % F0[i], markersize=3)
plt.yscale("log")
plt.legend(loc="best")
plt.show()
"""
plt.figure(3)
for i in range(len(kappa)):
	plt.plot(F0, U[i, :], label=r"$\kappa = %3.2f$" % kappa[i])

plt.ylabel(r"$\langle u^2 \rangle $")
plt.xlabel(r"viscosity $\nu$")
plt.legend(loc="best")
plt.show()
"""