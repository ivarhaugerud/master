import os
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import scipy.interpolate as scp
import scipy.signal as scs

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	
root = "../../../master_latex/results/"



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



#run 1
D_parallels = np.load("D_para_vary_kappa_D_F0_tau_nu1000.npy")
D_parallels = np.load("D_para_vary_kappa_D_F0_tau_nu1000_rerun.npy")
D_parallels = np.load("D_para_vary_kappa_D_F0_tau_nu1000_rerun.npy")[:, :, 1:, 2:]

taus      = np.array([10, 20, 30, 40])
omegas    = 2*np.pi/taus
nu        = 1000
F0s       = np.linspace(500, 2500, 5)[1:]/nu
Ds        = np.logspace(-2.5, 0, 6)[2:]
Ds_logval = np.linspace(-2.5, 0, 6)[2:]
kappas    = np.arange(0.1, 1.501, 0.10)
print(np.shape(D_parallels))
#D_parallels = np.zeros((len(kappas), len(taus), len(F0s), len(Ds)))


kapppa_cont = np.linspace(min(kappas), max(kappas), int(1e4))
kappa_res = np.zeros((len(taus), len(Ds), len(F0s)))

U = np.zeros((len(taus), len(F0s)))
Re = np.zeros(np.shape(U))
xi = np.linspace(-1, 1, int(1e5))
eps = 0.3


for i in range(len(omegas)):
	for j in range(len(F0s)):
		omega = 2*np.pi/taus[i]
		gamma = np.sqrt(1j*omega/nu)
		F0 = F0s[j]

		ux0 = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
		U[i, j] = 0.5*sci.trapz(np.real(ux0*np.conjugate(ux0)), xi)
		U[i, j] = np.sqrt(U[i, j])

"""
for i in range(len(taus)):
	omega = 2*np.pi/taus[i]
	gamma   = np.sqrt(1j*omega/nu)
	for j in range(len(kappas)):
		kappa = kappas[j]
		kappa_p = np.sqrt(gamma*gamma + kappa*kappa)
		P_1     = ((F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p))) )

		#analytic results
		ux0 = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)

		ux2 = 0.5*(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))
		uy1 = 0.5*(kappa*P_1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))
		ux1 = 0.5*((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p*xi)/np.cosh(kappa_p))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))


		U[i, j] = 0.5*sci.trapz(np.real(ux0*np.conjugate(ux0)), xi)
		U[i, j] += 0.5*eps*eps* ( sci.trapz(   2*np.real( ux1*(np.conjugate(ux1))/4 +ux0*(np.conjugate(ux1)/2 + np.conjugate(ux2)) ), xi) + sci.trapz(   np.real(uy1*np.conjugate(uy1)), xi))
		U[i, j] = np.sqrt(U[i, j])


for i in range(len(taus)):
	plt.plot(kappas, U[i, :])
plt.show()
"""
for i in range(len(taus)):
	for j in range(len(Ds)):
		for f in range(len(F0s)):
			interpoo = scp.interp1d(kappas, D_parallels[:, i, f, j], "cubic")(kapppa_cont)
			if len(scs.argrelmax(interpoo)[0]) == 1:
				kappa_res[i, j, f] = kapppa_cont[scs.argrelmax(interpoo)]
				plt.plot(kapppa_cont, interpoo)
				plt.plot(kappa_res[i, j, f], interpoo[scs.argrelmax(interpoo)], "ro")
				plt.plot(kappas, D_parallels[:, i, f, j], "o", markersize=3)
plt.show()


marker = ["o", "+", "s", "^", "x"]
"""
plt.figure(1)
for i in range(len(taus)):
	for f in range(len(F0s)):
		indx = np.where(kappa_res[i, :, f] != 0)[0]
		if len(indx) > 0:
			plt.plot(np.sqrt(omegas[i]/Ds[indx]), kappa_res[i, indx, f], marker[f], markersize=5, color="C"+str(f), label=r"$\tau=%3.2f, F_0=%3.2f$" % (taus[i], F0s[f]))
			plt.plot(np.sqrt(omegas[i]/Ds[indx]), kappa_res[i, indx, f], color="C"+str(f), markersize=3)
	plt.legend(loc="best", fontsize=8)
plt.show()
"""
"""
#plt.xscale("log")
plt.xlabel(r"Diffusive Womersley number $\sqrt{\omega a^2/D_m}$", fontsize=8)
plt.ylabel(r"Resonance Wave number $\kappa_{res}$", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/resonance_wavelength.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


plt.figure(1)
for i in range(len(Ds)):
	for f in range(len(F0s)):
		indx = np.where(kappa_res[:, i, f] != 0)[0]
		if len(indx) > 0:
			plt.plot(F0s[f]/omegas[indx], 1/kappa_res[indx, i, f], marker[f], markersize=5, color="C"+str(f), label=r"$D_m=%3.2f, F_0=%3.2f$" % (Ds[i], F0s[f]))
			plt.plot(F0s[f]/omegas[indx], 1/kappa_res[indx, i, f], color="C"+str(f), markersize=3)
	plt.legend(loc="best", fontsize=8)
	plt.title(str(i))
	plt.show()
"""

plt.figure(1)
for i in range(len(Ds)):
	for f in range(len(F0s)):
		indx = np.where(kappa_res[:, i, f] != 0)[0]
		if len(indx) > 0:
			plt.plot(U[indx, f]*taus[indx], 2*np.pi/kappa_res[indx, i, f], marker[f], markersize=4, color="C"+str(i))
			plt.plot(U[indx, f]*taus[indx], 2*np.pi/kappa_res[indx, i, f], color="C"+str(i), markersize=3)
			print(np.mean(np.sqrt(Ds[i]*taus[indx])))
			#plt.title(str(Ds[i]))
			if len(indx) > 2:
				m, c, delta_m, delta_c = linear_regresion(F0s[f]/omegas[indx], 2*np.pi/kappa_res[indx, i, f])
				#print(m, delta_m)
#plt.show()

for i in range(len(F0s)):
	plt.plot(0, 0, marker[i], color="k", label=r"$F_0=%2.1f$" % (F0s[i]))
for i in range(len(Ds)):
	plt.plot(0, 0, "-", color="C"+str(i), label=r"$D_m=10^{%2.1f}$" % (Ds_logval[i]))
	#plt.plot(0, 0, "-", color="C"+str(i), label=r"$\rho=%2.1f$" % (np.sqrt(Ds[i])))

plt.axis([1.25, 30, 3.0, 19.75])
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Inverse Strouhal number $U\tau/a$", fontsize=8)
plt.ylabel(r"Resonance Wave length $\lambda_{res}$ [$a$]", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/resonance_wavelength_largerun.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()





plt.figure(1)
for i in range(len(Ds)):
	for f in range(len(F0s)):
		indx = np.where(kappa_res[:, i, f] != 0)[0]
		if len(indx) > 0:
			Lambda = 2*np.pi/kappa_res[indx, i, f]
			tau    = 2*np.pi/omegas[indx]
			x = F0s[f]*tau/(Lambda)
			plt.plot(x, np.ones(len(x)), marker[f], markersize=4, color="C"+str(i))
			plt.plot(x, np.ones(len(x)), color="C"+str(i), markersize=3)
			if len(indx) > 2:
				m, c, delta_m, delta_c = linear_regresion(F0s[f]/omegas[indx], 2*np.pi/kappa_res[indx, i, f])
				print(m, delta_m)
for i in range(len(F0s)):
	plt.plot(0, 0, marker[i], color="k", label=r"$F_0=%2.1f$" % (F0s[i]))
for i in range(len(Ds)):
	plt.plot(0, 0, "-", color="C"+str(i), label=r"$D_m=10^{%2.1f}$" % (Ds_logval[i]))

plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Inverse Strouhal number $U/(a\omega)$", fontsize=8)
plt.ylabel(r"Resonance Wave length $\lambda_{res}$", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/resonance_wavelength_largerun_2.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()







plt.figure(1)
for i in range(len(taus)):
	for f in range(len(F0s)):
		indx = np.where(kappa_res[i, :, f] != 0)[0]
		if len(indx) > 0:
			plt.plot(Ds[indx], kappa_res[i, indx, f], marker[f], color="C"+str(i), markersize=5)
			plt.plot(Ds[indx], kappa_res[i, indx, f], color="C"+str(i))


#plt.xscale("log")
plt.xlabel(r"Diffusive Womersley number $\sqrt{\omega a^2/D_m}$", fontsize=8)
plt.ylabel(r"Resonance Wave number $\kappa_{res}$", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/resonance_wavelength.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()




"""
for j in range(len(Ds)):
	for i in range(len(taus)):
		plt.plot(U[i, :], D_parallels[:, i, j], label=r"$D=%3.2f$" % Ds[j])
plt.legend(loc="best", fontsize=8)
plt.show()
"""
"""
for i in range(len(taus)-1):
	i += 1
	idxs = (np.where(kappa_res[i, :] != 0)[0])
	plt.plot(np.sqrt(omegas[i]/Ds[idxs]), kappa_res[i, idxs], label=r"$\omega=%3.2f$" % omegas[i], color="C"+str(i-1), markersize=3)
	#plt.plot(np.sqrt(omegas[i]/Ds[idxs]), np.sqrt(omegas[i]/(2*Ds[idxs])), color="C"+str(i-1))
plt.legend(loc="best", fontsize=8)
plt.xscale("log")
#plt.yscale("log")
plt.xlabel("D")
plt.ylabel("kappa res")
plt.show()

"""
"""
#kappa_res = taus, Ds
#U         = taus, kappas
U_res = np.zeros((len(taus), len(Ds)))

for i in range(len(taus)):
	interpoo = scp.interp1d(kappas, U[i, :], "cubic")(kapppa_cont)
	for j in range(len(Ds)):
		U_res[i, j] = interpoo[np.argmin(abs(kapppa_cont-kappa_res[i, j]))]
		plt.plot(interpoo[np.argmin(abs(kapppa_cont-kappa_res[i, j]))], kappa_res[i, j],  "o", color="C"+str(i))

plt.legend(loc="best", fontsize=8)
plt.xlabel("U")
plt.ylabel("kappa res")
plt.show()



for i in range(len(taus)):
	#i += 1
	idxs = (np.where(kappa_res[i, :] != 0)[0])
	if len(idxs) != 0:
		plt.plot(U_res[i, idxs], kappa_res[i, idxs], alpha=0.5, label=r"$\omega=%3.2f$" % omegas[i])

plt.legend(loc="best")
plt.show()


"""