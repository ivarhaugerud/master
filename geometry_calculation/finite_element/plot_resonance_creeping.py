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
D_parallels = np.load("data/D_para_vary_kappa_D_omega_nu1000_F0500.npy")
taus      = np.logspace(-2, 2, 10)
omegas    = 2*np.pi/taus
nu        = 1000
F0        = 1000/nu
Ds        = np.logspace(-2, 1, 10)
kappas   = np.array([0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6, 4.4, 5.0])
kapppa_cont = np.linspace(min(kappas), max(kappas), int(1e4))
kappa_res = np.zeros((len(taus), len(Ds)))

#run 2
D_parallels2 = np.load("data/D_para_vary_kappa_D_tau20_nu1000_F01000_more.npy")
taus2      = np.array([1, 10, 20, 30, 40])
omegas2    = 2*np.pi/taus
nu2        = 1000
F02        = 500/nu
Ds2        = np.logspace(-1.5, 0, 10)
kappas2   = np.arange(0.1, 1.801, 0.1)

kapppa_cont2 = np.linspace(min(kappas2), max(kappas2), int(1e4))
kappa_res2 = np.zeros((len(taus), len(Ds)))


U = np.zeros((len(taus), len(kappas)))
Re = np.zeros(np.shape(U))
xi = np.linspace(-1, 1, int(1e4))
eps = 0.3

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


U2 = np.zeros((len(taus2), len(kappas2)))
Re2 = np.zeros(np.shape(U2))
F0 = F02
for i in range(len(taus2)):
	omega = 2*np.pi/taus2[i]
	gamma   = np.sqrt(1j*omega/nu)
	for j in range(len(kappas2)):
		kappa = kappas2[j]
		kappa_p = np.sqrt(gamma*gamma + kappa*kappa)
		P_1     = ((F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p))) )

		#analytic results
		ux0 = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)

		ux2 = 0.5*(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))
		uy1 = 0.5*(kappa*P_1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))
		ux1 = 0.5*((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p*xi)/np.cosh(kappa_p))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))


		U2[i, j] = 0.5*sci.trapz(np.real(ux0*np.conjugate(ux0)), xi)
		U2[i, j] += 0.5*eps*eps* ( sci.trapz(   2*np.real( ux1*(np.conjugate(ux1))/4 +ux0*(np.conjugate(ux1)/2 + np.conjugate(ux2)) ), xi) + sci.trapz(   np.real(uy1*np.conjugate(uy1)), xi))
		U2[i, j] = np.sqrt(U2[i, j])

for i in range(len(taus2)):
	plt.plot(kappas2, U2[i, :])
plt.show()
for i in range(len(taus)):
	plt.plot(kappas, U[i, :])
plt.show()

for i in range(len(taus)):
	for j in range(len(Ds)):
		interpoo = scp.interp1d(kappas, D_parallels[:, i, j], "cubic")(kapppa_cont)
		if len(scs.argrelmax(interpoo)[0]) == 1:
			kappa_res[i, j] = kapppa_cont[scs.argrelmax(interpoo)]

for i in range(len(taus2)):
	for j in range(len(Ds2)):
		interpoo = scp.interp1d(kappas2, D_parallels2[:, i, j], "cubic")(kapppa_cont2)
		if len(scs.argrelmax(interpoo)[0]) == 1:
			kappa_res2[i, j] = kapppa_cont2[scs.argrelmax(interpoo)]





plt.figure(1)
for i in range(len(taus2)):
	indx = np.where(kappa_res2[i, :] != 0)[0]
	if len(indx) > 0:
		print(Ds2[indx], kappa_res2[i, indx])
		plt.plot(np.sqrt(omegas2[i]/Ds2[indx]), kappa_res2[i, indx], "x", color="C"+str(i), markersize=3, label=r"$\tau=%3.2f, F_0/\nu=1$" % taus2[i])


for i in range(len(taus)):
	indx = np.where(kappa_res[i, :] != 0)[0]
	if len(indx) > 0:
		plt.plot(np.sqrt(omegas[i]/Ds[indx]), kappa_res[i, indx], "o", color="C"+str(i), markersize=3, label=r"$\tau=%3.2f, F_0/\nu=0.5$" % taus[i])

plt.xscale("log")
plt.xlabel(r"Diffusive Womersley number $\sqrt{\omega a^2/D_m}$", fontsize=8)
plt.ylabel(r"Resonance Wave number $\kappa_{res}$", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/resonance_wavelength.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

F0        = 1000/nu


plt.figure(1)
for i in range(len(taus2)):
	indx = np.where(kappa_res2[i, :] != 0)[0]
	if len(indx) > 0:
		print(Ds2[indx], kappa_res2[i, indx])
		plt.plot(Ds2[indx], kappa_res2[i, indx]*F02, "x", color="C"+str(i), markersize=3, label=r"$\tau=%3.2f, F_0/\nu=1$" % taus2[i])


for i in range(len(taus)):
	indx = np.where(kappa_res[i, :] != 0)[0]
	if len(indx) > 0:
		plt.plot(Ds2[indx], kappa_res[i, indx]*F0, "o", color="C"+str(i), markersize=3, label=r"$\tau=%3.2f, F_0/\nu=0.5$" % taus[i])

plt.xscale("log")
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