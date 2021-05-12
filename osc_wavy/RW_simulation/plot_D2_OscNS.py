import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

T = 3
omega = 2*np.pi/T
epsilon = 0.25
xi      = np.linspace(-1, 1, int(1e5))
datafiles = 500
nu = 3.6	
kappas  = np.array([0.5 , 0.9, 1.1, 1.3, 1.7, 2.1])

datas = np.load("data/tdatas_zeroeps_Pe6.npy") 
U = np.zeros(len(kappas))
D_ana = np.zeros(len(kappas))
num_D_para = np.zeros(len(kappas))
Dm = 1

for i in range(len(kappas)):
	Sc = nu#/Dm
	F0 = 12/(nu)
	gamma   = np.sqrt(1j*omega/Sc)
	gamma_c = np.conj(gamma)
	a       = np.sqrt(1j*omega)
	a_c     = np.conj(a)
	rho = a 
	rho_c = a_c

	U[i] = np.sqrt(sci.trapz(datas[i, -datafiles:, 4], datas[i, -datafiles:, 0])/T)
	plt.figure(1)
	plt.plot(datas[i, :, 0],datas[i, :, 4] )
	plt.plot(datas[i, -datafiles:, 0], datas[i, -datafiles:, 4], "--" )

	plt.figure(2)
	plt.plot(datas[i, :, 0],datas[i, :, 8] )
	plt.plot(datas[i, -datafiles:, 0],datas[i, -datafiles:, 8], "--" )
	Pe = 6#U[i]/Dm
	factor   = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
	D_ana[i] = (1 + np.real(factor * 0.5 * sci.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)))
	num_D_para[i] = sci.trapz(datas[i, -datafiles:, 8], datas[i, -datafiles:, 0])/T

plt.show()

plt.plot(kappas, num_D_para, "ko", label="numerisk")
plt.plot(kappas, num_D_para, "k-")
plt.plot(kappas, D_ana, "ro", label="analytisk")
plt.plot(kappas, D_ana, "r-")
plt.legend(loc="best")
plt.xlabel("viskositet")
plt.ylabel("D_eff")
#plt.show()




from numpy import *
import numpy as np  
import matplotlib.pyplot as plt 
import os
import scipy.integrate as sci

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'

xi = np.linspace(-1, 1, int(1e4))

rho = np.sqrt(1j*omega)
gamma = np.sqrt(1j*omega/Sc)

#kappa = np.linspace(0.01, 5, int(1e3))
D = np.zeros(len(kappas))

for i in range(len(kappas)):
	kappa = kappas[i]
	rho = np.sqrt(1j*omega)
	gamma = np.sqrt(1j*omega/Sc)
	kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
	P_1  = (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p)))

	B0_deriv       = Pe*F0*(sinh(rho*xi)/sinh(rho)-xi)/(2*rho*rho)
	B0_deriv_deriv = Pe*F0*(rho*cosh(rho*xi)/sinh(rho)-1)/(2*rho*rho)
	B_plus         = F0*(1-xi*xi - 2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(4*rho*rho) + Pe*F0*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)*(1+rho/tanh(rho))/(2*rho*sinh(rho)))/(rho*rho)

	sol1 = xi*xi*cosh(rho*xi)*F0/(sinh(rho)*16*rho*rho*rho)*(kappa*kappa+Pe*(2*rho*rho-kappa*kappa)-Pe*rho*kappa*kappa/tanh(rho)) + xi*sinh(rho*xi)/sinh(rho)*F0/(4*rho*rho)*(1-Pe/2-Pe*rho/tanh(rho) - (kappa*kappa+Pe*(2*rho*rho-kappa*kappa)-Pe*rho*kappa*kappa/(tanh(rho)))/(4*rho*rho)) - F0*kappa*kappa*(2*Pe-1)*(xi*xi+2/(rho*rho))/(4*rho*rho*rho*rho) - F0*(3*Pe/2-1)/(2*rho*rho*rho*rho)
	sol2 = Pe*F0*gamma*gamma*(24/(rho*rho*rho*rho) + 12*xi*xi/(rho*rho) - 4/(rho*rho) + xi*xi*xi*xi -2*xi*xi + 1)/(16*rho*rho) + Pe*(2/(rho*rho) + xi*xi - 1)*(F0 + P_1*kappa*(2-gamma*gamma))/(8*rho*rho) - Pe*(F0*gamma*gamma/15 - P_1*kappa*(1-gamma*gamma/2)/3)/(2*rho*rho)
	B2 = sol1+sol2 -cosh(rho*xi)*(-F0*Pe*kappa**2*rho**2 + 2*F0*Pe*rho**4*tanh(rho)**2 - 4*F0*Pe*rho**4 + F0*kappa**2*rho*tanh(rho) - F0*rho**3*(Pe*kappa**2 + 4*Pe - 4)*tanh(rho) + F0*(24*Pe*gamma**2 - 15*Pe*kappa**2 + 7*kappa**2)*tanh(rho)**2 + rho**2*(-F0*Pe*kappa**2 + F0*kappa**2 + 4*F0 - 4*P_1*Pe*gamma**2*kappa + 8*P_1*Pe*kappa)*tanh(rho)**2)/(sinh(rho)*16*rho**5*tanh(rho)**2)


	D_para =  B0_deriv*np.conj(B0_deriv)*(1+kappa*kappa*xi*xi)/2 + (np.conj(B_plus)*B_plus*kappa*kappa/2 + np.gradient(np.conj(B_plus), xi)*np.gradient(B_plus, xi)/2) #- np.conj(B0_deriv)*(np.gradient(B_plus, xi)+kappa*kappa*xi*B_plus -2*np.gradient(B2, xi))
	D_para += np.conj(D_para)
	D_para /= 2

	#plt.plot(xi, D_para)
	D[i] = sci.trapz(D_para, xi)/2
#plt.show()

plt.plot(kappas, D+D_ana[0])
#np.save("data/analytic_Dpara", np.array([kappas, D]))
#plt.xscale("log")
#plt.yscale("log")

plt.show()
