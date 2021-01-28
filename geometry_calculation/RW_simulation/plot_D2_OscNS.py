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
	Pe = 1#U[i]/Dm
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
plt.show()