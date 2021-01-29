import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

T = 3
omega = 2*np.pi/T
epsilon = 0.0
xi      = np.linspace(-1, 1, int(1e5))
datafiles = 500

nus = np.array([1.5, 3, 5])
datas = np.load("data/tdatas_zeroeps.npy") 
U = np.zeros(len(nus))
D_ana = np.zeros(len(nus))
num_D_para = np.zeros(len(nus))
Dm = 1

for i in range(len(nus)):
	Sc = nus[i]#/Dm
	F0 = 12/(nus[i])
	gamma   = np.sqrt(1j*omega/Sc)
	gamma_c = np.conj(gamma)
	a       = np.sqrt(1j*omega)
	a_c     = np.conj(a)
	rho = a 
	rho_c = a_c

	U[i] = np.sqrt(sci.trapz(datas[i, -datafiles:, 4], datas[i, -datafiles:, 0])/T)
	"""
	plt.figure(1)
	plt.plot(datas[i, :, 0],datas[i, :, 4] )
	plt.plot(datas[i, -datafiles:, 0], datas[i, -datafiles:, 4], "--" )

	plt.figure(2)
	plt.plot(datas[i, :, 0],datas[i, :, 8] )
	plt.plot(datas[i, -datafiles:, 0],datas[i, -datafiles:, 8], "--" )
	"""
	Pe = 1#U[i]/Dm
	factor   = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
	D_ana[i] = (1 + np.real(factor * 0.5 * sci.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)))
	num_D_para[i] = sci.trapz(datas[i, -datafiles:, 8], datas[i, -datafiles:, 0])/T

plt.show()

plt.plot(nus, num_D_para, "ko", label="numerisk")
plt.plot(nus, num_D_para, "k-")
plt.plot(nus, D_ana, "ro", label="analytisk")
plt.plot(nus, D_ana, "r-")
plt.legend(loc="best")
plt.xlabel("viskositet")
plt.ylabel("D_eff")
plt.show()

plt.plot(nus, (D_ana-1)/(num_D_para-1))
plt.show()