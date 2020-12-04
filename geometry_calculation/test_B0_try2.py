import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci 

Pe = 1
F0 = 3
pi = np.pi
N = int(1000)
xi = np.linspace(-1, 1, int(1e4))

Omega   = np.logspace(-3, 3, 60)
Schmidt = np.logspace(-3, 3, 60)
#paramters = np.zeros((2, len(Omega)))
#paramters[0, :] = Omega
#paramters[1, :] = Schmidt

D_eff   = np.zeros((len(Omega), len(Schmidt)), dtype="complex")

for s in range(len(Schmidt)):
	Sc = Schmidt[s]
	for o in range(len(Omega)):
		omega = Omega[o]

		gamma = np.sqrt(1j*omega/Sc)
		gamma_c = np.conj(gamma)
		a = np.sqrt(1j*omega)
		a_c = np.sqrt(-1j*omega)

		factor = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
		integral = 0.5*sci.trapz(np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)
		D_eff[o, s] = factor*integral#*105/2


plt.plot(Omega, np.real(D_eff[:, 40]))
plt.plot(Omega, np.imag(D_eff[:, 40]))

np.save("data_test/D_eff_2", D_eff)
prev = np.load("data_test/D_eff.npy")
plt.plot(Omega, np.real(prev[:, 40]), "--")
plt.xscale("log")
#plt.yscale("log")
plt.show()




