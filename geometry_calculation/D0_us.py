import numpy as np 
import matplotlib.pyplot as plt 

Pe = 1
F0 = 3
pi = np.pi
N = int(100)
xi = np.linspace(-1, 1, int(1e3))

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


		theta = 0
		for n_prime in range(N):
			n       = 2*n_prime+1
			alpha_n = n*n*pi*pi*Sc/4
			a_n     = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
			psi_n   = (a_n*n*pi/2)/(1j*omega + n*n*pi*pi/4)
			theta += psi_n*((-1)**(0.5*(n-1)))/(np.sinh(np.sqrt(1j*omega)))

		phi_n = np.zeros(N, dtype="complex")
		psi_n = np.zeros(N, dtype="complex")
		rho_n       = np.zeros(N, dtype="complex")

		for n_prime in range(N):
			n = 2*n_prime + 1
			phi_n[n_prime] = 2*np.sqrt(1j*omega)*((-1)**(0.5*(n-1)))*np.cosh(np.sqrt(1j*omega))/(1j*omega +pi*pi*n*n/4)

			alpha_n        = n*n*pi*pi*Sc/4
			a_n            = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
			psi_n[n_prime] = (a_n*n*pi/2)/(1j*omega + n*n*pi*pi/4)

			rho_n[n_prime] = theta*phi_n[n_prime] - psi_n[n_prime]

		"""
		###COMPLEX CONJUGATE
		theta_cc = 0
		for n_prime in range(N):
			n       = 2*n_prime+1
			alpha_n = n*n*pi*pi*Sc/4
			a_n     = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n+1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
			psi_n   = (a_n*n*pi/2)/(-1j*omega + n*n*pi*pi/4)
			theta_cc += psi_n*((-1)**(0.5*(n-1)))/(np.sinh(np.sqrt(-1j*omega)))


		phi_n_cc = np.zeros(N, dtype="complex")
		psi_n_cc = np.zeros(N, dtype="complex")
		rho_n_cc       = np.zeros(N, dtype="complex")

		for n_prime in range(N):
			n = 2*n_prime + 1
			phi_n_cc[n_prime] = 2*np.sqrt(-1j*omega)*((-1)**(0.5*(n-1)))*np.cosh(np.sqrt(-1j*omega))/(-1j*omega +pi*pi*n*n/4)

			alpha_n        = n*n*pi*pi*Sc/4
			a_n            = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n+1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
			psi_n_cc[n_prime] = (a_n*n*pi/2)/(-1j*omega + n*n*pi*pi/4)

			rho_n_cc[n_prime] = theta_cc*phi_n_cc[n_prime] - psi_n_cc[n_prime]
		"""
		for i in range(N):
			D_eff[o, s] += rho_n[i]*np.conj(rho_n[i])

		D_eff[o, s] *= 0.5*Pe*Pe

plt.plot(Omega, np.real(D_eff[:, 40]))
plt.plot(Omega, np.imag(D_eff[:, 40]))

plt.xscale("log")
#plt.yscale("log")
plt.show()


np.save("data_test/D_eff", D_eff)
np.save("data_test/D_eff_params_Omega_Schmidt", D_eff)


