import numpy as np 
import matplotlib.pyplot as plt 

omega = 1
Sc = 1
Pe = 1
F0 = 1
pi = np.pi
N = int(1000)
xi = np.linspace(-1, 1, int(1e4))

Omega = np.logspace(-3, 3, 40)
D_eff = np.zeros(len(Omega), dtype="complex")

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

	# TEST IF COMPLEX CONJUGATION ACTUALLY WORKED
	"""
	plt.plot(np.real(rho_n))
	plt.plot(np.real(rho_n_cc), "--")
	plt.show()

	plt.plot(np.imag(rho_n))
	plt.plot(-np.imag(rho_n_cc), "--")
	plt.show()
	"""

	for i in range(N):
		D_eff[o] += rho_n[i]*rho_n_cc[i]

	D_eff[o] *= Pe*Pe*0.5

np.save("D_eff_Ivar", D_eff)
print("Saved")
plt.plot(Omega, np.real(D_eff), "o")
print(np.mean(np.imag(D_eff)))

plt.xscale("log")
plt.show()




