import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci 


Omega = np.logspace(-4, 4, 70)
Sc = 1
Pe = 1
F0 = 3
pi = np.pi
M = int(30)

a_m_squared = np.zeros(M, dtype="complex")
D_eff = np.zeros(len(Omega), dtype="complex")

for o in range(len(Omega)):
	omega = Omega[o]
	a_m_squared = np.zeros(M, dtype="complex")
	for m in range(1, M):
		for k_prime in range(M):
			k = 2*k_prime + 1
			alpha_k = k*k*pi*pi*Sc/4

			for n_prime in range(M):
				n = 2*n_prime + 1
				alpha_n = n*n*pi*pi*Sc/4
				a_m_squared[m] += (alpha_n-1j*omega)*(alpha_k+1j*omega)/((alpha_n*alpha_n + omega*omega)*(alpha_k*alpha_k + omega*omega)*(m*m-n*n/4)*(m*m-k*k/4))

	g_tilde = 0
	for m in range(M):
		g_tilde += m*m*a_m_squared[m]/(m*m*m*m*pi*pi*pi*pi + omega*omega)

	g_tilde *= 16*Pe*Pe*Sc*Sc*F0*F0/(pi*pi)
	D_eff[o] = 0.5*g_tilde #0.5 due to decomposisition of u in eq 2.6 in Bruus

np.save("data/D_eff_Bruus", D_eff)
print("saved")
plt.plot(Omega, np.real(D_eff), "o")
plt.xscale("log")
plt.show()
