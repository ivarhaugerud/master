import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci 


Omega = np.logspace(-3, 3, 10)
Sc = 1
Pe = 1
F0 = 1
pi = np.pi
M = int(100)

a_m_squared = np.zeros(M, dtype="complex")
D_eff = np.zeros(len(Omega), dtype="complex")

for o in range(len(Omega)):
	omega = Omega[o]
	a_m_squared = np.zeros(M, dtype="complex")
	for m in range(M):
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
	D_eff[o] = g_tilde

plt.plot(Omega, D_eff, "o")
plt.xscale("log")
plt.show()

N = int(500)
t = np.linspace(0, 2*np.pi/omega, int(1e4))

xi = np.linspace(-1, 1, 1e4)
#Omega = np.logspace(-3, 3, 20)
D_eff_2 = np.zeros(len(Omega))
B_0_squared_averaged = np.zeros(len(t))
B_0 = np.zeros((len(t), len(xi)))

for o in  range(len(Omega)):
	omega = Omega[o]
	summ = 0
	for n_prime in range(N):
		n       = 2*n_prime+1
		alpha_n = n*n*pi*pi*Sc/4
		a_n     = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
		prefactor   = Pe*(a_n*n*pi/2)/(1j*omega + n*n*pi*pi/4)
		first_term  = np.sin(n*pi*xi/2)
		second_term = ((-1)**(0.5*(n-1)))*np.sinh(np.sqrt(1j*omega)*xi)/(np.sinh(np.sqrt(1j*omega)))

		summ += prefactor*(first_term-second_term)

	for i in range(len(t)):
		B_0[i,:] = np.real(summ*np.exp(1j*omega*t[i]))
		B_0_squared_averaged[i] = 0.5*sci.trapz(B_0[i,:]*B_0[i,:], xi)

	D_eff_2[o] = sci.trapz(B_0_squared_averaged, t)#/(2*pi/omega)
	
plt.plot(Omega, D_eff_2, "o", label="Meg")
plt.plot(Omega, D_eff, "o", label="Bruus")
plt.legend(loc="best")
plt.xscale("log")
plt.show()




plt.plot(Omega, 300*D_eff_2/D_eff, "o", label="Meg")
plt.xscale("log")
plt.show()