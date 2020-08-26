import numpy as np 
import matplotlib.pyplot as plt 

omega = 1
Sc = 1
Pe = 1
F0 = 3/2
pi = np.pi
N = int(100)
xi = np.linspace(-1, 1, int(1e4))

summ_before = 0
summ_after = 0
partial_summ = 0
for n_prime in range(N):
	n       = 2*n_prime+1
	alpha_n = n*n*pi*pi*Sc/4
	a_n     = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
	prefactor   = Pe*(a_n*n*pi/2)/(1j*omega + n*n*pi*pi/4)
	first_term  = np.sin(n*pi*xi/2)
	second_term = ((-1)**(0.5*(n-1)))*np.sinh(np.sqrt(1j*omega)*xi)/(np.sinh(np.sqrt(1j*omega)))

	summ_before += prefactor*(first_term-second_term)

total = 0
term1 = 0
term2 = 0
m_sum = 0

for n_prime in range(N):
	n       = 2*n_prime+1
	alpha_n = n*n*pi*pi*Sc/4
	a_n     = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
	prefactor   = Pe*(a_n*n*pi/2)/(1j*omega + n*n*pi*pi/4)

	term1 += prefactor*np.sin(n*pi*xi/2)

for m_prime in range(N):
	m = 2*m_prime+1
	m_sum += ((-1)**(0.5*(m-1)))*np.sin(m*pi*xi/2)/(1j*omega + m*m*pi*pi/4)

for n_prime in range(N):
	n = 2*n_prime+1
	alpha_n = n*n*pi*pi*Sc/4
	a_n     = Pe*4*((-1)**(0.5*(n-1)))*Sc*F0*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n+omega*omega))
	prefactor   = Pe*(a_n*n*pi/2)/(1j*omega + n*n*pi*pi/4)

	term2 += prefactor*(-1)**(0.5*(n-1))*2*np.sqrt(1j*omega)*m_sum/np.sinh(np.sqrt(1j*omega))
	

plt.plot(xi, term1 + term2)
plt.plot(xi, summ_before, "--")
plt.show()

"""
plt.plot(xi, summ)
plt.plot(xi, np.sinh(np.sqrt(1j*omega)*xi), "--")
plt.show()
"""