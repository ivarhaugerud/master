import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci 

omega = 1
Sc = 1
Pe = 1
F0 = 3/2
pi = np.pi
N = int(500)

t = np.linspace(0, 2*np.pi/omega, 1e3)
xi = np.linspace(-1, 1, 1e4)
B_0 = np.zeros((len(t), len(xi)))
B_0_squared_averaged = np.zeros(len(t))

summ = 0
"""
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

D_eff = np.mean(B_0_squared_averaged)
print("Using B_0, squaring then integrating: ", D_eff)
"""
"""
term1 = 0
term2 = 0
term3 = 0
N = 500
#
s2w = np.sqrt(2*omega)

for n_prime in range(N):
	n = 2*n_prime + 1
	alpha_n = Sc*n*n*pi*pi/4
	term1 += 1/((alpha_n*alpha_n + omega*omega)*(alpha_n*alpha_n/(Sc*Sc) + omega*omega))

	for m_prime in range(N):
		m = 2*m_prime + 1
		alpha_m = Sc*m*m*pi*pi/4

		term2 += ((-1)**(0.5*(m+n)-1))*(alpha_n*alpha_m+omega*omega)*np.sqrt(omega/2)*(m*m + n*n)/((alpha_n*alpha_n+omega*omega)*(alpha_m*alpha_m+omega*omega)*(alpha_n*alpha_n/(Sc*Sc)+omega*omega)*(alpha_m*alpha_m/(Sc*Sc)+omega*omega))
		term3 += ((-1)**(0.5*(m+n)-1))*np.real((alpha_n*alpha_m+omega*omega)/((alpha_n*alpha_n+omega*omega)*(alpha_m*alpha_m+omega*omega)*(alpha_m/Sc + 1j*omega)*(alpha_n/Sc-1j*omega)))

answer = 0.5*term1 + (np.sin(s2w)+np.sinh(s2w))*pi*pi*term2/(4*(np.cosh(s2w)-np.cos(s2w))) + (np.sin(s2w)-np.sinh(s2w))*term3/(8*s2w*(np.cos(s2w)-np.cosh(s2w)))
answer *= Pe*Pe*Sc*Sc*F0*F0
ans = Pe*Pe*Sc*Sc*F0*F0*np.array([0.5*term1, (np.sin(s2w)+np.sinh(s2w))*pi*pi*term2/(4*(np.cosh(s2w)-np.cos(s2w))), (np.sin(s2w)-np.sinh(s2w))*term3/(8*s2w*(np.cos(s2w)-np.cosh(s2w)))])
print("My calculation with analytic integratin and squaring: ", answer)
"""
"""
Omega = np.logspace(-3, 3, 20)
D_eff = np.zeros(len(Omega))
D_eff_bruus = np.zeros(len(Omega))

n = 0
alpha_n = 0 

term1 = 0
term2 = 0
term3 = 0
N = 1
Sc = 1
Pe = 1 
F0 = 3/2
pi = np.pi 


for o in range(len(Omega)):
	n = 0
	alpha_n = 0 

	term1 = 0
	term2 = 0
	term3 = 0
	answer= 0

	omega = Omega[o]
	s2w = np.sqrt(2*omega)
	for n_prime in range(N):
		n = 2*n_prime + 1
		alpha_n = Sc*n*n*pi*pi/4
		term1 += 1/((alpha_n*alpha_n + omega*omega)*(alpha_n*alpha_n/(Sc*Sc) + omega*omega))

		for m_prime in range(N):
			m = 2*m_prime + 1
			alpha_m = Sc*m*m*pi*pi/4

			term2 += (alpha_n*alpha_m+omega*omega)*np.sqrt(omega/2)*(m*m + n*n)/((alpha_n*alpha_n+omega*omega)*(alpha_m*alpha_m+omega*omega)*(alpha_n*alpha_n/(Sc*Sc)+omega*omega)*(alpha_m*alpha_m/(Sc*Sc)+omega*omega))
			term3 += np.real((alpha_n*alpha_m+omega*omega)/((alpha_n*alpha_n+omega*omega)*(alpha_m*alpha_m+omega*omega)*(alpha_m/Sc + 1j*omega)*(alpha_n/Sc-1j*omega)))


	answer = term1/2 + (np.sin(s2w)+np.sinh(s2w))*pi*pi*term2/(4*(np.cosh(s2w)-np.cos(s2w))) + (np.sin(s2w)-np.sinh(s2w))*term3/(8*s2w*(np.cos(s2w)-np.cosh(s2w)))
	answer *= Pe*Pe*Sc*Sc*F0*F0 # I dont have to divide by 2 pi omega since this would just cancel what I get from integrating
	D_eff[o] = answer
	D_eff_bruus[o] = term1*8

plt.plot(Omega, D_eff, "o", label="Meg")
plt.plot(Omega, D_eff_bruus, "o", label="Bruus")
plt.xscale("log")
plt.legend(loc="best")
plt.show()

N = int(500)
t = np.linspace(0, 2*np.pi/omega, 2000)

xi = np.linspace(-1, 1, 1e4)
Omega = np.logspace(-3, 3, 20)
D_eff = np.zeros(len(Omega))
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

	D_eff[o] = sci.trapz(B_0_squared_averaged, t)#/(2*pi/omega)
"""
#plt.plot(Omega, D_eff, "--")
#plt.xscale("log")
#plt.show()

