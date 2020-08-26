import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci 

omega = 1
Sc = 1
Pe = 1
F0 = 3/2
pi = np.pi
N = int(500)
### Vendel og Bruus
"""
g_tilde_VB = 0
for n_tilde in range(N):
	n = 2*n_tilde+1
	alpha_n = n*n*pi*pi*Sc/4

	g_tilde_VB += 1/((alpha_n*alpha_n/(Sc*Sc)+omega*omega)*(alpha_n*alpha_n+omega*omega))


g_tilde_VB *= 8*Pe*Pe*Sc*Sc*F0*F0
print("Vendel og Bruus: ", g_tilde_VB)
###
t = np.linspace(0, 2*np.pi/omega, 1e3)
summ = 0
xi = np.linspace(-1, 1, 1e4)
B_0_squared_averaged = np.zeros(len(t))

B_0 = np.zeros((len(t), len(xi)))
summ = 0
test_2 = 0

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

		term2 += (alpha_n*alpha_m+omega*omega)*np.sqrt(omega/2)*(m*m + n*n)/((alpha_n*alpha_n+omega*omega)*(alpha_m*alpha_m+omega*omega)*(alpha_n*alpha_n/(Sc*Sc)+omega*omega)*(alpha_m*alpha_m/(Sc*Sc)+omega*omega))
		term3 += ((-1)**(0.5*(m+n)-1))*np.real((alpha_n*alpha_m+omega*omega)/((alpha_n*alpha_n+omega*omega)*(alpha_m*alpha_m+omega*omega)*(alpha_m/Sc + 1j*omega)*(alpha_n/Sc-1j*omega)))

#print(term1, term2, term3)
answer = 0.5*term1 + (np.sin(s2w)+np.sinh(s2w))*pi*pi*term2/(4*(np.cosh(s2w)-np.cos(s2w))) + (np.sin(s2w)-np.sinh(s2w))*term3/(8*s2w*(np.cos(s2w)-np.cosh(s2w)))
answer *= Pe*Pe*Sc*Sc*F0*F0
ans = Pe*Pe*Sc*Sc*F0*F0*np.array([0.5*term1, (np.sin(s2w)+np.sinh(s2w))*pi*pi*term2/(4*(np.cosh(s2w)-np.cos(s2w))), (np.sin(s2w)-np.sinh(s2w))*term3/(8*s2w*(np.cos(s2w)-np.cosh(s2w)))])
print(ans[0], ans[1], ans[2])
print(14*ans[0]+ans[1]+ans[2])
#print(8*0.5*term1, 8*(np.sin(s2w)+np.sinh(s2w))*pi*pi*term2/(4*(np.cosh(s2w)-np.cos(s2w))), 8*(np.sin(s2w)-np.sinh(s2w))*term3/(8*s2w*(np.cos(s2w)-np.cosh(s2w))))
print("My calculation with analytic integratin and squaring: ", answer)
"""

Omega = np.logspace(0, 0, 1)
D_eff = np.zeros(len(Omega))

n = 0
alpha_n = 0 

term1 = 0
term2 = 0
term3 = 0
N = 1
Sc = 1
Pe = 3/2 
F0 = 1
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

plt.plot(Omega, D_eff, "o")
print(D_eff)
plt.xscale("log")
plt.show()

N = int(250)
t = np.linspace(0, 2*np.pi/omega, 10000)
summ = 0
xi = np.linspace(-1, 1, 1e4)
Omega = np.logspace(-2, 2, 50)
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

	D_eff[o] = sci.trapz(B_0_squared_averaged, t)

plt.plot(Omega, D_eff*105/2, "o")
plt.xscale("log")
plt.show()


"""
alpha = omega*1*1/1
alpha_bar = alpha/2 # mean of alpha over cross-section
epsilon = F0*Sc
S = Sc
D_eff = 0

for k in range(1, 1000):
	a_k = 0#((-1)**(k+1))*np.sqrt(2)/(k*k*pi*pi)
	a_k_prime = 0#a_k/(k*k*pi*pi)
	b_k = ((-1)**k)*1j*np.sqrt(2*1j)*epsilon*np.tanh(1j*alpha_bar)/((np.sqrt(alpha)*(1j*alpha+k*k*pi*pi)))
	b_k_prime = b_k/(1j*alpha*S + k*k*pi*pi)
	#print(k, b_k, b_k_prime)
	D_eff += np.exp(2*1j*alpha*S*t)*b_k*b_k_prime + np.exp(1j*alpha*S*t)*(a_k*b_k_prime+a_k_prime*b_k)

D_eff *= Pe*Pe/2
#plt.plot(D_eff)
#plt.show()
print(np.mean(np.real(D_eff)), np.max(np.real(D_eff)))
"""
