import numpy as np 
import matplotlib.pyplot as plt 
from scipy import integrate

nu = 1.2
Sc = nu 
F0 = 6
pi = np.pi 
Pe = 1

xi = np.linspace(-1, 1,  int(1e4))
t  = np.linspace(0, 1, int(500))

B0_grad = np.zeros(len(xi), dtype="complex")
B0_grad_2 = np.zeros(len(xi), dtype="complex")
B0_grad_real = np.zeros((len(t), len(xi)))
B0_grad_2_real = np.zeros((len(t), len(xi)))
Omega = np.logspace(-3, 3, 12)
D_para = np.zeros(len(Omega))
D_para_real = np.zeros((len(Omega), len(t)))
D_ana  = np.zeros(len(Omega))

for o in range(len(Omega)):
	omega = Omega[o]
	tau = 2*np.pi/omega
	gamma = np.sqrt(1j*omega/Sc)
	rho   = np.sqrt(1j*omega)
	gamma_c = np.conj(gamma)
	rho_c   = np.conj(rho)
	t  = np.linspace(0, 2*np.pi/omega, int(500))

	B0_grad[:]      = Pe*F0*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))
	
	for j in range(len(t)):
		B0_grad_real[j, :] = np.real(np.exp(1j*omega*t[j])*Pe*F0*np.tanh(gamma)/(gamma*(rho*rho-gamma*gamma))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma)))

	B0_grad_2[:]      = 2*B0_grad[:]*np.conj(B0_grad[:])
	B0_grad_2_real = B0_grad_real*B0_grad_real

	D_para[o]      = 105/2 * integrate.trapz(np.real(B0_grad_2), xi)/2
	D_para_real[o] = 105/2 * integrate.trapz(B0_grad_2_real, xi)/2
	D_para_real[o, 0] = integrate.trapz(D_para_real[o,:], t)/tau
	#D_ana[o]       = 2*105/2 * Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(rho*rho-gamma*gamma)*(rho_c*rho_c-gamma_c*gamma_c))*integrate.trapz((np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))*(np.sinh(rho_c*xi)/np.sinh(rho_c) - np.sinh(gamma_c*xi)/np.sinh(gamma_c)), xi)/2
	D_ana[o] = 105/2*np.real(Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*omega*omega*omega*(Sc*Sc-1))*(1/(gamma*np.tanh(gamma)) + 1j/(gamma*np.tanh(gamma_c)) - 1/(rho*np.tanh(rho)) - 1j/(rho*np.tanh(rho_c))))


plt.plot(Omega, D_para)
plt.plot(Omega, D_ana, "--")
plt.plot(Omega, D_para_real[:,0], "o")
plt.xscale("log")
plt.show()


#VELOCITY FIELD
"""
u_x_complex = np.zeros((len(t), len(xi)), dtype="complex")
u_x_real    = np.zeros((len(t), len(xi)))
N = 100
for i in range(len(t)):
	u_x_complex[i, :] = (F0/(2*gamma*gamma))*(1-np.cosh(gamma*xi)/np.cosh(gamma))*np.exp(1j*omega*t[i])
	for j in range(N):
		n = 2*j+1
		alpha_n = n*n*pi*pi*Sc/4
		u_x_real[i, :] += (4*(-1)**(int((n-1)/2))/(n*pi))*np.cos(n*pi*xi/2)*(alpha_n*np.cos(omega*t[i]) + omega*np.sin(omega*t[i]))/(alpha_n*alpha_n+omega*omega)
u_x_real *= Sc*F0

u_x_complex = u_x_complex + np.conj(u_x_complex)

plt.plot(xi, np.real(u_x_complex[50, :]))
plt.plot(xi, u_x_real[50, :], "--")
plt.plot(xi, np.real(u_x_complex[170, :]))
plt.plot(xi, u_x_real[170, :], "--")
plt.plot(xi, np.real(u_x_complex[300, :]))
plt.plot(xi, u_x_real[300, :], "--")
plt.plot(xi, np.real(u_x_complex[460, :]))
plt.plot(xi, u_x_real[460, :], "--")
plt.show()
"""