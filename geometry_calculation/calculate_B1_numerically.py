import numpy as np 
import matplotlib.pyplot as plt 

k = np.arange(-10, 10, 1)
xi = np.linspace(-1, 1, int(1e3))
b_k = np.zeros((len(k), len(xi)))
dx = xi[1]-xi[0]
q = np.zeros((len(k), len(xi)), dtype="complex")

kappa = 2
Sc = 1.2 
omega = 5
F0 = 3
Pe = F0/Sc 

gamma = np.sqrt(1j*omega/Sc )
kappa_p = np.sqrt(gamma*gamma+kappa*kappa)
rho = np.sqrt(kappa*kappa+1j*omega)
P_1 = (F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))

ux0  = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)
ux1 = (P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p *xi)/np.cosh(kappa_p)) + (F0*np.tanh(gamma)/gamma)*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma))
uy1 = (kappa*P_1*np.sinh(kappa)/(gamma*gamma))*(np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))

B0 = Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma)))
B0_deriv = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))
B0_deriv_deriv = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma))

RHS_sin  = kappa*kappa*xi*B0_deriv - B0_deriv_deriv + Pe*ux1
RHS_cos  = (ux0*kappa*xi - uy1)*B0_deriv

q[np.argmin(abs(k+2)), :] = kappa*xi*np.conj(ux0)*np.conj(B0_deriv)/4 - np.conj(uy1)*np.conj(B0_deriv)/4
q[np.argmin(abs(k+1)), :] = kappa*kappa*xi*np.conj(B0_deriv)/2 - 2*np.conj(B0_deriv_deriv)/2 + Pe*np.conj(ux1)/2
q[np.argmin(abs(k-0)), :] = kappa*xi*(ux0*np.conj(B0_deriv) + np.conj(ux0)*B0_deriv) - uy1*np.conj(B0_deriv) - np.conj(uy1)*B0_deriv
q[np.argmin(abs(k-1)), :] = kappa*kappa*xi*B0_deriv/2 - 2*B0_deriv_deriv/2 + Pe*ux1/2
q[np.argmin(abs(k-2)), :] = kappa*xi*ux0*B0_deriv/4 - uy1*B0_deriv/4


for x in range(1, int(len(xi)-1)):
	for i in range(len(k)):
		rho_kp_2 = k[i]*1j*omega + kappa*kappa
		if k[i] == 0: #then it's cosine, with boundary = +- kappa
			b_k[i, 0]  = 0
			b_k[i, 1]  = -kappa*dx
			b_k[i, -1] = 0
			b_k[i, -2] = kappa*dx		

		b_k[i, x+1] = 2*b_k[i, x] - b_k[i, x-1] + np.real(dx*dx*(rho_kp_2*b_k[i, x] - q[i, x]))
		if i != len(k)-1:
			if i != 0:
				b_k[i, x+1]	+= dx*dx*np.real(- ux0[x]*b_k[i-1, x]/2 - np.conj(ux0[x])*b_k[i+1, x]/2)

for i in range(len(k)):
	plt.plot(xi, b_k[i,:])
#plt.axis([-1, 1, -5, 5])
plt.show()

plt.plot(xi, np.gradient(b_k[10, :], xi))
plt.show()

plt.plot(xi, np.real(q[8,:]+q[12,:]), label=r"$2\omega$")
plt.plot(xi, np.real(q[9,:]+q[11,:]), label=r"$1\omega$")
plt.plot(xi, np.real(q[10,:]), label=r"$0\omega$")
plt.legend(loc="best")
plt.ylabel(r"Source term $q$")
plt.xlabel(r"Vertical position $\xi$")
#plt.show()
