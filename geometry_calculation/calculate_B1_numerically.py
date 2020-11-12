import numpy as np 
import matplotlib.pyplot as plt 

k = np.arange(-5, 5, 1)
xi = np.linspace(-1, 1, int(1e3))
b_k = np.zeros((len(k), len(xi)), dtype="complex")
q = np.zeros((len(k), len(xi)), dtype="complex")
dx = xi[1]-xi[0]

kappa = 0.1
Sc = 1.2 
omega = 1.3
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

neg_freq_indx = np.argmin(abs(k+1))
pos_freq_indx = np.argmin(abs(k-1))

sol1 = np.zeros((len(k),len(xi)), dtype="complex")
sol2 = np.zeros((len(k),len(xi)), dtype="complex")

"""
plt.plot(xi, np.real(ux0))
plt.plot(xi, np.imag(ux0))
plt.show()

plt.plot(xi, np.real(q[pos_freq_indx, :]))
plt.plot(xi, np.imag(q[pos_freq_indx, :]))
plt.show()

plt.plot(xi, np.real(q[neg_freq_indx, :]))
plt.plot(xi, np.imag(q[neg_freq_indx, :]))
plt.show()
"""
p_np2 = 1j*omega*k + kappa*kappa
dxdx = dx*dx

for j in range(1, len(xi)-2):
	for i in range(len(k)):
		sol1[i, j+1] = 2*sol1[i, j] - sol1[i, j-1] + dxdx*p_np2[i]*sol1[i, j]
		
		if i == pos_freq_indx:
			sol1[i, j+1] -= dxdx*q[pos_freq_indx, j]

		if i == neg_freq_indx:
			sol1[i, j+1] -= dxdx*q[neg_freq_indx, j]

		if i == 0:
			sol1[i, j+1] -= dxdx*(sol1[i+1, j]*kappa*np.conj(ux0[j])/2)
		elif i == len(k)-1:
			sol1[i, j+1] -= dxdx*(sol1[i-1, j]*kappa*ux0[j+1]/2)
		else:
			sol1[i, j+1] -= dxdx*(sol1[i+1, j]*kappa*np.conj(ux0[j+1])/2 + sol1[i-1, j]*kappa*ux0[j+1]/2)
	sol1[:, -1] = sol1[:, -2]


neg_freq_indx = np.argmin(abs(k+2))
nul_freq_indx = np.argmin(abs(k-0))
pos_freq_indx = np.argmin(abs(k-2))

for j in range(1, len(xi)-2):
	for i in range(len(k)):
		sol2[i, j+1] = 2*sol2[i, j] - sol2[i, j-1] + dxdx*p_np2[i]*sol2[i, j]
		
		if i == neg_freq_indx:
			sol2[i, j+1] -= dxdx*q[neg_freq_indx, j]

		if i == nul_freq_indx:
			sol2[i, j+1] -= dxdx*q[nul_freq_indx, j]

		if i == pos_freq_indx:
			sol2[i, j+1] -= dxdx*q[pos_freq_indx, j]

		if i == 0:
			sol2[i, j+1] = dxdx*(sol2[i+1, j]*kappa*np.conj(ux0[j])/2)
		elif i == len(k)-1:
			sol2[i, j+1] = dxdx*(sol2[i-1, j]*kappa*ux0[j+1]/2)
		else:
			sol2[i, j+1] = dxdx*(sol2[i+1, j]*kappa*np.conj(ux0[j+1])/2 + sol1[i-1, j]*kappa*ux0[j+1]/2)
	sol1[:, -1] = sol1[:, -2]
		
"""
for x in range(1, int(len(xi)-1)):
	for i in range(len(k)):
		rho_kp_2 = k[i]*1j*omega + kappa*kappa
		if k[i] == 0: #then it's cosine, with boundary = +- kappa
			b_k[i, 0]  = 0
			b_k[i, 1]  = -kappa*dx
	
		b_k[i, x+1] = 2*b_k[i, x] - b_k[i, x-1] + dx*dx*(b_k[i,x] - q[i, x])
		
		if i == len(k)-1:
			b_k[i, x+1]	+= (-1)**i * dx*dx*(- ux0[x]*b_k[i-1, x]/2)
		elif i == 0:
			b_k[i, x+1]	+= (-1)**i * dx*dx*(- np.conj(ux0[x])*b_k[i+1, x]/2)
		else:
			b_k[i, x+1]	+=  (-1)**i * dx*dx*(- ux0[x]*b_k[i-1, x]/2 - np.conj(ux0[x])*b_k[i+1, x]/2)
"""
total_sol_1 = np.zeros((len(xi)), dtype="complex")
total_sol_2 = np.zeros((len(xi)), dtype="complex")

for i in range(len(k)):
	#plt.title("frequency = " + str(k[i]))
	plt.plot(xi, np.real(sol1[i,:]))
	plt.plot(xi, np.imag(sol1[i,:]), "--")
	if i % 2 == 0:
		total_sol_1 += np.imag(sol1[i,:])
		total_sol_2 += np.real(sol1[i,:])

	plt.show()
plt.show()
plt.title("real part of sin and cos")
plt.plot(xi, np.real(total_sol_2))
plt.plot(xi, np.real(total_sol_1))
plt.show()

plt.show()
plt.plot(xi, np.real(q[np.argmin(abs(k+2)),:]+q[np.argmin(abs(k-2)),:]), label=r"$2\omega$")
plt.plot(xi, np.real(q[np.argmin(abs(k+1)),:]+q[np.argmin(abs(k-1)),:]), label=r"$1\omega$")
plt.plot(xi, np.real(q[np.argmin(abs(k-0)),:]), label=r"$0\omega$")
plt.legend(loc="best")
plt.ylabel(r"Source term $q$")
plt.xlabel(r"Vertical position $\xi$")
#plt.show()
