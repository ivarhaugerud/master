import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scp 

x = np.linspace(0, 1, int(pow(2, 10)))
N = 100
phi = np.zeros((N, len(x)))
phi_p = np.zeros((N, len(x)))
N_pos = np.linspace(min(x), max(x), N)
Delta_x = N_pos[1]-N_pos[0]

for i in range(N):
	sta_index = np.argmin(abs(x-(N_pos[i]-Delta_x)))
	top_index = np.argmin(abs(x-(N_pos[i]        )))
	end_index = np.argmin(abs(x-(N_pos[i]+Delta_x)))

	phi[i, sta_index:top_index] = np.linspace(0, 1, top_index-sta_index)
	phi[i, top_index:end_index+1] = np.linspace(1, 0, end_index+1-top_index)

	phi_p[i, :] =  np.gradient(phi[i, :], x)


A   = np.zeros((N, N))
A_p = np.zeros((N, N))

f = (1+np.pi*np.pi)*np.sin(np.pi*x)
b = np.zeros(N)

BC_cond_0 = 0
BC_cond_1 = 0


psi = BC_cond_0*phi[0, :] + BC_cond_1*phi[-1, :]

for i in range(N):
	for j in range(N):
		if abs(i-j)-0.1 < 1:
			A[i,j]    = 0#scp.trapz(phi[i,:] *    phi[j,:], x)
			A_p[i, j] = scp.trapz(phi_p[i,:] * phi_p[j,:], x)
	b[i] = scp.trapz(f*phi[i,:], x)


sol = np.linalg.solve(A+A_p, b)
u = np.zeros(len(x))
for i in range(N):
	u += sol[i]*phi[i, :]

plt.plot(x, u)
plt.plot(x, np.sin(np.pi*x))
plt.show()
"""
k = np.arange(-3, 3, 1)
xi = np.linspace(-1, 1, int(1e5))
q = np.zeros((len(k), len(xi)), dtype="complex")
dx = xi[1]-xi[0]

kappa = 0.5
Sc = 1.2 
omega = 1.35
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

p_np2 = 1j*omega*k + kappa*kappa
dxdx = dx*dx

for j in range(1, len(xi)-1):
	for i in range(len(k)):
		sol1[i, j+1] = 2*sol1[i, j] - sol1[i, j-1] -  dxdx*(p_np2[i]*sol1[i, j])
		
		if i == pos_freq_indx:
			sol1[pos_freq_indx, j+1] += -dxdx*xi[j]*xi[j]#dxdx*q[pos_freq_indx, j]

		if i == neg_freq_indx:
			sol1[neg_freq_indx, j+1] += -dxdx*xi[j]*xi[j]#dxdx*q[neg_freq_indx, j]
		
		#if i == 0:
		#	sol1[i, j+1] -= dxdx*(sol1[i+1, j]*kappa*np.conj(ux0[j])/2)
		#elif i == len(k)-1:
		#	sol1[i, j+1] -= dxdx*(sol1[i-1, j]*kappa*ux0[j]/2)
		#else:
		#	sol1[i, j+1] -= dxdx*(sol1[i+1, j]*kappa*np.conj(ux0[j])/2 + sol1[i-1, j]*kappa*ux0[j]/2)

"""