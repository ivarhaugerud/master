import numpy as np 
import matplotlib.pyplot as plt 

def finite_element_solver(N, x, f, a, Bc0, Bc1):
	#define general parameters
	phi = np.zeros((N, len(x)))
	N_pos = np.linspace(min(x), max(x), N)
	Delta_x = N_pos[1]-N_pos[0]

	A   = np.zeros((N, N), dtype="complex")
	A_p = np.zeros((N, N), dtype="complex")
	u   = np.zeros(len(x), dtype="complex")
	b   = np.zeros(N, dtype="complex")

	#find form of phi's given our interval
	for i in range(N):
		sta_index = np.argmin(abs(x-(N_pos[i]-Delta_x)))
		top_index = np.argmin(abs(x-(N_pos[i]        )))
		end_index = np.argmin(abs(x-(N_pos[i]+Delta_x)))
		phi[i, sta_index:top_index]   = np.linspace(0, 1, top_index-sta_index)
		phi[i, top_index:end_index]   = np.linspace(1, 0, end_index-top_index)
	
	phi[-1, -1] = 1 #some times the last element is not included, and is therefore set to zero if this is not done


	#calculate matrix elements using analytical results
	# A   = phi_i  * phi_j  matrix 
	# A_p = phi_i' * phi_j' matrix 

	for i in range(N-1):
		A_p[i, i]     += 1
		A_p[i+1, i]   -= 1
		A_p[i, i+1]   -= 1
		A_p[i+1, i+1] += 1

		A[i, i]     += 2
		A[i+1, i]   += 1
		A[i, i+1]   += 1
		A[i+1, i+1] += 2

	A_p *= 1/Delta_x
	A   *= a*Delta_x/6 # a=k^2

	#calculate source vector
	for i in range(len(b)):
		b[i] = -Delta_x*(f[np.argmin(abs(x-(N_pos[i]+Delta_x/2)))] + f[np.argmin(abs(x-(N_pos[i]-Delta_x/2)))])/2

	b[0]  = -Delta_x*(f[np.argmin(abs(x-(N_pos[0] + Delta_x/2)))])/2
	b[-1] = -Delta_x*(f[np.argmin(abs(x-(N_pos[-1]- Delta_x/2)))])/2

	#if derivative is non-zero at boundary
	b[0]   -= Bc0
	b[-1]  +=  Bc1

	sol = np.linalg.solve(A+A_p, b)

	#transfer back solution to regular basis
	for i in range(N):
		u += sol[i]*phi[i, :]

	return u

xi = np.linspace(-1, 1, int(pow(10, 5)))

F0 = 3
tau = 5
omega = 2*np.pi/tau
rho = np.sqrt(1j*omega)
Pe = 5

N = 500
kappa = 0.6
f = Pe*F0*(rho*np.cosh(rho*xi)/np.sinh(rho)-1)/(rho*rho) -F0*(1-xi*xi)/4
sol = finite_element_solver(N, xi, f, 1j*omega, 0, 0)
plt.plot(xi, sol, label="numerical solution")

#kappa = 0.25
#f = F0*(1-xi*xi)/4 - Pe*F0*(rho*np.cosh(rho*xi)/np.sinh(rho)-1)
#sol = finite_element_solver(N, xi, f, 1j*omega+kappa*kappa, 0, 0)
#plt.plot(xi, sol, "--", label="numerical solution")
ana = Pe*F0/(rho*rho)*(xi*np.sinh(rho*xi)/(2*np.sinh(rho))+1/(rho*rho)-np.cosh(rho*xi)/(2*rho*np.sinh(rho))*(1+rho/np.tanh(rho))) + F0*(1-xi*xi-2/(rho*rho)+2*np.cosh(rho*xi)/(rho*np.sinh(rho)))/(4*rho*rho) #+ 
plt.plot(xi, ana, "--")
plt.xlabel("x")
plt.legend(loc="best")
plt.show()



from numpy import *
b0_source = Pe*F0*kappa*(1-xi*xi)*ana/4
sol2 = Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(12*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) - 9*xi*sinh(rho*xi)/(rho*rho*sinh(rho))+3*xi*xi*cosh(rho*xi)/(rho*sinh(rho))-cosh(rho*xi)/(rho*sinh(rho))-xi*xi*xi*xi/12-xi*xi*xi*sinh(rho*xi)/(2*sinh(rho))+xi*xi/2+xi*sinh(rho*xi)/(2*sinh(rho)) + (1+rho/tanh(rho))/(2*rho*sinh(rho))*(6*cosh(rho*xi)/(rho*rho) - 4*xi*sinh(rho*xi)/rho + xi*xi*cosh(rho*xi)-cosh(rho*xi)))
sol1 = Pe*F0*F0*kappa/(16*rho*rho)*((xi*xi/2)*(1-2/(rho*rho))+xi*xi*xi*xi*(1/(rho*rho)-1)/6 + xi*xi*xi*xi*xi*xi/30 + 2*cosh(rho*xi)/(rho*rho*rho*sinh(rho))*(1-xi*xi-6/(rho*rho)) + 8*xi*sinh(rho*xi)/(rho*rho*rho*rho*sinh(rho)))
A = -kappa*(5*F0**2*Pe**2*rho**2/12 - F0**2*Pe**2*rho**2/tanh(rho)**2/4 - 3/4*F0**2*Pe**2*rho/tanh(rho) + F0**2*Pe**2 + F0**2*Pe*rho**4/30 - F0**2*Pe*rho**2/12 + F0**2*Pe*rho/tanh(rho)/4 - F0**2*Pe/4 + rho**6)/rho**6
sol2 += A*cosh(kappa*xi)/(kappa*sinh(kappa))
ana = sol1+sol2

sol = finite_element_solver(N, xi, b0_source, kappa*kappa, kappa, -kappa)
plt.plot(xi, sol-sol[0], label="numerical solution")
plt.plot(xi, ana-ana[0], "--")
plt.show()
"""
plt.plot(x[:-1], abs(ana[:-1]-sol[:-1])/abs(ana[:-1]))
plt.show()


N = 100
x = np.linspace(0, 1, int(pow(10, 5)))
f = (1+np.pi*np.pi)*np.cos(np.pi*x)+1.5
sol = finite_eleme"nt_solver(N, x, f, 1, 0, 0)

ana = -np.cos(np.pi*x)-1.5
plt.plot(x, ana, label="analytical solution")
plt.plot(x, sol, "--", label="numerical solution")
plt.xlabel("x")
plt.legend(loc="best")
plt.savefig("figures/test_FE_2.png")
plt.show()

plt.plot(x[:-1], abs(ana[:-1]-sol[:-1])/abs(ana[:-1]))
plt.show()
"""

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