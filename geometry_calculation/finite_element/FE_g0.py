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


N = 500
xi = np.linspace(-1, 1, int(pow(10, 5)))


#system parameters
kappa = 0.5
Sc = 1.2
omega = 2.3
F0 = 3
Pe = 1


#implicitly defined parameters
gamma   = np.sqrt(1j*omega/Sc)
kappa_p = np.sqrt(gamma*gamma+kappa*kappa)
rho     = np.sqrt(1j*omega)
rho_p   = np.sqrt(rho*rho + kappa*kappa)
P_1     = (F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))

#analytic results up to first order
ux0  = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)
ux1 = (P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p *xi)/np.cosh(kappa_p)) + (F0*np.tanh(gamma)/gamma)*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma))
uy1 = (kappa*P_1*np.sinh(kappa)/(gamma*gamma))*(np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))

B0             =  Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma)))
B0_deriv       = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))
B0_deriv_deriv = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma))

f = -Pe*ux1 - kappa*kappa*xi*B0_deriv + 2*B0_deriv_deriv -Pe*kappa*ux0*np.cosh(kappa*xi)/np.sinh(kappa)


sol = finite_element_solver(N, xi, f, rho_p*rho_p, 0, 0)

plt.plot(xi[:-1], np.real(sol[:-1]-sol[0]), label="numerical solution")
plt.xlabel("x")



from numpy import *
rho_p   = sqrt(1j*omega+kappa*kappa)
rho     = sqrt(1j*omega)
gamma   = sqrt(1j*omega/Sc)

analytic_f1_solution  = np.zeros(len(xi), dtype="complex")
analytic_f1_solution += kappa*kappa*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)*sinh(rho)) * (xi*sinh(rho*xi)/(kappa*kappa) + 2*rho*cosh(rho*xi)/(kappa*kappa*kappa*kappa))
analytic_f1_solution += -kappa*kappa*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)*sinh(gamma))*(xi*sinh(gamma*xi)/(rho_p*rho_p-gamma*gamma) + 2*gamma*cosh(gamma*xi)/(rho_p*rho_p-gamma*gamma)**2)
analytic_f1_solution += F0*Pe*kappa*cosh(kappa*xi)/(gamma*gamma*rho*rho*sinh(kappa))
analytic_f1_solution +=  F0*Pe*kappa/(2*gamma*gamma*cosh(gamma)*sinh(kappa))*(cosh(xi*(gamma+kappa))/(gamma*gamma+2*kappa*gamma-rho*rho) + cosh(xi*(gamma-kappa))/(gamma*gamma-2*kappa*gamma-rho*rho) )
analytic_f1_solution += -2*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*kappa*kappa*(Sc-1))*(rho*cosh(rho*xi)/sinh(rho) - gamma*kappa*kappa*cosh(gamma*xi)/((rho_p*rho_p - gamma*gamma)*sinh(gamma)))
analytic_f1_solution += Pe*P_1*kappa*cosh(kappa)/(gamma*gamma)*(cosh(kappa*xi)/(cosh(kappa)*rho*rho) + cosh(kappa_p*xi)/(cosh(kappa_p)*(gamma*gamma-rho*rho)))
analytic_f1_solution += (Pe*F0*tanh(gamma)/gamma)*(xi*sinh(gamma*xi)*(gamma*gamma-rho_p*rho_p) - 2*gamma*cosh(gamma*xi))/(sinh(gamma)*(rho_p*rho_p-gamma*gamma)**2) - Pe*F0*tanh(gamma)*cosh(kappa_p*xi)/(gamma*(gamma*gamma-rho*rho)*cosh(kappa_p))

A = F0*Pe*kappa_p*sinh(kappa_p)*tanh(gamma)/(gamma*(gamma**2 - rho**2)*cosh(kappa_p)) - F0*Pe*(-2*gamma**2*sinh(gamma) + gamma*(gamma**2 - rho_p**2)*cosh(gamma) + (gamma**2 - rho_p**2)*sinh(gamma))*tanh(gamma)/(gamma*(gamma**2 - rho_p**2)**2*sinh(gamma)) - F0*Pe*kappa**2/(gamma**2*rho**2) + F0*Pe*kappa*((gamma - kappa)*sinh(gamma - kappa)/(-gamma**2 + 2*gamma*kappa + rho**2) - (gamma + kappa)*sinh(gamma + kappa)/(gamma**2 + 2*gamma*kappa - rho**2))/(2*gamma**2*sinh(kappa)*cosh(gamma)) - F0*Pe*kappa**2*(-2*gamma**2*sinh(gamma)/(gamma**2 - rho_p**2)**2 + gamma*cosh(gamma)/(gamma**2 - rho_p**2) + sinh(gamma)/(gamma**2 - rho_p**2))*tanh(gamma)/(gamma**3*(Sc - 1)*sinh(gamma)) - F0*Pe*kappa**2*(rho*cosh(rho)/kappa**2 + sinh(rho)/kappa**2 + 2*rho**2*sinh(rho)/kappa**4)*tanh(gamma)/(gamma**3*(Sc - 1)*sinh(rho)) + 2*F0*Pe*(gamma**2*kappa**2/(gamma**2 - rho_p**2) + rho**2)*tanh(gamma)/(gamma**3*kappa**2*(Sc - 1)) - P_1*Pe*kappa*(kappa*sinh(kappa)/(rho**2*cosh(kappa)) + kappa_p*sinh(kappa_p)/((gamma**2 - rho**2)*cosh(kappa_p)))*cosh(kappa)/gamma**2
#A = np.gradient(analytic_f1_solution, xi)[0]
my_sol_homo = cosh(rho_p*xi)*A/(rho_p*sinh(rho_p))
my_sol      = np.real(analytic_f1_solution+my_sol_homo-analytic_f1_solution[0]-my_sol_homo[0])

plt.plot(xi, np.real(my_sol-my_sol[0]), "--", label="analytic solution")
plt.legend(loc="best")
plt.show()