import numpy as np 
import matplotlib.pyplot as plt 

def coupled_finite_element_solver(N, n, x, alpha, couple_forward, couple_backward, f, Bc0, Bc1):
	#define general parameters
	N_pos = np.linspace(min(x), max(x), N)
	Delta_x = N_pos[1]-N_pos[0]

	Nn = int(N*n)
	A   = np.zeros((Nn, Nn), dtype="complex")
	A_p = np.zeros((Nn, Nn), dtype="complex")
	b = np.zeros(Nn, dtype="complex")
	u = np.zeros((n, len(x)), dtype="complex")

	phi = np.zeros((N, len(x)), dtype="complex")

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

	for j in range(n):
		for i in range(N-1):
			A_p[j*N+i,     j*N+i] += 1
			A_p[j*N+i+1,   j*N+i] -= 1
			A_p[j*N+i,   j*N+i+1] -= 1
			A_p[j*N+i+1, j*N+i+1] += 1

			A[j*N+ i, j*N+ i]     += 2*alpha[j]
			A[j*N+ i+1, j*N+i]    += 1*alpha[j]
			A[j*N+ i, j*N+ i+1]   += 1*alpha[j]
			A[j*N+ i+1, j*N+i+1]  += 2*alpha[j]

	A_p *= 1/Delta_x
	A   *= Delta_x/6

	for i in range(n-1):
		A[(i+1)*N : (i+2)*N, i*N    :(i+1)*N] += couple_forward  * np.identity(N)
		A[i*N     : (i+1)*N, (i+1)*N:(i+2)*N] += couple_backward * np.identity(N)

	#calculate source vector
	for j in range(n):
		for i in range(N):
			b[j*N+i] = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[i]+Delta_x/2)))] + f[j, np.argmin(abs(x-(N_pos[i]-Delta_x/2)))])/2

		b[j*N+0]   = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[0] + Delta_x/2)))])/2
		b[j*N+N-1] = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[-1]- Delta_x/2)))])/2

		#if derivative is non-zero at boundary
		b[N*j+0]    -= Bc0[j]
		b[N*j+N-1]  += Bc1[j]

	sol = np.linalg.solve(A+A_p, b)

	#transfer back solution to regular basis
	for j in range(n):
		for i in range(N):
			u[j,:] += sol[j*N+i]*phi[i, :]
	return u

tol = 1e-6
k = np.array([-3, -2, -1, 0, 1, 2, 3])
xi = np.linspace(-1, 1, int(1e5))
Delta = xi[1]-xi[0]

#system parameters
kappa = 0.5
Sc = 1.2
omega = 2.3
F0 = 1
Pe = 3


#implicitly defined parameters
gamma = np.sqrt(1j*omega/Sc)
kappa_p = np.sqrt(gamma*gamma+kappa*kappa)
rho = np.sqrt(kappa*kappa+1j*omega)
P_1 = (F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))
p_np2 = 1j*omega*k + kappa*kappa

#analytic results up to first order
ux0  = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)
ux1 = (P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p *xi)/np.cosh(kappa_p)) + (F0*np.tanh(gamma)/gamma)*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma))
uy1 = (kappa*P_1*np.sinh(kappa)/(gamma*gamma))*(np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))

B0 = Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma)))
B0_deriv = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))
B0_deriv_deriv = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma))

#define source terms
q = np.zeros((len(k), len(xi)), dtype="complex")
q[np.argmin(abs(k+1)), :] = 0#Pe*np.conj(ux1)/2 + kappa*kappa*xi*np.conj(B0_deriv)/2- 2*np.conj(B0_deriv_deriv)/2
q[np.argmin(abs(k-1)), :] = 2*(Pe*ux1/2 + kappa*kappa*xi*B0_deriv/2 - 2*B0_deriv_deriv/2) - Pe*kappa*ux0*np.cosh(kappa*xi)/np.sinh(kappa)


#works for differential equation with constant terms, now just need coupeling to work as well
n = len(k) #number of vectors
N = 500 #length of each vector 
N_pos = np.linspace(min(xi), max(xi), N)

#coupleing between vectors
couple_forward  =  np.zeros(N, dtype="complex")
couple_backward =  np.zeros(N, dtype="complex")

for i in range(N):
	couple_backward[i] = 0#Delta*(ux0[np.argmin(abs(xi-(N_pos[i]+Delta/2)))] + ux0[np.argmin(abs(xi-(N_pos[i]-Delta/2)))])/2

couple_backward[0]   = 0#Delta*(ux0[np.argmin(abs(xi-(N_pos[0] + Delta/2)))])/2
couple_backward[-1]  = 0#Delta*(ux0[np.argmin(abs(xi-(N_pos[-1]- Delta/2)))])/2
couple_backward *= kappa 

couple_forward = np.conj(couple_backward)

#since integral over last and first basis function is half the value of the others
couple_backward[0]  *= 0.5
couple_backward[-1] *= 0.5 

couple_forward[0]  *= 0.5
couple_forward[-1] *= 0.5 

#boundary conditions
Bc0 = np.zeros(n, dtype="complex") #BC at xi = -1
Bc1 = np.zeros(n, dtype="complex") #BC at xi =  1

Bc0[np.argmin(abs(k-0))] =  -kappa
Bc1[np.argmin(abs(k-0))] =  kappa

sol = coupled_finite_element_solver(N, n, xi, p_np2, couple_backward, couple_forward, q, Bc0, Bc1)


for i in range(int(len(k)/2)+1):
	if k[i] < 1e-3:
		#plt.plot(xi, np.real(sol[i,:]+sol[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))
		
		if np.max(abs(np.imag(sol[i,:]+sol[-i-1,:]) - np.imag(sol[i-1,0]+sol[-i,0]))) > tol:
			print("LARGE IMAGINARY VALUE FOR " + str(k[i]) + "-OMEGA = ", np.max(abs(np.imag(sol[i-1,:]+sol[-i-1,:]) - np.imag(sol[i-1,0]+sol[-i-1,0]))))
	else:
		#plt.plot(xi, np.real(sol[0,:]), label=str(k[i]))

		if np.max(abs(np.imag(sol[0,:]-sol[0,0]))) > tol:
			print("LARGE IMAGINARY VALUE FOR 0-OMEGA = ", np.max(abs(np.imag(sol[0,:]))))

#plt.figure(1)
#plt.legend(loc="best")
#plt.show()

f0_g1 = sol/2


###

#reset source term for new solution
q = np.zeros((len(k), len(xi)), dtype="complex")
q[np.argmin(abs(k+2)), :] = kappa*xi*np.conj(ux0)*np.conj(B0_deriv)/4 - np.conj(uy1)*np.conj(B0_deriv)/4
q[np.argmin(abs(k-0)), :] = kappa*xi*(ux0*np.conj(B0_deriv) + np.conj(ux0)*B0_deriv) - uy1*np.conj(B0_deriv) - np.conj(uy1)*B0_deriv
q[np.argmin(abs(k-2)), :] = kappa*xi*ux0*B0_deriv/4 - uy1*B0_deriv/4

#boundary conditions
Bc0 = np.zeros(n, dtype="complex")
Bc1 = np.zeros(n, dtype="complex")

sol = coupled_finite_element_solver(N, n, xi, p_np2, couple_backward, couple_forward, q, Bc0, Bc1)

for i in range(int(len(k)/2)+1):
	if k[i] < 1e-3:
		#plt.plot(xi, -np.real(sol[i,0]+sol[-i-1,0]) + np.real(sol[i,:]+sol[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))
		
		if np.max(abs(np.imag(sol[i,:]+sol[-i-1,:]) - np.imag(sol[i,0]+sol[-i-1,0]))) > tol:
			print("LARGE IMAGINARY VALUE FOR " + str(k[i]) + "-OMEGA = ", np.max(abs(np.imag(sol[i,:]+sol[-i-1,:]) - np.imag(sol[i,0]+sol[-i-1,0]))))
	else:
		#plt.plot(xi, np.real(sol[0,:]), label=str(k[i]))

		if np.max(abs(np.imag(sol[0,:]-sol[0,0]))) > tol:
			print("LARGE IMAGINARY VALUE FOR 0-OMEGA = ", np.max(abs(np.imag(sol[0,:]))))

#plt.figure(1)
#plt.legend(loc="best")
#plt.show()

g0_f1 = sol/2
g0_f1 = np.zeros(np.shape(f0_g1))

full_sol = np.zeros((len(xi), 2)) #[xi, sin, cos]
for i in range(int(len(k)/2)+1):
	if abs(k[i]) % 2 == 0:
		print("even:", k[i], "and", k[-i-1])
		plt.figure(1) #cos figure
		plt.plot(xi, -np.real(f0_g1[i,0]+f0_g1[-i-1,0]) + np.real(f0_g1[i,:]+f0_g1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))
		full_sol[:, 1] += -np.real(f0_g1[i,0]+f0_g1[-i-1,0]) + np.real(f0_g1[i,:]+f0_g1[-i-1,:])

		plt.figure(2) #sin figure
		plt.plot(xi, -np.real(g0_f1[i,0]+g0_f1[-i-1,0]) + np.real(g0_f1[i,:]+g0_f1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))
		full_sol[:, 0] += -np.real(g0_f1[i,0]+g0_f1[-i-1,0]) + np.real(g0_f1[i,:]+g0_f1[-i-1,:])

	else:
		print("odd:", k[i], "and", k[-i-1])
		plt.figure(2) #sin figure
		plt.plot(xi, -np.real(f0_g1[i,0]+f0_g1[-i-1,0]) + np.real(f0_g1[i,:]+f0_g1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))
		full_sol[:, 0] += -np.real(f0_g1[i,0]+f0_g1[-i-1,0]) + np.real(f0_g1[i,:]+f0_g1[-i-1,:])

		plt.figure(1) #cos figure
		plt.plot(xi, -np.real(g0_f1[i,0]+g0_f1[-i-1,0]) + np.real(g0_f1[i,:]+g0_f1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))	
		full_sol[:, 1] += -np.real(g0_f1[i,0]+g0_f1[-i-1,0]) + np.real(g0_f1[i,:]+g0_f1[-i-1,:])

from numpy import *
rho_p   = sqrt(1j*omega+kappa*kappa)
rho     = sqrt(1j*omega)
gamma   = sqrt(1j*omega/Sc)

analytic_f1_solution  = np.zeros(len(xi), dtype="complex")
analytic_f1_solution +=  kappa*kappa*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)*sinh(rho)) * (xi*sinh(rho*xi)/(kappa*kappa) + 2*rho*cosh(rho*xi)/(kappa*kappa*kappa*kappa))
analytic_f1_solution += -kappa*kappa*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)*sinh(gamma))*(xi*sinh(gamma*xi)/(rho_p*rho_p-gamma*gamma) + 2*gamma*cosh(gamma*xi)/(rho_p*rho_p-gamma*gamma)**2)
analytic_f1_solution +=  F0*Pe*kappa*cosh(kappa*xi)/(gamma*gamma*rho*rho*sinh(kappa))
analytic_f1_solution +=  F0*Pe*kappa/(2*gamma*gamma*cosh(gamma)*sinh(kappa))*(cosh(xi*(gamma+kappa))/(gamma*gamma+2*kappa*gamma-rho*rho) + cosh(xi*(gamma-kappa))/(gamma*gamma-2*kappa*gamma-rho*rho) )
analytic_f1_solution += -2*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*kappa*kappa*(Sc-1))*(rho*cosh(rho*xi)/sinh(rho) - gamma*kappa*kappa*cosh(gamma*xi)/((rho_p*rho_p - gamma*gamma)*sinh(gamma)))
analytic_f1_solution +=  Pe*P_1*kappa*cosh(kappa)/(gamma*gamma)*(cosh(kappa*xi)/(cosh(kappa)*rho*rho) + cosh(kappa_p*xi)/(cosh(kappa_p)*(gamma*gamma-rho*rho)))
analytic_f1_solution += (Pe*F0*tanh(gamma)/gamma)*(xi*sinh(gamma*xi)*(gamma*gamma-rho_p*rho_p) - 2*gamma*cosh(gamma*xi))/(sinh(gamma)*(rho_p*rho_p-gamma*gamma)**2) - Pe*F0*tanh(gamma)*cosh(kappa_p*xi)/(gamma*(gamma*gamma-rho*rho)*cosh(kappa_p))

A = F0*Pe*kappa_p*sinh(kappa_p)*tanh(gamma)/(gamma*(gamma**2 - rho**2)*cosh(kappa_p)) - F0*Pe*(-2*gamma**2*sinh(gamma) + gamma*(gamma**2 - rho_p**2)*cosh(gamma) + (gamma**2 - rho_p**2)*sinh(gamma))*tanh(gamma)/(gamma*(gamma**2 - rho_p**2)**2*sinh(gamma)) - F0*Pe*kappa**2/(gamma**2*rho**2) + F0*Pe*kappa*((gamma - kappa)*sinh(gamma - kappa)/(-gamma**2 + 2*gamma*kappa + rho**2) - (gamma + kappa)*sinh(gamma + kappa)/(gamma**2 + 2*gamma*kappa - rho**2))/(2*gamma**2*sinh(kappa)*cosh(gamma)) - F0*Pe*kappa**2*(-2*gamma**2*sinh(gamma)/(gamma**2 - rho_p**2)**2 + gamma*cosh(gamma)/(gamma**2 - rho_p**2) + sinh(gamma)/(gamma**2 - rho_p**2))*tanh(gamma)/(gamma**3*(Sc - 1)*sinh(gamma)) - F0*Pe*kappa**2*(rho*cosh(rho)/kappa**2 + sinh(rho)/kappa**2 + 2*rho**2*sinh(rho)/kappa**4)*tanh(gamma)/(gamma**3*(Sc - 1)*sinh(rho)) + 2*F0*Pe*(gamma**2*kappa**2/(gamma**2 - rho_p**2) + rho**2)*tanh(gamma)/(gamma**3*kappa**2*(Sc - 1)) - P_1*Pe*kappa*(kappa*sinh(kappa)/(rho**2*cosh(kappa)) + kappa_p*sinh(kappa_p)/((gamma**2 - rho**2)*cosh(kappa_p)))*cosh(kappa)/gamma**2
#A = np.gradient(analytic_f1_solution, xi)[0]
my_sol_homo = cosh(rho_p*xi)*A/(rho_p*sinh(rho_p))
my_sol      = np.real(analytic_f1_solution+my_sol_homo-analytic_f1_solution[0]-my_sol_homo[0])


plt.figure(1)
plt.legend(loc="best")
plt.title(r"$\cos{\kappa\eta}$-solution")
plt.xlabel(r"x-axis $\xi$")
plt.ylabel(r"Brenner field $B(\xi)$")
plt.plot(xi, -np.cosh(-kappa)/np.sinh(kappa)+np.cosh(kappa*xi)/np.sinh(kappa), "--")

plt.savefig("figures/Brennerfield_cos.pdf")
plt.figure(2)
plt.legend(loc="best")
plt.title(r"$\sin{\kappa\eta}$-solution")
plt.xlabel(r"x-axis $\xi$")
plt.ylabel(r"Brenner field $B(\xi)$")
plt.plot(xi, my_sol, "--")
plt.plot(xi, my_sol*np.max(-np.real(f0_g1[np.argmin(abs(k-1)),0]+f0_g1[-np.argmin(abs(k-1))-1,0]) + np.real(f0_g1[np.argmin(abs(k-1)),:]+f0_g1[-np.argmin(abs(k-1))-1,:]))/np.max(my_sol), "k--")

plt.savefig("figures/Brennerfield_sin.pdf")

plt.figure(3)
plt.plot(xi, full_sol[:, 0])
plt.plot(xi, full_sol[:, 1])
plt.title(r"Full solution")
plt.xlabel(r"x-axis $\xi$")
plt.ylabel(r"Brenner field $B(\xi)$")
plt.savefig("figures/Brennerfield_full.pdf")
plt.show()