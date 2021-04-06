import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
import scipy.integrate as scpi 

def finite_element_solver(N, x, alpha, f):
	#define general parameters
	phi = np.zeros((N, len(x)))
	N_pos = np.linspace(min(x), max(x), N)
	Delta_x = N_pos[1]-N_pos[0]

	A   = np.zeros((N, N),  dtype="complex")
	A_p = np.zeros((N, N),  dtype="complex")
	b = np.zeros(N,          dtype="complex")
	u = np.zeros(len(x), dtype="complex")

	phi = np.zeros((N, len(x)))

	#find form of phi's given our interval

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
		A_p[i,     i] += 1
		A_p[i+1,   i] -= 1
		A_p[i,   i+1] -= 1
		A_p[i+1, i+1] += 1

		A[i, i]      += 2
		A[i+1, i]    += 1
		A[i, i+1]    += 1
		A[i+1, i+1]  += 2

	A_p *= 1/Delta_x
	A   *= alpha*Delta_x/6

	#calculate source vector
	for i in range(N):
		b[i] = -Delta_x*(2*f[np.argmin(abs(x-(N_pos[i]+Delta_x/2)))] + 2*f[np.argmin(abs(x-(N_pos[i])))] + 2*f[np.argmin(abs(x-(N_pos[i]-Delta_x/2)))])/6

	b[0]   = -Delta_x*(f[np.argmin(abs(x-(N_pos[0 ])))] + 2*f[np.argmin(abs(x-(N_pos[0] + Delta_x/2)))])/6
	b[-1]  = -Delta_x*(f[np.argmin(abs(x-(N_pos[-1])))] + 2*f[np.argmin(abs(x-(N_pos[-1]- Delta_x/2)))])/6

	sol = np.linalg.solve(A+A_p, b)
	#transfer back solution to regular basis
	for i in range(N):
		u[:] += sol[i]*phi[i, :]
	return u, sol

def coupled_finite_element_solver(N, n, x, alpha, couple_forward, couple_backward, f, Bc0, Bc1, k):
	#define general parameters
	phi = np.zeros((N, len(x)))
	N_pos = np.linspace(min(x), max(x), N)
	Delta_x = N_pos[1]-N_pos[0]

	Nn = int(N*n)
	A   = np.zeros((Nn, Nn),  dtype="complex")
	A_p = np.zeros((Nn, Nn),  dtype="complex")
	b = np.zeros(Nn,          dtype="complex")
	u = np.zeros((n, len(x)), dtype="complex")

	phi = np.zeros((N, len(x)))

	#find form of phi's given our interval

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
		A[(i+1)*N : (i+2)*N, i*N    :(i+1)*N] = couple_forward  * np.power(-1, abs(k[i]+1))     #was originally not a minus sign there
		A[i*N     : (i+1)*N, (i+1)*N:(i+2)*N] = couple_backward * np.power(-1, abs(k[i]))       #was originally not a minus sign there

	#print(np.real(A))
	#print(np.real(np.dot(A, np.ones(N*n))))
	#calculate source vector
	for j in range(n):
		for i in range(N):
			b[j*N+i] = -Delta_x*(2*f[j, np.argmin(abs(x-(N_pos[i]+Delta_x/2)))] + 2*f[j, np.argmin(abs(x-(N_pos[i])))] + 2*f[j, np.argmin(abs(x-(N_pos[i]-Delta_x/2)))])/6

		b[j*N+0]   = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[0 ])))] + 2*f[j, np.argmin(abs(x-(N_pos[0] + Delta_x/2)))])/6
		b[j*N+N-1] = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[-1])))] + 2*f[j, np.argmin(abs(x-(N_pos[-1]- Delta_x/2)))])/6

		#if derivative is non-zero at boundary
		b[N*j+0]    -= Bc0[j]
		b[N*j+N-1]  += Bc1[j]

	sol = np.linalg.solve(A+A_p, b)
	#transfer back solution to regular basis
	for j in range(n):
		for i in range(N):
			u[j,:] += sol[j*N+i]*phi[i, :]
	return u, sol

tol   = 1e-6
k     = np.arange(-10, 10+0.01, 1)
xi    = np.linspace(-1, 1, int(1e5))
Lx = 4.48
kappa = 2*np.pi/Lx
nus = np.logspace(-2, 2, 10)[:9]
#kappas   = np.arange(0.2, 2.201, 0.1) #np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])

#system parameters
#nu  = 1.2
#D   = 1.0
#F0  = 12/nu
#Pe = 1/D

nu = 1.2
D = 0.1
F0 = 12/nu 
Pe = 1/D
tau = np.arange(2.0, 4.0, 0.25)
tau = np.array([3.0])
kappas = np.arange(0.5, 1.5, 0.05)
kappas   = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2]) #0.2

#Lx = np.array([125.66, 62.83, 41.89])
#kappas = 2*np.pi/Lx
#kappas    = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.5]) #0.2
omegas = 2*np.pi/tau #np.logspace(-1.5, 2.5, 10)
D_parallels = np.zeros((len(kappas), len(tau)))

"""
D_k_o = np.load("data/vary_kappa_omega.npy")
print(np.shape(D_k_o))
contin_kappa = np.linspace(min(kappas), max(kappas), int(1e4))

D0 = np.zeros(len(omegas))

from numpy import *
for i in range(len(omegas)):
	omega = omegas[i]
	gamma   = np.sqrt(1j*omega/nu)
	rho     = np.sqrt(1j*omega/D)

	D0[i]    = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))
	#plt.plot(kappas, D_k_o[:, i], label=r"$\omega=%3.2f$" % omegas[i])
	res = np.sqrt(omegas[i]/(2*D))
	contin_D = interp1d(kappas, D_k_o[:, i], kind="cubic")(contin_kappa)
	plt.plot(res, D0[i]+contin_D[np.argmin(abs(contin_kappa-res))], "ko")
	plt.plot(contin_kappa, D0[i]+contin_D)
	plt.legend(loc="best")
plt.show()
"""

for K in range(len(tau)):
	for O in range(len(kappas)):
		#kappa = kappas[K]
		#omega = omegas[K]
		omega = omegas[K]
		kappa = kappas[O]

		#implicitly defined parameters
		gamma   = np.sqrt(1j*omega/nu)
		rho     = np.sqrt(1j*omega/D)
		kappa_p = np.sqrt(gamma*gamma + kappa*kappa)
		P_1     = ((F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p))) )
		p_np2   = rho*rho*k + kappa*kappa

		#analytic results
		ux0 = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
		#ux1 = 0.5*((  (P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p*xi)/np.cosh(kappa_p))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma))  ))
		#uy1 = 0.5*(kappa*P_1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))
		#ux2 = 0.5*(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))
		ux2 = 0.5*(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))
		uy1 = 0.5*(kappa*P_1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))
		ux1 = 0.5*((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p*xi)/np.cosh(kappa_p))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))
		ux2 -= scpi.trapz(ux2, xi)/2 + scpi.trapz(ux1, xi)/4

		B0             = (Pe*F0*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma)))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma))) + Pe*F0*np.tanh(gamma)/(2*gamma*gamma*gamma*rho*rho)
		B0_deriv       = (Pe*F0*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))
		B0_deriv_deriv = (Pe*F0*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma))

		#define source terms
		q = np.zeros((len(k), len(xi)), dtype="complex")
		q[np.argmin(abs(k+2)), :] = Pe*kappa*xi*np.conj(ux0)*np.conj(B0_deriv) - Pe*np.conj(uy1)*np.conj(B0_deriv)
		q[np.argmin(abs(k+1)), :] = Pe*np.conj(ux1) + kappa*kappa*xi*np.conj(B0_deriv) - 2*np.conj(B0_deriv_deriv)
		q[np.argmin(abs(k-0)), :] = Pe*kappa*xi*(ux0*np.conj(B0_deriv) + np.conj(ux0)*B0_deriv)- Pe*uy1*np.conj(B0_deriv) - Pe*np.conj(uy1)*B0_deriv
		q[np.argmin(abs(k-1)), :] = Pe*ux1          + kappa*kappa*xi*B0_deriv          - 2*B0_deriv_deriv 
		q[np.argmin(abs(k-2)), :] = Pe*kappa*xi*ux0*B0_deriv                   - Pe*uy1*B0_deriv

		#works for differential equation with constant terms, now just need coupeling to work as well
		n = len(k) #number of vectors
		N = 350
		N_pos = np.linspace(-1, 1, N)
		Delta = N_pos[1]-N_pos[0]

		#coupleing between vectors
		couple_backward =  np.zeros((N, N), dtype="complex")

		for i in range(N-1):
			couple_backward[i, i]   = Delta*(ux0[np.argmin(abs(xi-(N_pos[i]-Delta/2)))] + 2*ux0[np.argmin(abs(xi-(N_pos[i])))] + ux0[np.argmin(abs(xi-(N_pos[i]+Delta/2)))])/6
			couple_backward[i, i+1] = Delta*ux0[np.argmin(abs(xi-(N_pos[i]+Delta/2)))]/6
			couple_backward[i+1, i] = Delta*ux0[np.argmin(abs(xi-(N_pos[i]+Delta/2)))]/6

		couple_backward[0, 0]    = Delta*(ux0[np.argmin(abs(xi-(N_pos[0])))]  + ux0[np.argmin(abs(xi-(N_pos[0]  + Delta/2)))])/6
		couple_backward[-1, -1]  = Delta*(ux0[np.argmin(abs(xi-(N_pos[-1])))] + ux0[np.argmin(abs(xi-(N_pos[-1] - Delta/2)))])/6
		couple_backward         *= Pe*kappa # unsure if factor of 1/2 is correct here
		couple_forward = np.conj(couple_backward)

		#boundary conditions
		Bc0 = np.zeros(n, dtype="complex") #BC at xi = -1
		Bc1 = np.zeros(n, dtype="complex") #BC at xi =  1

		Bc0[np.argmin(abs(k-0))] =  kappa
		Bc1[np.argmin(abs(k-0))] = -kappa

		f0_g1, coeff_f0g1 = coupled_finite_element_solver(N, n, xi, p_np2, couple_backward, couple_forward, -q, Bc0, Bc1, k)

		for i in range(int(len(k))):
			if np.max(abs(np.imag(f0_g1[i,:]+f0_g1[-i-1,:]) - np.imag(f0_g1[i-1,0]+f0_g1[-i,0])))/2 > tol:
				print("LARGE IMAGINARY VALUE FOR " + str(k[i]) + "-OMEGA = ", np.max(abs(np.imag(f0_g1[i-1,:]+f0_g1[-i-1,:]) - np.imag(f0_g1[i-1,0]+f0_g1[-i-1,0]))))
		
		g0_f1 = np.zeros(np.shape(f0_g1))
		coeff_g0f1 = np.zeros(np.shape(coeff_f0g1))

		B_plus   = np.zeros((len(xi),  len(k)), dtype="complex")  #sin-solution
		B_minus  = np.zeros((len(xi),  len(k)), dtype="complex")  #cos-solution
		B_plus_coeff    = np.zeros((N, len(k)), dtype="complex")  #sin-solution
		B_minus_coeff   = np.zeros((N, len(k)), dtype="complex")  #cos-solution

		for i in range(len(k)):
			if abs(k[i]) % 2 == 0:
				B_minus[:, i] = f0_g1[i,:]
				B_plus[:, i]  = g0_f1[i,:]
				B_minus_coeff[:, i] = coeff_f0g1[i*N:(i+1)*N]
				B_plus_coeff[:, i]  = coeff_g0f1[i*N:(i+1)*N]


			else:
				B_plus[:, i]  = f0_g1[i,:]
				B_minus[:, i] = g0_f1[i,:]
				B_minus_coeff[:, i] = coeff_g0f1[i*N:(i+1)*N]
				B_plus_coeff[:, i]  = coeff_f0g1[i*N:(i+1)*N]

		#define source terms

		q = np.zeros((len(k), len(xi)), dtype="complex")
		q[np.argmin(abs(k-1)), :] = 0.5*kappa*kappa*xi*B0_deriv          + 0.5*(3+kappa*kappa*xi*xi)*B0_deriv_deriv          + Pe*ux2
		q[np.argmin(abs(k+1)), :] = 0.5*kappa*kappa*xi*np.conj(B0_deriv) + 0.5*(3+kappa*kappa*xi*xi)*np.conj(B0_deriv_deriv) + Pe*np.conj(ux2)

		for i in range(n):
			B_min_deriv        = interp1d(N_pos, np.gradient(B_minus_coeff[:,i], N_pos), kind='cubic')(xi)
			B_plus_deriv       = interp1d(N_pos, np.gradient(B_plus_coeff[:,i], N_pos), kind='cubic')(xi)
			B_plus_deriv_deriv = interp1d(N_pos, np.gradient(np.gradient(B_plus_coeff[:,i], N_pos), N_pos), kind='cubic')(xi)

			q[i, :] += -0.5*kappa*kappa*xi*B_plus_deriv - B_plus_deriv_deriv
			if i != int(len(k)-1):
				q[i+1, :] += Pe*0.5*(ux0*kappa*xi*B_min_deriv + ux1*kappa*B_minus[:, i] - uy1*B_min_deriv)

			if i != 0:
				q[i-1, :] += Pe*0.5*(np.conj(ux0)*kappa*xi*B_min_deriv + np.conj(ux1)*kappa*B_minus[:, i] - np.conj(uy1)*B_min_deriv)

		sol       = np.zeros((len(xi), n), dtype="complex")
		sol_coeff = np.zeros((N, n), dtype="complex")

		for i in range(n):
			if abs(k[i]) > 1e-4:
				#plt.plot(xi, -q[i, :])
				sol[:,i], sol_coeff[:, i] = finite_element_solver(N, xi, rho*rho*k[i], -q[i,:])
		#plt.show()
		for i in range(n):
			if np.max(abs(np.imag(sol[:, i]+sol[:,-i-1]))) > tol:
				print("LARGE IMAGINARY VALUE FOR " + str(k[i]) + "-OMEGA = ", np.max(abs(np.imag(sol[:, i]+sol[:, -i-1]))))

		
		t        = np.linspace(0, 2*np.pi/omega, int(1e4))
		new_k    = np.arange(-2*max(k), 2*max(k)+1e-3, 1)
		D_eff_xi = np.zeros((len(xi), len(new_k)), dtype="complex")
		total_D  = np.zeros(len(t),                dtype="complex")
		D_eff    = np.zeros(len(new_k),            dtype="complex")

		b0                  = np.zeros((len(xi), n), dtype="complex")
		b0_deriv            = np.zeros((len(xi), n), dtype="complex")
		b0_deriv_deriv      = np.zeros((len(xi), n), dtype="complex")
		B1_min_deriv        = np.zeros((len(xi), n), dtype="complex")
		B1_plus_deriv       = np.zeros((len(xi), n), dtype="complex")
		B1_plus_deriv_deriv = np.zeros((len(xi), n), dtype="complex")
		B2_0_deriv          = np.zeros((len(xi), n), dtype="complex")

		b0[:, np.argmin(abs(k-1))]    = B0
		b0[:, np.argmin(abs(k+1))]    = np.conj(B0)

		b0_deriv[:, np.argmin(abs(k-1))] = B0_deriv
		b0_deriv[:, np.argmin(abs(k+1))] = np.conj(B0_deriv)

		b0_deriv_deriv[:, np.argmin(abs(k-1))] = B0_deriv_deriv 
		b0_deriv_deriv[:, np.argmin(abs(k+1))] = np.conj(B0_deriv_deriv)

		for i in range(n):
			B1_min_deriv[:, i]        = interp1d(N_pos, np.gradient(B_minus_coeff[:,i], N_pos), kind='cubic')(xi)
			B1_plus_deriv[:, i]       = interp1d(N_pos, np.gradient(B_plus_coeff[:,i],  N_pos), kind='cubic')(xi)
			B1_plus_deriv_deriv[:, i] = interp1d(N_pos, np.gradient(np.gradient(B_plus_coeff[:,i], N_pos), N_pos), kind='cubic')(xi)
			B2_0_deriv[:, i]          = interp1d(N_pos, np.gradient(sol_coeff[:,i], N_pos),     kind='cubic')(xi)

		"""
		for i in range(len(k)):
			D_eff_xi[:, np.argmin(abs(new_k -k[i]))] += (kappa*xi*B1_min_deriv[:, i] + kappa*B_minus[:, i])
			for j in range(len(k)):
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_plus[:,j]*B_plus[:,i]    + B1_plus_deriv[:, i]*B1_plus_deriv[:, j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_minus[:,j]*B_minus[:,i]  + B1_min_deriv[:, i] *B1_min_deriv[:, j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += b0_deriv[:, i]*(2*B2_0_deriv[:,j] - B1_plus_deriv[:,j] - kappa*kappa*xi*B_plus[:, j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5*(1+kappa*kappa*xi*xi)*b0_deriv[:, i]*b0_deriv[:, j]
		#	for i in range(len(k)):
	#		D_eff_xi[:, np.argmin(abs(new_k -k[i]))] += 0.5*(kappa*xi*B1_min_deriv[:, i] + kappa*B_minus[:, i] + kappa*xi*np.conjugate(B1_min_deriv[:, i]) + kappa*np.conjugate(B_minus[:, i])) 
		
			for j in range(len(k)):
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5*b0_deriv[:,i]*np.conjugate(b0_deriv[:,j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_plus[:,j] *np.conjugate(B_plus[:,i])   + B1_plus_deriv[:, i]*np.conjugate(B1_plus_deriv[:, j]))
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_minus[:,j]*np.conjugate(B_minus[:,i])  + B1_min_deriv[:, i] *np.conjugate(B1_min_deriv[:, j]))
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += b0_deriv[:, i]*np.conjugate((B2_0_deriv[:,j] - 0.5*B1_plus_deriv[:,j]))
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += np.conjugate(b0_deriv[:, i])*(B2_0_deriv[:,j] - 0.5*B1_plus_deriv[:,j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5*(kappa*kappa*xi*xi)*b0_deriv[:, i]*np.conjugate(b0_deriv[:, j]) 
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += -0.5*kappa*kappa*xi*xi*(b0_deriv[:,i]*np.conjugate(B1_plus_deriv[:, j]) + np.conjugate(b0_deriv[:,i])*B1_plus_deriv[:,j]) 
		"""
		
		for i in range(len(k)):
			D_eff_xi[:, np.argmin(abs(new_k -k[i]))] += (kappa*xi*B1_min_deriv[:, i] + kappa*B_minus[:, i])
			for j in range(len(k)):
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_plus[:,j]*B_plus[:,i]    + B1_plus_deriv[:, i]*B1_plus_deriv[:, j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_minus[:,j]*B_minus[:,i]  + B1_min_deriv[:, i] *B1_min_deriv[:, j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += b0_deriv[:, i]*(2*B2_0_deriv[:,j]-kappa*kappa*xi*B_plus[:,j] - B1_plus_deriv[:,j])
				D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5*(1+kappa*kappa*xi*xi)*b0_deriv[:, i]*b0_deriv[:, j]
		
		for i in range(len(new_k)):
			D_eff[i] = scpi.trapz(D_eff_xi[:, i], xi)/2
			total_D += np.exp(1j*omega*new_k[i]*t)*D_eff[i]

		#print("HERE:", np.max(np.imag(total_D)))

		D_parallels[O, K] = scpi.trapz(np.real(total_D), t)/(2*np.pi/(omega))
		print(omegas[K], kappas[O], D_parallels[O, K])
		np.save("data/total_D_kappa"+str(kappa)[:4], D_eff)

		plt.figure(1)
		plt.plot(new_k, np.real(D_eff))
		plt.xlabel(r"frequency [$\omega$]", fontsize=12)
		plt.ylabel(r"Amplitude", fontsize=12)
		plt.savefig("figures/fourier_D.pdf")

		plt.figure(2)
		plt.plot(t, np.real(total_D))
		plt.xlabel(r"time [$t$]", fontsize=12)
		plt.ylabel(r"Brenner field", fontsize=12)
		plt.savefig("figures/Brenner_field_vs_t.pdf")
		#plt.show()
np.save("data/vary_kappa_D0.1", D_parallels)
plt.show()