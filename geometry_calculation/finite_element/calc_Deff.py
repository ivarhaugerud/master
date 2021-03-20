import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scpi
from scipy.interpolate import interp1d

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
		A[(i+1)*N : (i+2)*N, i*N    :(i+1)*N] += couple_forward  * np.power(-1, abs(k[i]+1))
		A[i*N     : (i+1)*N, (i+1)*N:(i+2)*N] += couple_backward * np.power(-1, abs(k[i]))

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
k     = np.array([-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
xi    = np.linspace(-1, 1, int(1e5))

#system parameters
omega = (2*np.pi)/3.6
D = 0.1
nu = 1.2
F0 = 12
Pe = 1/D
kappa = 0.5

gamma   = np.sqrt(1j*omega/nu)
rho     = np.sqrt(1j*omega/D)
kappa_p = np.sqrt(gamma*gamma + kappa*kappa)
P_1     = ((F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p))) )
p_np2   = rho*rho*k + kappa*kappa

#analytic results
ux0 = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
ux1 = ((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p *xi)/np.cosh(kappa_p)) + (F0*np.tanh(gamma)/gamma)*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))/2
uy1 = ((kappa*P_1*np.sinh(kappa)/(gamma*gamma))*(np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa)))/2
ux2 = P_1*np.sinh(kappa)*(kappa*kappa*xi*np.sinh(kappa*xi)/np.sinh(kappa) - kappa_p*kappa_p*xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p) + gamma*gamma*np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma))
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
N = 200
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
couple_backward         *= Pe*kappa
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

"""
plt.plot(N_pos, coeff_f0g1[N*Z:N*(1+Z)])
plt.show()

plt.plot(N_pos, np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos))
plt.show()
"""
for i in range(int(len(k)/2)+1):
	plt.plot(xi, f0_g1[i, :]+f0_g1[-i-1, :], label=(str(k[-i-1])))
plt.legend(loc="best")
plt.show()

Z = np.argmin(abs(k-0))
plt.figure(0)
plt.title(str(0))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] + Pe*kappa*ux0*f0_g1[Z-1, :] + Pe*kappa*np.conj(ux0)*f0_g1[Z+1, :])
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)


Z = np.argmin(abs(k-1))
plt.figure(1)
plt.title(str(1))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] - Pe*kappa*ux0*f0_g1[Z-1, :] - Pe*kappa*np.conj(ux0)*f0_g1[Z+1, :])
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

Z = np.argmin(abs(k-2))
plt.figure(2)
plt.title(str(2))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] + Pe*kappa*ux0*f0_g1[Z-1, :] + Pe*kappa*np.conj(ux0)*f0_g1[Z+1, :])
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

Z = np.argmin(abs(k-3))
plt.figure(3)
plt.title(str(3))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] - Pe*kappa*ux0*f0_g1[Z-1, :] - Pe*kappa*np.conj(ux0)*f0_g1[Z+1, :])
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

Z = np.argmin(abs(k-4))
plt.figure(4)
plt.title(str(4))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] + Pe*kappa*ux0*f0_g1[Z-1, :] + Pe*kappa*np.conj(ux0)*f0_g1[Z+1, :])
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

Z = np.argmin(abs(k-5))
plt.figure(5)
plt.title(str(5))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] - Pe*kappa*ux0*f0_g1[Z-1, :] - Pe*kappa*np.conj(ux0)*f0_g1[Z+1, :])
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)
plt.show()

for i in range(int(len(k)/2)+1):
	if abs(k[i]) % 2 == 0:
		print("even:", k[i], "and", k[-i-1])
		plt.figure(1) #cos figure
		plt.plot(xi, np.real(f0_g1[i,:]+f0_g1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))

	else:
		print("odd:", k[i], "and", k[-i-1])
		plt.figure(2) #sin figure
		plt.plot(xi, np.real(f0_g1[i,:]+f0_g1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))


B_plus   = np.zeros((len(xi),  len(k)), dtype="complex")  #sin-solution
B_minus  = np.zeros((len(xi),  len(k)), dtype="complex")  #cos-solution
B_plus_coeff    = np.zeros((N, len(k)), dtype="complex") #sin-solution
B_minus_coeff   = np.zeros((N, len(k)), dtype="complex") #sin-solution

for i in range(len(k)):
	if abs(k[i]) % 2 == 0:
		B_minus[:, i] = f0_g1[i,:]
		B_minus_coeff[:, i] = coeff_f0g1[i*N:(i+1)*N]

	else:
		B_plus[:, i]  = f0_g1[i,:]
		B_plus_coeff[:,i]  = coeff_f0g1[i*N:(i+1)*N]

plt.figure(1)
plt.legend(loc="best")
plt.title(r"$\cos{\kappa\eta}$-solution")
plt.xlabel(r"x-axis $\xi$")
plt.ylabel(r"Brenner field $B(\xi)$")
#plt.plot(xi, np.cosh(kappa*xi)/np.sinh(kappa)-1/np.tanh(kappa), "--")

plt.savefig("figures/Brennerfield_cos.pdf")
plt.figure(2)
plt.legend(loc="best")
plt.title(r"$\sin{\kappa\eta}$-solution")
plt.xlabel(r"x-axis $\xi$")
plt.ylabel(r"Brenner field $B(\xi)$")
plt.savefig("figures/Brennerfield_sin.pdf")

np.save("data/B_plus", B_plus)
np.save("data/B_minus", B_minus)
np.save("data/B_minus_coeff", B_minus_coeff)
np.save("data/B_plus_coeff", B_plus_coeff)
np.save("data/k", k)
plt.show()

def finite_element_solver(N, x, f, double_deriv, a, Bc0, Bc1, laplace):
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

	#if laplace:
	#	A_p[0,-1]    += 1
	#	A_p[-1,0]    += 1

	A_p *= 1/Delta_x
	A   *= a*Delta_x/6 # a=k^2

	#calculate source vector
	for i in range(len(b)):
		b[i] = -Delta_x*(2*f[np.argmin(abs(x-(N_pos[i]+Delta_x/2)))] + 2*f[np.argmin(abs(x-(N_pos[i])))] + 2*f[np.argmin(abs(x-(N_pos[i]-Delta_x/2)))])/6

	b[0]  = -Delta_x*(f[np.argmin(abs(x-(N_pos[0 ])))] + 2*f[np.argmin(abs(x-(N_pos[0] + Delta_x/2)))])/6
	b[-1] = -Delta_x*(f[np.argmin(abs(x-(N_pos[-1])))] + 2*f[np.argmin(abs(x-(N_pos[-1]- Delta_x/2)))])/6

	#if derivative is non-zero at boundary
	#plt.plot(N_pos, b)
	b[0]   += -Bc0
	b[-1]  +=  Bc1
	sol = np.linalg.solve(A+A_p, b-double_deriv)

	#transfer back solution to regular basis
	for i in range(N):
		u += sol[i]*phi[i, :]
	return u, sol



ux0            = np.zeros((len(xi), len(k)), dtype="complex")
ux1            = np.zeros((len(xi), len(k)), dtype="complex")
uy1            = np.zeros((len(xi), len(k)), dtype="complex")
ux2            = np.zeros((len(xi), len(k)), dtype="complex")
B0             = np.zeros((len(xi), len(k)), dtype="complex")
B0_deriv       = np.zeros((len(xi), len(k)), dtype="complex")
B0_deriv_deriv = np.zeros((len(xi), len(k)), dtype="complex")
f              = np.zeros((len(xi), len(k)), dtype="complex")

#analytic results up to first order
ux0[:, np.argmin(abs(k-1))]  =         F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
ux0[:, np.argmin(abs(k+1))]  = np.conj(F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma))

ux1[:, np.argmin(abs(k-1))]  =        ((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p *xi)/np.cosh(kappa_p)) + (F0*np.tanh(gamma)/gamma)*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))/2
ux1[:, np.argmin(abs(k+1))]  = np.conj((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p *xi)/np.cosh(kappa_p)) + (F0*np.tanh(gamma)/gamma)*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))/2

uy1[:, np.argmin(abs(k-1))]  =        ((kappa*P_1*np.sinh(kappa)/(gamma*gamma))*(np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa)))/2
uy1[:, np.argmin(abs(k+1))]  = np.conj((kappa*P_1*np.sinh(kappa)/(gamma*gamma))*(np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa)))/2

ux2[:, np.argmin(abs(k-1))]  =        (P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))/2
ux2[:, np.argmin(abs(k+1))]  = np.conj(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))/2
ux2[:, np.argmin(abs(k-1))] -= scpi.trapz(ux2[:, np.argmin(abs(k-1))], xi)/2
ux2[:, np.argmin(abs(k+1))] -= scpi.trapz(ux2[:, np.argmin(abs(k+1))], xi)/2

B0[:, np.argmin(abs(k-1))]             =  (Pe*F0*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma)))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma))) + Pe*F0*np.tanh(gamma)/(2*gamma*gamma*gamma*rho*rho)
B0_deriv[:, np.argmin(abs(k-1))]       = (Pe*F0*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma)))*(np.sinh(rho*xi)/(np.sinh(rho)) - np.sinh(gamma*xi)/(np.sinh(gamma)))
B0_deriv_deriv[:, np.argmin(abs(k-1))] = (Pe*F0*np.tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma)))*(rho*np.cosh(rho*xi)/(np.sinh(rho)) - gamma*np.cosh(gamma*xi)/(np.sinh(gamma)))

B0[:, np.argmin(abs(k+1))]             = np.conj(B0[:, np.argmin(abs(k-1))])
B0_deriv[:, np.argmin(abs(k+1))]       = np.conj(B0_deriv[:, np.argmin(abs(k-1))])
B0_deriv_deriv[:, np.argmin(abs(k+1))] = np.conj(B0_deriv_deriv[:, np.argmin(abs(k-1))])

o1  = np.argmin(abs(k-1))
om1 = np.argmin(abs(k+1))

B_minus_deriv = np.zeros(np.shape(B_minus), dtype="complex")
B_plus_deriv  = np.zeros(np.shape(B_plus),  dtype="complex")
B_plus_deriv_deriv  = np.zeros(np.shape(B_plus),  dtype="complex")

term1 = np.zeros((N, len(k)), dtype="complex")
term2 = np.zeros((N, len(k)), dtype="complex")
term3 = np.zeros((N, len(k)), dtype="complex")
factor1 = np.zeros((len(xi), len(k)), dtype="complex")

for i in range(len(k)):
	factor1[:, i] = -uy1[:, i] + kappa*xi*ux1[:, i]

for i in range(len(k)):
	f[:, i]  = 0.5*kappa*kappa*xi*B0_deriv[:, i] + 0.5*(3+kappa*kappa*xi*xi)*B0_deriv_deriv[:, i] + Pe*ux2[:, i]
	if i != 0:
		f[:, i] += Pe*0.5*(ux1[:,  o1]*kappa*B_minus[:, i-1])

		for j in range(1, N-1):
			term1[j, i] += B_minus_coeff[j,   i-1]*(factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), o1] - factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), o1])/3
			term1[j, i] += B_minus_coeff[j+1, i-1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), o1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), o1])/6
			term1[j, i] -= B_minus_coeff[j-1, i-1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), o1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), o1])/6
		term1[0, i]  += (-B_minus_coeff[0,  i-1]+B_minus_coeff[1, i-1])*(factor1[np.argmin(abs(xi-(N_pos[0]))),  o1]   + 2*factor1[np.argmin(abs(xi-(N_pos[0]+Delta/2))),  o1])/6
		term1[-1, i] += (-B_minus_coeff[-2, i-1]+B_minus_coeff[-1,i-1])*(factor1[np.argmin(abs(xi-(N_pos[-1]))), o1]   + 2*factor1[np.argmin(abs(xi-(N_pos[-1]-Delta/2))), o1])/6
	
	if i != len(f[0,:])-1:
		f[:, i] += Pe*0.5*(ux1[:, om1]*kappa*B_minus[:, i+1])

		for j in range(1, N-1):
			term1[j, i] +=  B_minus_coeff[j,i+1]*(factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), om1] - factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), om1])/3
			term1[j, i] += B_minus_coeff[j+1, i+1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), om1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), om1])/6
			term1[j, i] -= B_minus_coeff[j-1, i+1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), om1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), om1])/6
		term1[0,  i] += (-B_minus_coeff[0,  i+1]+B_minus_coeff[1, i+1])*(factor1[np.argmin(abs(xi-(N_pos[0]))), om1]   + 2*factor1[np.argmin(abs(xi-(N_pos[0]+Delta/2))),  om1])/6
		term1[-1, i] += (-B_minus_coeff[-2, i+1]+B_minus_coeff[-1,i+1])*(factor1[np.argmin(abs(xi-(N_pos[-1]))), om1]  + 2*factor1[np.argmin(abs(xi-(N_pos[-1]-Delta/2))), om1])/6

for j in range(1, N-1):
	#term 2 and 3 have frequency equal to their
	term2[j, :] = -Delta*B_plus_coeff[j, :]/3 + 0.5*(N_pos[j]+Delta/3)*B_plus_coeff[j+1, :] - 0.5*(N_pos[j]-Delta/3)*B_plus_coeff[j-1, :] 
	term3[j, :] = (B_plus_coeff[j+1, :]+B_plus_coeff[j-1, :]-2*B_plus_coeff[j, :])/Delta 

term2[0,  :] = (B_plus_coeff[1,  :] - B_plus_coeff[0,  :])*Delta/3
term2[-1, :] = (B_plus_coeff[-1, :] - B_plus_coeff[-2, :])*Delta/3

term3[0,  :] = 2*(-B_plus_coeff[0,  :]+B_plus_coeff[1, :])/Delta
term3[-1, :] = 2*(-B_plus_coeff[-1, :]+B_plus_coeff[-2, :])/Delta

derivatives = term1 + term2 + term3
sol = np.zeros((len(xi), len(k)), dtype="complex")
BC0 = np.zeros(len(k))
BC1 = np.zeros(len(k))
sol_coeff = np.zeros((len(N_pos), len(k)))

for i in range(len(k)):
	if abs(k[i]) > 1e-4:
		sol[:,i], sol_coeff[:,i] = finite_element_solver(N, xi, f[:,i], derivatives[:, i], 1j*omega*k[i], BC0[i], BC1[i], False)
	else:
		#when looking at the source terms for n=0, we see that there are actually no contributions..., so we just take the homogeneous solution to satisfy the BCs
		sol[:,i] =  np.zeros(len(xi))

for i in range(int(len(k)/2+1)):
	plt.plot(xi, np.real(sol[:, i]+sol[:, -i-1]), label=r"$omega=$"+str(abs(k[i])))
plt.legend(loc="best")
plt.xlabel(r"Vertical position $xi$")
plt.ylabel(r"Second order Brenner field $B^{(2)}$")
plt.savefig("figures/second_order_B.pdf")
plt.show()

for i in range(len(k)):
	plt.plot(xi, np.imag(sol[:, i]+sol[:,-i-1]))
plt.show()

np.save("data/B2_new", sol)




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

b0    = B0
b0_deriv = B0_deriv
b0_deriv_deriv = B0_deriv_deriv 

for i in range(n):
	B1_min_deriv[:, i]        = interp1d(N_pos, np.gradient(B_minus_coeff[:,i], N_pos), kind='cubic')(xi)
	B1_plus_deriv[:, i]       = interp1d(N_pos, np.gradient(B_plus_coeff[:,i],  N_pos), kind='cubic')(xi)
	B1_plus_deriv_deriv[:, i] = interp1d(N_pos, np.gradient(np.gradient(B_plus_coeff[:,i], N_pos), N_pos), kind='cubic')(xi)
	B2_0_deriv[:, i]          = interp1d(N_pos, np.gradient(sol_coeff[:,i], N_pos),     kind='cubic')(xi)


for i in range(len(k)):
	D_eff_xi[:, np.argmin(abs(new_k -k[i]))] += (kappa*xi*B1_min_deriv[:, i] + kappa*B_minus[:, i])
	for j in range(len(k)):
		D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_plus[:,j]*B_plus[:,i]    + B1_plus_deriv[:, i]*B1_plus_deriv[:, j])
		D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5* (kappa*kappa*B_minus[:,j]*B_minus[:,i]  + B1_min_deriv[:, i] *B1_min_deriv[:, j])
		D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += b0_deriv[:, i]*(- B1_plus_deriv[:,j] - kappa*kappa*xi*B_plus[:, j] + 2*B2_0_deriv[:,j] )
		D_eff_xi[:, np.argmin(abs(new_k - (k[i]+k[j])))] += 0.5*(1+kappa*kappa*xi*xi)*b0_deriv[:, i]*b0_deriv[:, j]

for i in range(len(new_k)):
	D_eff[i] = scpi.trapz(D_eff_xi[:, i], xi)/2
	total_D += np.exp(1j*omega*new_k[i]*t)*D_eff[i]

#print("HERE:", np.max(np.imag(total_D)))

D_parallels = scpi.trapz(np.real(total_D), t)/(2*np.pi/(omega))
print("D_eff = ", D_parallels)

plt.figure(1)
plt.plot(new_k, np.real(D_eff))
plt.xlabel(r"frequency [$omega$]", fontsize=12)
plt.ylabel(r"Amplitude", fontsize=12)
plt.savefig("figures/fourier_D.pdf")

plt.figure(2)
plt.plot(t, np.real(total_D))
plt.xlabel(r"time [$t$]", fontsize=12)
plt.ylabel(r"Brenner field", fontsize=12)
plt.savefig("figures/Brenner_field_vs_t.pdf")

np.save("data/D_parallels_kappa_D1", D_parallels)
plt.show()