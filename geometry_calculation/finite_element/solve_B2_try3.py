import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
import scipy.integrate as scpi 

def coupled_finite_element_solver(N, x, alpha, f):
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

	#if derivative is non-zero at boundary
	#b[N*j+0]    -= Bc0[j]
	#b[N*j+N-1]  += Bc1[j]

	sol = np.linalg.solve(A+A_p, b)
	#transfer back solution to regular basis
	for i in range(N):
		u[:] += sol[i]*phi[i, :]
	return u

B_plus  	   = np.load("data/B_plus.npy")
B_minus 	   = np.load("data/B_minus.npy")
k       	   = np.load("data/k.npy")
B_minus_coeff  = np.load("data/B_minus_coeff.npy")
B_plus_coeff   = np.load("data/B_plus_coeff.npy")

tol   = 1e-6
k     = np.array([-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
N 	  = len(B_plus_coeff[:, 0])
N_pos = np.linspace(-1, 1, N)
xi    = np.linspace(-1, 1, int(1e5))
n     = len(k) #number of vectors
Delta = N_pos[1]-N_pos[0]

#system parameters
omega = 5/(2*np.pi)
nu = 16
Sc = nu
F0 = 10
U = 1
Pe = 3
kappa = 10

#implicitly defined parameters
gamma   = np.sqrt(1j*omega/Sc)
kappa_p = np.sqrt(gamma*gamma + kappa*kappa)
rho     = np.sqrt(1j*omega)
P_1     = (F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))
p_np2   = 1j*omega*k + kappa*kappa

#analytic results up to first order
ux0  = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
ux1  = ((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p *xi)/np.cosh(kappa_p)) + (F0*np.tanh(gamma)/gamma)*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))/2
uy1  = ((kappa*P_1*np.sinh(kappa)/(gamma*gamma))*(np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa)))/2
ux2  = (P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))/2
ux2 -= scpi.trapz(ux2, xi)/2

B0             = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma))) + Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*rho*rho)
B0_deriv       = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))
B0_deriv_deriv = (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma))

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

for i in range(int(n/2)+1):
	plt.plot(xi, q[i,:]+q[-1-i,:], label=(str(k[i])))
plt.legend(loc="best")
plt.show()
#works for differential equation with constant terms, now just need coupeling to work as well
sol = np.zeros((len(xi), n), dtype="complex")

for i in range(n):
	if abs(k[i]) > 1e-4:
		sol[:,i] = coupled_finite_element_solver(N, xi, 1j*omega*k[i], -q[i,:])

for i in range(int(n/2)+1):
	plt.plot(xi, sol[:,i]+sol[:,-1-i], label=(str(k[i])))
plt.legend(loc="best")
plt.show()
for i in range(n):
	if np.max(abs(np.imag(sol[i,:]+sol[-i-1,:])))/2 > tol:
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
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] + kappa*ux0*f0_g1[Z-1, :]/2 + kappa*np.conj(ux0)*f0_g1[Z+1, :]/2)
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)


Z = np.argmin(abs(k-1))
plt.figure(1)
plt.title(str(1))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] - kappa*ux0*f0_g1[Z-1, :]/2 - kappa*np.conj(ux0)*f0_g1[Z+1, :]/2)
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

Z = np.argmin(abs(k-2))
plt.figure(2)
plt.title(str(2))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] + kappa*ux0*f0_g1[Z-1, :]/2 + kappa*np.conj(ux0)*f0_g1[Z+1, :]/2)
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)


Z = np.argmin(abs(k-3))
plt.figure(3)
plt.title(str(3))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] - kappa*ux0*f0_g1[Z-1, :]/2 - kappa*np.conj(ux0)*f0_g1[Z+1, :]/2)
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

Z = np.argmin(abs(k-4))
plt.figure(4)
plt.title(str(4))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] + kappa*ux0*f0_g1[Z-1, :]/2 + kappa*np.conj(ux0)*f0_g1[Z+1, :]/2)
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

Z = np.argmin(abs(k-5))
plt.figure(5)
plt.title(str(5))
RHS = np.real(p_np2[Z]*f0_g1[Z, :] -q[Z, :] - kappa*ux0*f0_g1[Z-1, :]/2 - kappa*np.conj(ux0)*f0_g1[Z+1, :]/2)
plt.plot(N_pos, np.real(np.gradient(np.gradient(coeff_f0g1[N*Z:N*(1+Z)], N_pos), N_pos)))
plt.plot(xi, RHS)

plt.show()

print(A)
#reset source term for new solution
q = np.zeros((len(k), len(xi)), dtype="complex")
couple_backward *= -1
couple_forward  *= -1

#boundary conditions
Bc0 = np.zeros(n, dtype="complex")
Bc1 = np.zeros(n, dtype="complex")

sol, coeff_g0f1 = coupled_finite_element_solver(N, n, xi, p_np2, couple_backward, couple_forward, q, Bc0, Bc1, k)

for i in range(int(len(k)/2)+1):
	if np.max(abs(np.imag(sol[i,:]+sol[-i-1,:]) - np.imag(sol[i,0]+sol[-i-1,0])))/2 > tol:
		print("LARGE IMAGINARY VALUE FOR " + str(k[i]) + "-OMEGA = ", np.max(abs(np.imag(sol[i,:]+sol[-i-1,:]) - np.imag(sol[i,0]+sol[-i-1,0]))))

sol = np.zeros(np.shape(sol))
g0_f1 = np.zeros(np.shape(sol))
coeff_g0f1 = np.zeros(np.shape(coeff_f0g1))

B_plus   = np.zeros((len(xi), len(k)), dtype="complex")  #sin-solution
B_minus  = np.zeros((len(xi), len(k)), dtype="complex")  #cos-solution
B_plus_coeff    = np.zeros((N, len(k)), dtype="complex") #sin-solution
B_minus_coeff   = np.zeros((N, len(k)), dtype="complex") #sin-solution

for i in range(int(len(k)/2)+1):
	if abs(k[i]) % 2 == 0:
		print("even:", k[i], "and", k[-i-1])
		plt.figure(1) #cos figure
		plt.plot(xi, np.real(f0_g1[i,:]+f0_g1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))
		plt.figure(2) #sin figure
		plt.plot(xi, np.real(g0_f1[i,:]+g0_f1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))

	else:
		print("odd:", k[i], "and", k[-i-1])
		plt.figure(2) #sin figure
		plt.plot(xi, np.real(f0_g1[i,:]+f0_g1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))
		plt.figure(1) #cos figure
		plt.plot(xi, np.real(g0_f1[i,:]+g0_f1[-i-1,:]), label=r"$\omega=\omega$"+str(abs(k[i])))	

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