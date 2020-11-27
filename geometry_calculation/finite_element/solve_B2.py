import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scpi 

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

B_plus  = np.load("data/B_plus.npy")
B_minus = np.load("data/B_minus.npy")
k       = np.load("data/k.npy")

N = 500
xi = np.linspace(-1, 1, len(B_plus[:,0]))

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

B0[:, np.argmin(abs(k-1))]             =  (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma))))/2
B0_deriv[:, np.argmin(abs(k-1))]       = ((Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma)))/2
B0_deriv_deriv[:, np.argmin(abs(k-1))] = ((Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma)))/2

B0[:, np.argmin(abs(k-1))]             =  np.conj(Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma))))/2
B0_deriv[:, np.argmin(abs(k-1))]       = np.conj((Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma)))/2
B0_deriv_deriv[:, np.argmin(abs(k-1))] = np.conj((Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma)))/2

o1  = np.argmin(abs(k-1))
om1 = np.argmin(abs(k-1))

from scipy.signal import savgol_filter
for i in range(len(k)):
	B_minus[:,i] = savgol_filter(B_minus[:,i], 1001, 5) # window size 51, polynomial order 3
	B_plus[:, i] = savgol_filter(B_plus[:, i], 1001, 5) # window size 51, polynomial order 3

for i in range(len(k)):
	f[:, i]  = 0.5*kappa*kappa*xi*B0_deriv[:, i] + 0.5*(3+kappa*kappa*xi*xi)*B0_deriv_deriv[:, i] - 0.5*kappa*kappa*xi*np.gradient(B_plus[:, i], xi) - np.gradient(savgol_filter(np.gradient(B_plus[:,i], xi), 1001, 5), xi) + Pe*ux2[:, i]
	if i != 0:
		f[:, i] += Pe*0.5*(ux0[:,  o1]*kappa*xi*np.gradient(B_minus[:, i-1]) + ux1[:,  o1]*kappa*B_minus[:, i-1] - uy1[:, o1]*np.gradient(B_minus[:, i-1]) )
	if i != len(f[0,:])-1:
		f[:, i] += Pe*0.5*(ux0[:, om1]*kappa*xi*np.gradient(B_minus[:, i+1]) + ux1[:, om1]*kappa*B_minus[:, i+1] - uy1[:, om1]*np.gradient(B_minus[:, i+1]))

sol = np.zeros((len(xi), len(k)))
for i in range(len(k)):
	print(i)
	if abs(k[i]) > 1e-4:
		sol[:,i] = finite_element_solver(N, xi, f[:,i], 1j*omega*k[i], 0, 0)

sol[:, np.argmin(abs(k-0))] = -kappa*kappa*xi*xi/8

plt.plot(xi[:-1], np.real(sol[:-1]), label="numerical solution")
plt.xlabel("x")
plt.show()