import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scpi 
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

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

	if laplace:
		A_p[0,0] += 1
		A_p[-1, -1] += 1

		#Bc0 *= 1/Delta_x
		#Bc1 *= 1/Delta_x

	A_p *= 1/Delta_x
	A   *= a*Delta_x/6 # a=k^2

	#calculate source vector
	for i in range(len(b)):
		b[i] = -Delta_x*(2*f[np.argmin(abs(x-(N_pos[i]+Delta_x/2)))] + 2*f[np.argmin(abs(x-(N_pos[i])))] + 2*f[np.argmin(abs(x-(N_pos[i]-Delta_x/2)))])/6

	b[0]  = -Delta_x*(f[np.argmin(abs(x-(N_pos[0 ])))] + 2*f[np.argmin(abs(x-(N_pos[0] + Delta_x/2)))])/6
	b[-1] = -Delta_x*(f[np.argmin(abs(x-(N_pos[-1])))] + 2*f[np.argmin(abs(x-(N_pos[-1]- Delta_x/2)))])/6

	#if derivative is non-zero at boundary
	b[0]   +=  Bc0
	b[-1]  += -Bc1
	sol = np.linalg.solve(A+A_p, b+double_deriv)

	#transfer back solution to regular basis
	for i in range(N):
		u += sol[i]*phi[i, :]

	return u

B_plus  = np.load("data/B_plus.npy")
B_plus_coeff  = np.load("data/B_plus_coeff.npy")

B_minus = np.load("data/B_minus.npy")
B_minus_coeff  = np.load("data/B_minus_coeff.npy")
k       = np.load("data/k.npy")

N = len(B_minus_coeff[:, 0])
xi = np.linspace(-1, 1, len(B_plus[:,0]))
N_pos = np.linspace(-1, 1, N)
Delta = N_pos[1]-N_pos[0]

#system parameters
kappa = 2
Sc = 1.2
omega = 1
F0 = 3
Pe = F0*Sc

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

B0[:, np.argmin(abs(k+1))]             =  np.conj(Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(np.cosh(rho*xi)/(rho*np.sinh(rho)) - np.cosh(gamma*xi)/(gamma*np.sinh(gamma))))/2
B0_deriv[:, np.argmin(abs(k+1))]       = np.conj((Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma)))/2
B0_deriv_deriv[:, np.argmin(abs(k+1))] = np.conj((Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(rho*np.cosh(rho*xi)/np.sinh(rho) - gamma*np.cosh(gamma*xi)/np.sinh(gamma)))/2

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
	factor1[:, i] = uy1[:, i] -kappa*xi*ux1[:, i]

for i in range(len(k)):
	f[:, i]  = 0.5*kappa*kappa*xi*B0_deriv[:, i] + 0.5*(3+kappa*kappa*xi*xi)*B0_deriv_deriv[:, i] + Pe*ux2[:, i]
	
	if i != 0:
		f[:, i] += Pe*0.5*(ux1[:,  o1]*kappa*B_minus[:, i-1])

		for j in range(1, N-1):
			term1[j, i] +=  B_minus_coeff[j,   i-1]*(factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), o1] - factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), o1])/3
			term1[j, i] += B_minus_coeff[j+1, i-1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), o1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), o1])/6
			term1[j, i] -= B_minus_coeff[j-1, i-1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), o1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), o1])/6
		term1[0, i]  += (-B_minus_coeff[0,  i-1]+B_plus_coeff[1, i-1])*(factor1[np.argmin(abs(xi-(N_pos[0]))),  o1]   + 2*factor1[np.argmin(abs(xi-(N_pos[0]+Delta/2))),  o1])/6
		term1[-1, i] += (-B_minus_coeff[-2, i-1]+B_plus_coeff[-1,i-1])*(factor1[np.argmin(abs(xi-(N_pos[-1]))), o1]   + 2*factor1[np.argmin(abs(xi-(N_pos[-1]-Delta/2))), o1])/6
	
	if i != len(f[0,:])-1:
		f[:, i] += Pe*0.5*(ux1[:, om1]*kappa*B_minus[:, i+1])

		for j in range(1, N-1):
			term1[j, i] +=  B_minus_coeff[j,i+1]*(factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), om1] - factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), om1])/3
			term1[j, i] += B_minus_coeff[j+1, i+1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), om1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]+Delta/2))), om1])/6
			term1[j, i] -= B_minus_coeff[j-1, i+1]*(factor1[np.argmin(abs(xi-(N_pos[j]))), om1] + 2*factor1[np.argmin(abs(xi-(N_pos[j]-Delta/2))), om1])/6
		term1[0,  i] += (-B_minus_coeff[0,  i+1]+B_plus_coeff[1, i+1])*(factor1[np.argmin(abs(xi-(N_pos[0]))), om1]   + 2*factor1[np.argmin(abs(xi-(N_pos[0]+Delta/2))),  om1])/6
		term1[-1, i] += (-B_minus_coeff[-2, i+1]+B_plus_coeff[-1,i+1])*(factor1[np.argmin(abs(xi-(N_pos[-1]))), om1]  + 2*factor1[np.argmin(abs(xi-(N_pos[-1]-Delta/2))), om1])/6

for j in range(1, N-1):
	#term 2 and 3 have frequency equal to their
	term2[j, :] = -Delta*B_plus_coeff[j, :]/3 + 0.5*(N_pos[j]+Delta/3)*B_plus_coeff[j+1, :] - 0.5*(N_pos[j]-Delta/3)*B_plus_coeff[j-1, :] 
	term3[j, :] = (B_plus_coeff[j+1, :]+B_plus_coeff[j-1, :]-2*B_plus_coeff[j, :])/Delta 

term2[0,  :] = (B_plus_coeff[1,  :] - B_plus_coeff[0,  :])*Delta/3
term2[-1, :] = (B_plus_coeff[-1, :] - B_plus_coeff[-2, :])*Delta/3

term3[0,  :] = (-B_plus_coeff[0,  :]+B_plus_coeff[1, :])/Delta
term3[-1, :] = (-B_plus_coeff[-1, :]+B_plus_coeff[-2, :])/Delta

derivatives = term1 + term2 + term3
sol = np.zeros((len(xi), len(k)), dtype="complex")
BC0 = np.zeros(len(k))
BC1 = np.zeros(len(k))

BC0[np.argmin(abs(k-0))] = -kappa*kappa/4
BC1[np.argmin(abs(k-0))] =  kappa*kappa/4

for i in range(len(k)):
	if abs(k[i]) > 1e-4:
		sol[:,i] = finite_element_solver(N, xi, f[:,i], derivatives[:, i], 1j*omega*k[i], BC0[i], BC1[i], False)
	else:
		sol[:,i] = finite_element_solver(N, xi, f[:,i], derivatives[:, i], 0,             BC0[i], BC1[i], True)
		plt.plot(xi, np.real(sol[:, i]))
		print(BC0[i], BC1[i])
		plt.show()
		plt.plot(xi, np.gradient(np.real(sol[:, i]), xi))
		plt.show()

for i in range(len(k)):
	plt.plot(xi, np.real(sol[:, i]))
plt.show()

for i in range(len(k)):
	plt.plot(xi, np.imag(sol[:, i]+sol[:,-i-1]))
plt.show()

np.save("data/B2", sol)

