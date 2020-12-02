import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scpi 
from scipy.signal import savgol_filter

B_plus  = np.load("data/B_plus.npy")
B_minus = np.load("data/B_minus.npy")
B_2     = np.load("data/B2_new.npy")
k       = np.load("data/k.npy")
xi = np.linspace(-1, 1, len(B_plus[:,0]))

#system parameters
kappa = 1
Sc = 1.2
omega = 0.5
F0 = 3
Pe = F0*Sc


#implicitly defined parameters
gamma   = np.sqrt(1j*omega/Sc)
kappa_p = np.sqrt(gamma*gamma + kappa*kappa)
rho     = np.sqrt(1j*omega)
P_1     = (F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))
B0_deriv      = np.zeros((len(xi), len(k)), dtype="complex")
B_plus_grad   = np.zeros((len(xi), len(k)), dtype="complex")
B_minus_grad  = np.zeros((len(xi), len(k)), dtype="complex")
B_2_grad      = np.zeros((len(xi), len(k)), dtype="complex")

t        = np.linspace(0, 2*np.pi/omega, int(1e4))
new_k    = np.arange(-2*max(k), 2*max(k)+1e-3, 1)
D_eff_xi = np.zeros((len(xi), len(new_k)), dtype="complex")
total_D  = np.zeros(len(t),                dtype="complex")
D_eff    = np.zeros(len(new_k),            dtype="complex")

for i in range(len(k)):
	B_plus_grad[:, i]  = np.gradient(np.real(B_plus[:,  i]), xi) +1j*np.gradient(np.imag(B_plus[:,  i]), xi)
	B_minus_grad[:, i] = np.gradient(np.real(B_minus[:, i]), xi) +1j*np.gradient(np.imag(B_minus[:, i]), xi)
	B_2_grad[:, i]     = np.gradient(np.real(B_2[:, i]),     xi) +1j*np.gradient(np.imag(B_2[:, i]),     xi)

B0_deriv[:, np.argmin(abs(k-1))] =         (Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma))/2
B0_deriv[:, np.argmin(abs(k+1))] = np.conj((Pe*F0*np.tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(np.sinh(rho*xi)/np.sinh(rho) - np.sinh(gamma*xi)/np.sinh(gamma)))/2

for i in range(len(k)):
	D_eff_xi[:, np.argmin(abs(new_k-k[i]))] -=  kappa*B_minus[:, i]/2 + kappa*xi*B_minus_grad[:,i]/2

for i in range(len(k)):
	for j in range(len(k)):
		D_eff_xi[:, np.argmin(abs(new_k -k[i]-k[j]))] += B0_deriv[:,i]*B0_deriv[:,j]*(5+kappa*kappa*xi*xi)/2
		D_eff_xi[:, np.argmin(abs(new_k -k[i]-k[j]))] += (kappa*kappa*B_plus[:, i]*B_plus[:, j] + B_plus_grad[:,i]*B_plus_grad[:,j] + kappa*kappa*B_minus[:, i]*B_minus[:, j] + B_minus_grad[:,i]*B_minus_grad[:,j])/2
		D_eff_xi[:, np.argmin(abs(new_k -k[i]-k[j]))] += 2*B0_deriv[:, i]*B_2_grad[:,j] + B0_deriv[:,i]*(B_plus_grad[:,j] - kappa*kappa*xi*B_plus[:, j])

for i in range(len(new_k)):
	D_eff[i] = scpi.trapz(D_eff_xi[:, i], xi)/2
	total_D += np.exp(1j*omega*new_k[i]*t)*D_eff[i]


plt.plot(new_k, np.imag(D_eff), "ro")
plt.plot(new_k, -np.flip(np.imag(D_eff)), "r--")
plt.show()
plt.plot(new_k, np.real(D_eff), "ko")
x = np.linspace(0, max(new_k), int(1e4))
plt.show()

plt.plot(t, np.real(total_D))
plt.plot(t, np.imag(total_D), "--")
plt.show()

D_eff = scpi.trapz(np.real(total_D), t)/(max(t)-min(t))
print(D_eff)