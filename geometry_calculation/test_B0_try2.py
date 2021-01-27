import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci 

Pe = 6
nu = 1
omega = 2*np.pi/3
pi = np.pi
xi = np.linspace(-1, 1, int(1e4))
M = int(30)
Nu       = np.logspace(-3, 3, 12)
Ds       = np.logspace(-3, 3, 13)
#paramters = np.zeros((2, len(Omega)))
#paramters[0, :] = Omega
#paramters[1, :] = Schmidt

D_eff   = np.zeros((len(Nu), len(Ds), 2), dtype="complex")

for s in range(len(Ds)):
	D = Ds[s]
	Pe = 1/D
	for o in range(len(Nu)):
		F0 = 3#/Nu[o] 
		Sc = Nu[o]
		gamma   = np.sqrt(1j*omega/Sc)
		gamma_c = np.conj(gamma)
		rho     = np.sqrt(1j*omega/D)
		rho_c   = np.conj(rho)

		D_ana = ((F0*F0*np.tanh(gamma)*np.tanh(gamma_c))/(gamma*gamma_c*(rho*rho-gamma*gamma)*(rho_c*rho_c-gamma_c*gamma_c)))*( (rho/np.tanh(rho) - rho_c/np.tanh(rho_c) )/(2*rho*rho) + (gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c) )/(2*gamma*gamma) + (rho_c/np.tanh(rho_c) - rho/np.tanh(rho) + gamma_c/np.tanh(gamma_c) - gamma/np.tanh(gamma) )/(gamma*gamma+rho*rho))
		D_eff[o, s, 0] = D_ana
		
		"""
		omega = omega/D
		a_m_squared = np.zeros(M, dtype="complex")
		for m in range(1, M):
			for k_prime in range(M):
				k = 2*k_prime + 1
				alpha_k = k*k*pi*pi*Sc/4

				for n_prime in range(M):
					n = 2*n_prime + 1
					alpha_n = n*n*pi*pi*Sc/4
					a_m_squared[m] += (alpha_n-1j*omega)*(alpha_k+1j*omega)/((alpha_n*alpha_n + omega*omega)*(alpha_k*alpha_k + omega*omega)*(m*m-n*n/4)*(m*m-k*k/4))

		g_tilde = 0
		for m in range(M):
			g_tilde += m*m*a_m_squared[m]/(m*m*m*m*pi*pi*pi*pi + omega*omega)

		g_tilde *= 16*Pe*Pe*Sc*Sc*F0*F0/(pi*pi)
		D_eff[o, s, 1] = 0.5*g_tilde #0.5 due to decomposisition of u in eq 2.6 in Bruus
		"""
		#print(np.real(D_eff[o, s, 0]), np.real(D_eff[o, s, 1]))


fig = plt.figure(1)
x_, y_ = np.meshgrid(Nu, Ds)

ax1 = plt.contourf(x_,y_, np.transpose((105/2*D_eff[:,:,0])) )#, levels=np.linspace(0, 1, 101))
cbar = fig.colorbar(ax1)
#cbar.ax.locator_params(nbins=10)
cbar.ax.set_ylabel(r'Geometric factor $\tilde{g}$', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r"Diffusion coefficient $D$", fontsize=14)
plt.xlabel(r"Kinematic viscosity $\nu$", fontsize=14)
plt.savefig("figures/D_0_eff.pdf", bbox_inches="tight")
#plt.savefig("../figures/g_tilde_D0.pdf")
plt.show()


print(Ds[3])
plt.plot(Nu, 105/2 * np.real(D_eff[:, 3, 0]), label="mine")
plt.plot(Nu, 105/2 * np.real(D_eff[:, 3, 1]), label="bruus")
plt.xscale("log")
plt.legend(loc="best")
plt.show()

print(Ds[6])
plt.plot(Nu, 105/2 * np.real(D_eff[:, 6, 0]), label="mine")
plt.plot(Nu, 105/2 * np.real(D_eff[:, 6, 1]), label="brus")
plt.xscale("log")
plt.legend(loc="best")
plt.show()

print(Ds[1])
plt.plot(Nu, 105/2 * np.real(D_eff[:, 1, 0]), label="mine")
plt.plot(Nu, 105/2 * np.real(D_eff[:, 1, 1]), label="brus")
plt.xscale("log")
plt.legend(loc="best")
plt.show()

print(Ds[2])
plt.plot(Nu, 105/2 * np.real(D_eff[:, 1, 0]), label="mine")
plt.plot(Nu, 105/2 * np.real(D_eff[:, 1, 1]), label="brus")
plt.xscale("log")
plt.legend(loc="best")
plt.show()


np.save("data_test/D_eff_2", D_eff)
