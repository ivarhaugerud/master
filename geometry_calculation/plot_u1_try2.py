import numpy as np 
import matplotlib.pyplot as plt 

Sc = 1
F0 = 1
pi = np.pi 
omega = 1
kappa = 1
Nt = 4
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(0, 2*np.pi/kappa, 300)
xi = np.linspace(-1, 1, 299)
epsilon = 0.3

gamma = np.sqrt(1j*omega/Sc)

u_x  = np.zeros((len(T), len(xi), len(eta)), dtype="complex")
u_y  = np.zeros((len(T), len(xi), len(eta)), dtype="complex")
u_x2 = np.zeros((len(T), len(xi), len(eta)), dtype="complex")
u_y2 = np.zeros((len(T), len(xi), len(eta)), dtype="complex")

kappa_prime = np.sqrt(1j*omega/Sc + kappa*kappa)

psi = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
Ay = kappa*np.sinh(kappa)*psi/(gamma*gamma*np.sinh(kappa_prime))
Bx = -kappa_prime*Ay/kappa

for t in range(len(T)):
	for x in range(len(xi)):
		u_x2[t, x, :] = -F0*xi[x]*np.sinh(gamma*xi[x])*np.sin(kappa*eta)/(gamma*np.cosh(gamma))
		u_x2[t, x, :] += kappa*np.cosh(kappa*xi[x])*psi*np.sin(kappa*eta)/(gamma*gamma)
		u_x2[t, x, :] += np.cosh(kappa_prime*xi[x])*Bx*np.sin(kappa*eta)

		u_x2[t,x,:] *= epsilon 
		u_x2[t,x,:] += F0*(1-np.cosh(gamma*xi[x])/np.cosh(gamma))/(gamma*gamma)
		u_x2[t, x, :] *= np.exp(1j*omega*T[t])

		u_y2[t, x, :] += -Sc*kappa*np.sinh(kappa*xi[x])*psi*np.cos(kappa*eta)/(1j*omega)
		u_y2[t, x, :] += np.sinh(kappa_prime*xi[x])*Ay*np.cos(kappa*eta)
		u_y2[t, x, :] *= np.exp(1j*omega*T[t])

P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))

for t in range(len(T)):
	for x in range(len(xi)):
		u_x[t,x,:]    = epsilon*np.sin(kappa*eta)*(  (P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi[x])/np.cosh(kappa) - np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime) - xi[x]*np.sinh(gamma*xi[x])/np.sinh(gamma))  )
		u_x[t,x,:]   += F0*(1-np.cosh(gamma*xi[x])/np.cosh(gamma))/(gamma*gamma)
		u_x[t,x,:]   *= np.exp(1j*omega*T[t])

		u_y[t, x, :] = (np.exp(1j*omega*T[t])*np.cos(kappa*eta)*kappa*P1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_prime*xi[x])/np.sinh(kappa_prime) - np.sinh(kappa*xi[x])/np.sinh(kappa))

#test boundary conditions

for i in range(len(T)):
	plt.plot(xi, np.real(u_y[i, :, 150]))
	plt.plot(xi, np.real(u_y2[i, :, 150]))
	plt.show()

for i in range(len(T)):
	plt.plot(xi, np.real(u_x[i, :, 150]))
	plt.plot(xi, np.real(u_x2[i, :, 150]), "--")
	plt.show()

for i in range(len(T)):
	plt.clf()
	plt.title(r"$T=$"+str(T[i])[:4] + r"$[\omega^{-1}]$")
	X, Y = np.meshgrid(eta, xi)
	Z = np.sqrt( (u_x[i,:,:]*u_x[i,:,:]) + u_y[i,:,:]*u_y[i,:,:])

	cp = plt.contourf(X, Y, Z)#, levels=np.linspace(0, 0.7, 40))
	cb = plt.colorbar(cp)

	#lw = 2*Z/max_vel
	plt.streamplot(X, Y, np.real(u_x[i,:,:]), np.real(u_y[i,:,:]), density=0.5, color='k')#, linewidth=lw)

	plt.axis([min(eta), max(eta), -1.01, 1.01])

	plt.draw()
	plt.pause(0.1)
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	plt.savefig("figures/streamplot_epsilon1_nr"+str(i)+".pdf")
plt.show()	


