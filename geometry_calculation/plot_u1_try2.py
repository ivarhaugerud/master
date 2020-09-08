import numpy as np 
import matplotlib.pyplot as plt 

Sc = 1
F0 = 1
pi = np.pi 
omega = 1
kappa = 1
Nt = 4
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(0, 6*np.pi/kappa, 300)
xi = np.linspace(-1, 1, 299)

gamma = np.sqrt(1j*omega/Sc)

u_x = np.zeros((len(T), len(xi), len(eta)), dtype="complex")
u_y = np.zeros((len(T), len(xi), len(eta)), dtype="complex")

kappa_prime = np.sqrt(1j*omega/Sc + kappa*kappa)
Ay          = (-F0*np.tanh(gamma)*np.sinh(kappa)*np.cosh(kappa)/(gamma*np.sinh(kappa_prime)*np.cosh(kappa_prime)*np.cosh(kappa_prime)))/(1+np.sinh(kappa)*np.cosh(kappa)*kappa_prime/(np.sinh(kappa_prime)*np.cosh(kappa_prime)*kappa))
psi_1       = gamma*gamma*np.sinh(kappa_prime)*Ay/(kappa*np.sinh(kappa))
#Bx          = -kappa_prime*Ay/kappa 
Bx = F0*np.tanh(gamma)/(gamma*np.cosh(kappa_prime)) + kappa*np.cosh(kappa)*psi_1/(gamma*gamma*np.cosh(kappa_prime))

for t in range(len(T)):
	for x in range(len(xi)):
		u_x[t, x, :] = -F0*xi[x]*np.sinh(gamma*xi[x])*np.sin(kappa*eta)/(gamma*np.cosh(gamma))
		u_x[t, x, :] += Sc*kappa*np.cosh(kappa*xi[x])*psi_1*np.sin(kappa*eta)/(1j*omega)
		u_x[t, x, :] += np.cosh(kappa_prime*xi[x])*Bx*np.sin(kappa*eta)
		u_x[t, x, :] *= np.exp(1j*omega*T[t])

for i in range(len(T)):
	#plt.imshow(np.real(u_x[i, :,:]))
	plt.plot(xi, np.real(u_x[i, :, 0]))
	plt.plot(xi, np.real(u_x[i, :, 100]))
	plt.plot(xi, np.real(u_x[i, :, 200]))
	plt.plot(xi, np.real(u_x[i, :, 299]))
	plt.show()


### plot u_y
for t in range(len(T)):
	for x in range(len(xi)):
		u_y[t, x, :] += -Sc*kappa*np.sinh(kappa*xi[x])*psi_1*np.cos(kappa*eta)/(1j*omega)
		u_y[t, x, :] += np.sinh(kappa_prime*xi[x])*Ay*np.cos(kappa*eta)
		u_y[t, x, :] *= np.exp(1j*omega*T[t])


for i in range(len(T)):
	plt.plot(xi, np.real(u_y[i, :, 0]))
	plt.plot(xi, np.real(u_y[i, :, 100]))
	plt.plot(xi, np.real(u_y[i, :, 200]))
	plt.plot(xi, np.real(u_y[i, :, 299]))
	"""
	plt.clf()
	plt.imshow(np.real(u_y[i, :,:]))
	plt.pause(0.01)
	plt.draw()
	"""

plt.show()







