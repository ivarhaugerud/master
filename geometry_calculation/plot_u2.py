import numpy as np 
import matplotlib.pyplot as plt 
from cmath import *
from scipy.interpolate import griddata

Ny = int(150)
Nx = 100
Nt = 10

epsilon = 0.3
omega = 1e-2
Sc = 1
F_0 = 1
kappa = 1
gamma = 1j*omega/Sc 
F0 = 3

kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
P_1 = (gamma*F_0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))
psi_2 = (F0*kappa_p*tanh(gamma)*tanh(kappa)/(2*gamma)-kappa*tanh(kappa_pp)*(F0/2 + P_1*sinh(kappa))/(kappa_pp))/(2*kappa*sinh(2*kappa)/(gamma*gamma) - 4*kappa*kappa*cosh(2*kappa)*tanh(kappa_pp)/(gamma*gamma*kappa_pp))
Ay = (2*kappa*psi_2*np.sinh(2*kappa)/(gamma*gamma) - F0*kappa_p*np.tanh(gamma)*np.tanh(kappa)/(2*gamma))/np.sinh(kappa_pp)
Ax = kappa_pp*Ay/(2*kappa)
"""
psi_2 = ((F0*gamma**2 - 2*P_1*kappa**2*np.sinh(kappa) + 2*P_1*kappa_p**2*np.sinh(kappa))*np.sinh(kappa_p)*np.sinh(np.sqrt(gamma**2 + 4*kappa**2)))/(4*(2*kappa*np.sinh(np.sqrt(gamma**2 + 4*kappa**2))*np.cosh(2*kappa) - np.sqrt(gamma**2 + 4*kappa**2)*np.sinh(2*kappa)*np.cosh(np.sqrt(gamma**2 + 4*kappa**2)))*np.sinh(kappa_p))
Ay = kappa*(F0*gamma**2*sinh(kappa)*cosh(kappa) + P_1*kappa**2*cosh(kappa) - 2*P_1*kappa*kappa_p*sinh(kappa)**3/tanh(kappa_p) - P_1*kappa*kappa_p*sinh(kappa)/tanh(kappa_p) + 2*P_1*kappa_p**2*sinh(kappa)**2*cosh(kappa))/(2*gamma**2*(2*kappa*sinh(kappa)**2*sinh(sqrt(gamma**2 + 4*kappa**2)) + kappa*sinh(sqrt(gamma**2 + 4*kappa**2)) - sqrt(gamma**2 + 4*kappa**2)*sinh(kappa)*cosh(kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))) 
Ax = kappa_pp*Ay/(2*kappa)
psi_2 = ((F0*gamma**2 - 2*P_1*kappa**2*sinh(kappa) + 2*P_1*kappa_p**2*sinh(kappa))*sinh(kappa_p)*sinh(sqrt(gamma**2 + 4*kappa**2)))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*sinh(kappa_p))
"""

T = np.linspace(0, 2*np.pi/omega, Nt)
eta = np.linspace(0, 4*pi/kappa, Nx)
xi = np.linspace(-1, 1, Ny)
u = np.zeros((2, Ny, Nx, Nt))

uy = P_1*kappa*kappa_p*np.sinh(kappa)*xi*np.cosh(kappa_p*xi)/(2*gamma*gamma*np.sinh(kappa_p)) - P_1*kappa*kappa*xi*np.cosh(kappa*xi)/(2*gamma*gamma) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/(gamma*gamma) + Ay*np.sinh(kappa_pp*xi)
ux = xi*xi*np.cosh(gamma*xi)*F0/(4*np.cosh(gamma)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/(gamma*gamma) + (P_1*kappa_p*kappa_p*np.sinh(kappa)/(gamma*gamma*np.sinh(kappa_p)))*(xi*np.sinh(kappa_p*xi)/2) - P_1*kappa*kappa*xi*np.sinh(kappa*xi)/(2*gamma*gamma) + Ax*np.cosh(kappa_pp*xi)

uy = (-P_1*kappa**2*xi*np.cosh(kappa*xi)/(2*gamma**2) + P_1*kappa*kappa_p*xi*np.sinh(kappa)*np.cosh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/gamma**2) + Ay*np.sinh(kappa_pp*xi)
ux = (F0*xi**2*np.cosh(gamma*xi)/(4*np.cosh(gamma)) - P_1*kappa**2*xi*np.sinh(kappa*xi)/(2*gamma**2) + P_1*kappa_p**2*xi*np.sinh(kappa)*np.sinh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/gamma**2) + Ax*np.cosh(kappa_pp*xi)
for i in range(Nt):
	for x in range(Nx):
		u[0, :, x, i] = np.real(np.exp(1j*omega*T[i])*ux*cos(2*kappa*eta[x]))
		u[1, :, x, i] = np.real(np.exp(1j*omega*T[i])*uy*sin(2*kappa*eta[x]))

for i in range(len(T)):
	plt.plot(u[0, :,27,i])
	#plt.plot(u[1, :,27,i])

plt.show()
x2 = eta
y2 = np.linspace(-1-epsilon, 1+epsilon, 110)
X, Y = np.meshgrid(x2, y2)

a = 1
x = np.zeros((len(xi), len(eta)))
y = np.zeros((len(xi), len(eta)))

for i in range(len(eta)):
	y[:, i] = xi*a*(1+epsilon*np.sin(kappa*eta[i]))

for i in range(len(x)):
	x[i, :] = a*eta

for i in range(len(T)):
	plt.clf()
	# Interpolate using three different methods and plot
	ux = griddata( (x.flatten(),  y.flatten()), u[0,:,:,i].flatten(), (X, Y), method='nearest')
	uy = griddata( (x.flatten(),  y.flatten()), u[1,:,:,i].flatten(), (X, Y), method='nearest')

	speed = np.sqrt(np.square(ux) + np.square(uy))

	ux[np.where(np.abs(uy)<1e-5)] = 0
	uy[np.where(np.abs(uy)<1e-5)] = 0

	plt.streamplot(X, Y, ux, uy, density=0.6, color='k')
	CS = plt.contourf(X, Y, speed, levels=np.linspace(0, 3, 40))
	cbar = plt.colorbar(CS)
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	plt.plot(x[0,:], y[0, :], "k", linewidth=3)
	plt.plot(x[0,:], y[-1,:], "k", linewidth=3)
	plt.fill_between(x[0,:], (-1-epsilon)*np.ones(len(x[0,:])), y[0,:], color="k")
	plt.fill_between(x[0,:], ( 1+epsilon)*np.ones(len(x[0,:])), -y[0,:], color="k")
	plt.draw()
	plt.pause(0.5)
	#plt.savefig("figures/streamplot_scaled_epsilon1_nr"+str(i)+".pdf")
	#plt.savefig("figures/test_"+str(i))
plt.show()

