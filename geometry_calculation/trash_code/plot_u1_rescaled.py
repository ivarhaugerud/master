import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
from scipy.interpolate import griddata
import matplotlib

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

N = int(10)
epsilon = 0.3

Sc = 1/50
F0 = 1
pi = np.pi 
omega = 1
kappa = 1
Nt = 40
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(-2*np.pi, 2*np.pi/kappa, 100)
xi = np.linspace(-1, 1, 500)
gamma = np.sqrt(1j*omega/Sc)

u_x = np.zeros((len(T), len(xi), len(eta)), dtype="complex")
u_y = np.zeros((len(T), len(xi), len(eta)), dtype="complex")


kappa_prime = np.sqrt(1j*omega/Sc + kappa*kappa)
P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))

for t in range(len(T)):
	for x in range(len(xi)):
		u_x[t,x,:]    = epsilon*np.sin(kappa*eta)*(  (P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi[x])/np.cosh(kappa) - np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime) - xi[x]*np.sinh(gamma*xi[x])/np.sinh(gamma))  )
		u_x[t,x,:]   += F0*(1-np.cosh(gamma*xi[x])/np.cosh(gamma))/(gamma*gamma)
		u_x[t,x,:]   *= np.exp(1j*omega*T[t])

		u_y[t, x, :] = (np.exp(1j*omega*T[t])*np.cos(kappa*eta)*kappa*P1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_prime*xi[x])/np.sinh(kappa_prime) - np.sinh(kappa*xi[x])/np.sinh(kappa))

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
	# Interpolate using three different methods and plot
	ux = griddata( (x.flatten(),  y.flatten()), np.real(u_x[i,:,:].flatten()), (X, Y), method='nearest')
	uy = griddata( (x.flatten(),  y.flatten()), np.real(u_y[i,:,:].flatten()), (X, Y), method='nearest')

	speed = np.sqrt(np.square(ux) + np.square(uy))

	ux[np.where(np.abs(uy)<1e-5)] = 0
	uy[np.where(np.abs(uy)<1e-5)] = 0

	plt.streamplot(X, Y, ux, uy, density=0.6, color='k')
	CS = plt.contourf(X, Y, speed, levels=np.linspace(0, 3.5*1e-2, 40))
	cbar = plt.colorbar(CS)
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	plt.plot(x[0,:], y[0, :], "k", linewidth=3)
	plt.plot(x[0,:], y[-1,:], "k", linewidth=3)
	plt.draw()
	plt.pause(0.5)
	plt.clf()
	#plt.savefig("figures/streamplot_scaled_epsilon1_nr"+str(i)+".pdf")
	#plt.savefig("figures/test_"+str(i))
	plt.show()




