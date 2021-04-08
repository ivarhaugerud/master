import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
from scipy.interpolate import griddata
import matplotlib
import matplotlib.ticker as tick 
import os

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
root = "../../../master_latex/results/figures/"

epsilon = 0.2

Sc = 0.2
F0 = 12
pi = np.pi 
tau = 3.0
omega = 2*pi/tau
kappa = 1.0
Nt = 10
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(-np.pi/(2*kappa), np.pi/kappa+np.pi/(2*kappa), 100)
xi = np.linspace(-1, 1, 500)
gamma = np.sqrt(1j*omega/Sc)

u_x = np.zeros((len(T), len(xi), len(eta)), dtype="complex")
u_y = np.zeros((len(T), len(xi), len(eta)), dtype="complex")

kappa_prime = np.sqrt(1j*omega/Sc + kappa*kappa)
P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
P_1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))

for t in range(len(T)):
	for x in range(len(xi)):
		u_x[t,x,:]    = epsilon*np.sin(kappa*eta)*(  (P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi[x])/np.cosh(kappa) - np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime) - xi[x]*np.sinh(gamma*xi[x])/np.sinh(gamma))  )
		u_x[t,x,:]   += F0*(1-np.cosh(gamma*xi[x])/np.cosh(gamma))/(gamma*gamma)
		u_x[t,x,:]   *= np.exp(1j*omega*T[t])

		u_y[t, x, :] = epsilon*(np.exp(1j*omega*T[t])*np.cos(kappa*eta)*kappa*P1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_prime*xi[x])/np.sinh(kappa_prime) - np.sinh(kappa*xi[x])/np.sinh(kappa))

u_x = np.real(u_x)
u_y = np.real(u_y)

from cmath import *
Ax = sqrt(gamma**2 + 4*kappa**2)*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(4*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
Ay = kappa*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(2*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
psi_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

uy = (-P_1*kappa**2*xi*np.cosh(kappa*xi)/(2*gamma**2) + P_1*kappa*kappa_p*xi*np.sinh(kappa)*np.cosh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/gamma**2) + Ay*np.sinh(kappa_pp*xi)
ux = F0*xi**2*np.cosh(gamma*xi)/(4*np.cosh(gamma)) - P_1*kappa**2*xi*np.sinh(kappa*xi)/(2*gamma**2) + P_1*kappa_p*kappa_p*np.sinh(kappa)*xi*np.sinh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/gamma**2 + Ax*np.cosh(kappa_pp*xi)
ux_no_eta = P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma)
U2 = np.zeros(len(T))

for i in range(len(T)):
	for x in range(len(eta)):
		u_x[i, :, x] += epsilon*epsilon*np.real(np.exp(1j*omega*T[i])*ux*np.cos(2*kappa*eta[x])) + epsilon*epsilon*np.real(np.exp(1j*omega*T[i])*ux_no_eta)
		u_y[i, :, x] += epsilon*epsilon*np.real(np.exp(1j*omega*T[i])*uy*np.sin(2*kappa*eta[x]))
	U2[i] = np.max(u_x)

U = np.max(U2)

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

tot_height = max(eta)-(min(eta))
print(tot_height)

D = tot_height/2 - 1
print(D)

stream_points   = np.array(list(zip( pi/(2*kappa)*np.ones(20), np.linspace(-1-epsilon, 1+epsilon, 20))))
stream_points2  = np.array(list(zip( pi/(2*kappa)*np.ones(20)+2*np.pi/kappa, np.linspace(-1-epsilon, 1+epsilon, 20))))

#fig = plt.figure(1)
for i in range(len(T)):
	#ax = fig.add_subplot(111)
	#ax.plot([min(eta),max(eta)],[-1-epsilon-0.1,1+epsilon+0.1])
	#ax.set_aspect(1)
	#ax.set_xlim(min(eta), max(eta))
	#fig = plt.figure(1, figsize=(16.8, 2.95), edgecolor="white")
	plt.clf()
	Map = matplotlib.cm.get_cmap('Spectral_r')
	# Interpolate using three different methods and plot
	ux = griddata( (x.flatten(),  y.flatten()), u_x[i,:,:].flatten()/U, (X, Y), method='nearest')
	uy = griddata( (x.flatten(),  y.flatten()), u_y[i,:,:].flatten()/U, (X, Y), method='nearest')
	speed = np.sqrt(np.square(ux) + np.square(uy))
	ux[np.where(np.abs(uy)<1e-5)] = 0
	uy[np.where(np.abs(uy)<1e-5)] = 0

	plt.streamplot(X, Y, ux, uy, color='k', start_points=stream_points, density=35)
	CS = plt.contourf(X, Y, speed, levels=np.linspace(0, 1.0, 15), cmap=Map)
	plt.streamplot(X+2*np.pi/kappa, Y, ux, uy, color='k', start_points=stream_points2, density=35)
	CS = plt.contourf(X+2*np.pi/kappa, Y, speed, levels=np.linspace(0, 1.0, 15), cmap=Map)
	cbar = plt.colorbar(CS)
	cbar.set_label(r"Velocity $[U_{max}]$", fontsize=12)
	cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter("%.1f"))
	plt.xlabel(r"Horizontal position $x$   $[a]$  ", fontsize=12)
	plt.ylabel(r"Vertical position   $y$   $[a]$  ", fontsize=12)
	plt.plot(x[0,:], y[0, :], "k", linewidth=3)
	plt.plot(x[0,:], y[-1,:], "k", linewidth=3)
	plt.fill_between(x[0,:], (-1.25)*np.ones(len(x[0,:])), y[0,:], color="k")
	plt.fill_between(x[0,:], ( 1.25)*np.ones(len(x[0,:])), -y[0,:], color="k")
	plt.fill_between(x[0,:]+2*np.pi/kappa, (-1.25)*np.ones(len(x[0,:])), y[0,:], color="k")
	plt.fill_between(x[0,:]+2*np.pi/kappa, ( 1.25)*np.ones(len(x[0,:])), -y[0,:], color="k")
	#plt.axis([min(eta), max(eta), -1-epsilon-0.05, 1+epsilon+0.05])
	plt.axis([min(eta), max(eta)+2*np.pi/kappa, -1-epsilon-0.05, 1+epsilon+0.05])
	#plt.axis("equal")
	plt.draw()
	#fig.tight_layout()
	#plt.axis("scaled")
	plt.pause(0.5)

	filename = root + "streamplot_equal_nr"+str(i)+".pdf"
	plt.savefig(filename, bbox_inches="tight")
	os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
