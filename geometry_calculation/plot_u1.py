import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd

#plt.style.use("bmh")
#sns.color_palette("hls", 1)
import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

N = int(500)
epsilon = 0.10

Sc = 10
F0 = 1
pi = np.pi 
omega = 1
kappa = 1
Nt = 10
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(0, 2*np.pi/kappa, 100)
xi = np.linspace(-1, 1, 100)

u_x = np.zeros(len(xi), dtype="complex")
u_y = np.zeros(len(xi), dtype="complex")
flow_field = np.zeros((len(xi), len(eta), len(T), 2), dtype="complex")


for i in range(len(T)):
	t = T[i]
	psi1 = 0
	phi1 = 0

	for n in range(N):
		n = 2*n + 1
		alpha_n = n*n*pi*pi*Sc/4

		an = 4*Sc*F0*((-1)**(0.5*(n-1)))*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n + omega*omega))
		bn = -2*kappa*kappa*F0*((-1)**(0.5*(n-1)))/(kappa*kappa + n*n*pi*pi/4)

		dn = bn - an*n*pi*kappa*kappa/2 
		cn = an*n*n*pi*pi/2 + bn*n*pi/(kappa*kappa + n*n*pi*pi/4)

		theta_n = -2*kappa*F0*((-1)**(0.5*(n-1)))/(kappa*kappa + n*n*pi*pi/4)
		denominator = 1j*omega/Sc + kappa*kappa + n*n*pi*pi/4
		lambda_n = (1-n*n*pi*pi*0.5*(1/denominator + 1/(kappa*kappa+n*n*pi*pi/4)))

		psi1 += ((-1)**(0.5*(n-1)))*(theta_n*lambda_n*np.sinh(kappa) - dn*np.cosh(kappa))/denominator
		phi1 += ((-1)**(0.5*(n-1)))*(theta_n*lambda_n*np.cosh(kappa) + dn*np.sinh(kappa))/denominator

		u_x += cn*np.cos(n*pi*xi/2)/(1j*omega/Sc + kappa*kappa + n*n*pi*pi/4)
		u_x += dn*(xi*np.sin(n*pi*xi/2)/denominator + n*pi*np.cos(n*pi*xi/2)/(denominator*denominator))


		u_y += theta_n*np.sin(n*pi*xi/2)*(1-n*n*pi*pi*0.5/(kappa*kappa+n*n*pi*pi/4))/denominator
		u_y += theta_n*n*pi*0.5*(xi*np.cos(n*pi*xi/2)/denominator - n*pi*np.sin(n*pi*xi/2)/(denominator*denominator))
	
	psi1 *= 1j*omega/(Sc*kappa*np.cosh(2*kappa))*np.exp(1j*omega*t)
	phi1 *= 1j*omega/(Sc*kappa*np.cosh(2*kappa))*np.exp(1j*omega*t)

	u_x *= np.exp(1j*omega*t)
	u_y *= np.exp(1j*omega*t)

	for j in range(len(xi)):
		flow_field[j, :, i, 0] = np.real(np.sin(kappa*eta)*u_x[j] + Sc*kappa*np.cosh(kappa*xi[j])*(psi1*np.sin(kappa*eta) - phi1*np.cos(kappa*eta))/(1j*omega))
		flow_field[j, :, i, 1] = np.real(np.cos(kappa*eta)*u_y[j] - Sc*kappa*np.sinh(kappa*xi[j])*(psi1*np.cos(kappa*eta) + phi1*np.sin(kappa*eta))/(1j*omega))

	u_x = np.zeros(len(xi), dtype="complex")
	u_y = np.zeros(len(xi), dtype="complex")
"""
for t in range(len(T)):
	X, Y = np.meshgrid(eta, xi)
	plt.quiver(X, Y, flow_field[:,:,t,0], flow_field[:,:,t,1])
	plt.draw()
	plt.pause(0.01)
	plt.clf()
plt.show()
"""
u_x = np.zeros(len(xi), dtype=complex)
u_x0 = np.zeros((len(xi), len(eta), len(T)), dtype=complex)

for i in range(len(T)):
	t = T[i]
	for n in range(N):
		n = 2*n +1
		alpha_n = n*n*pi*pi*Sc/4
		u_x += (-1)**(0.5*(n-1))*np.cos(n*np.pi*xi/2)*(alpha_n*np.cos(omega*t) + omega*np.sin(omega*t))/(n*(alpha_n**2 + omega**2))

	u_x *= 4/np.pi

	for j in range(len(xi)):
		u_x0[j, :, i] = u_x[j]

	u_x = np.zeros(len(xi))
"""
for t in range(len(T)):
	X, Y = np.meshgrid(eta, xi)
	Z = np.sqrt((epsilon*flow_field[:,:,t,0]+u_x0[:, :, t])*(epsilon*flow_field[:,:,t,0]+u_x0[:, :, t]) + epsilon*flow_field[:,:,t,1]*epsilon*flow_field[:,:,t,1])

	cp = plt.contourf(X, Y, Z, levels=np.linspace(0, 0.55, 30))
	cb = plt.colorbar(cp)

	plt.quiver(X, Y, epsilon*flow_field[:,:,t,0]+u_x0[:, :, t], epsilon*flow_field[:,:,t,1], linewidths=Z.flatten()/60)
	plt.draw()
	#plt.axis([min(eta), max(eta), -1.01, 1.1])
	plt.pause(0.01)
	plt.clf()
plt.show()
"""
max_vel = np.amax(epsilon*flow_field[:,:,:,0]+u_x0[:, :, :])

for i in range(len(T)):
	plt.clf()
	plt.title(r"$T=$"+str(T[i])[:4] + r"$[\omega^{-1}]$")
	X, Y = np.meshgrid(eta, xi)
	Z = np.sqrt((epsilon*flow_field[:,:,i,0]+u_x0[:, :, i])*(epsilon*flow_field[:,:,i,0]+u_x0[:, :, i]) + epsilon*flow_field[:,:,i,1]*epsilon*flow_field[:,:,i,1])

	cp = plt.contourf(X, Y, Z)#, levels=np.linspace(0, 0.6, 40))
	cb = plt.colorbar(cp)

	lw = 2*Z/max_vel
	plt.streamplot(X, Y, np.real(epsilon*flow_field[:,:,i,0]+u_x0[:, :, i]), np.real(epsilon*flow_field[:,:,i,1]), density=0.5, color='k', linewidth=lw)

	plt.axis([min(eta), max(eta), -1.01, 1.01])

	plt.draw()
	plt.pause(0.2)
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	plt.savefig("figures/streamplot_epsilon1_nr"+str(i)+".pdf")
plt.show()
