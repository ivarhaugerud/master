import numpy as np 
import matplotlib.pyplot as plt 
#import seaborn as sns
import pandas as pd

#plt.style.use("bmh")
#sns.color_palette("hls", 1)
import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

N = int(100)
epsilon = 0.3

Sc = 1
F0 = 1
pi = np.pi 
omega = 1
kappa = 1
Nt = 15
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(-2*np.pi, 2*np.pi/kappa, 75)
xi = np.linspace(-1, 1, 50)

u_x = np.zeros(len(xi), dtype="complex")
u_y = np.zeros(len(xi), dtype="complex")
flow_field = np.zeros((len(xi), len(eta), len(T), 2), dtype="complex")

a = 1
x = np.zeros((len(xi), len(eta)))
y = np.zeros((len(xi), len(eta)))

for i in range(len(eta)):
	y[:, i] = xi*a*(1+epsilon*np.sin(kappa*eta[i]))

for i in range(len(x)):
	x[i, :] = a*eta

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
		flow_field[j, :, i, 1] = -(np.real(np.cos(kappa*eta)*u_y[j] - Sc*kappa*np.sinh(kappa*xi[j])*(psi1*np.cos(kappa*eta) + phi1*np.sin(kappa*eta))/(1j*omega)))

	u_x = np.zeros(len(xi), dtype="complex")
	u_y = np.zeros(len(xi), dtype="complex")

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

max_vel = np.amax(np.sqrt(np.real( (epsilon*flow_field[:,:,:,0]+u_x0[:, :, :])*(epsilon*flow_field[:,:,:,0]+u_x0[:, :, :]) + flow_field[:,:,:,1]*flow_field[:,:,:,1])))
from scipy.interpolate import griddata
total_flow = np.zeros((len(xi), len(eta), len(T), 2))

for t in range(len(T)):
	for i in range(len(eta)):
		total_flow[:, i, t, 0] = np.real(epsilon*flow_field[:,i,t,0]+u_x0[:, i, t])
		total_flow[:, i, t, 1] = np.real(epsilon*flow_field[:,i,t,1])

x2 = eta
y2 = np.linspace(-1-epsilon, 1+epsilon, 110)
print(x2, y2)
X, Y = np.meshgrid(x2, y2)

for i in range(len(T)):
	# Interpolate using three different methods and plot
	print(np.shape(total_flow[:,:, i, 0]), np.shape(total_flow[:,:, i, 1]), np.shape(xi), np.shape(eta), np.shape( total_flow[:,:, i, 0].flatten()))
	ux = griddata( (x.flatten(),  y.flatten()), total_flow[:,:, i, 0].flatten(), (X, Y), method='nearest')
	uy = griddata( (x.flatten(),  y.flatten()), total_flow[:,:, i, 1].flatten(), (X, Y), method='nearest')

	#plt.figure(2)
	speed = np.sqrt(np.square(ux) + np.square(uy)) # / speed.max()
	plt.streamplot(X, Y, ux, uy, density=1.5, color='k')
	CS = plt.contourf(X, Y, speed, 15)
	cbar = plt.colorbar(CS)
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	plt.draw()
	plt.pause(0.5)
	plt.clf()
	#plt.savefig("figures/streamplot_scaled_epsilon1_nr"+str(i)+".pdf")
	#plt.savefig("figures/test_"+str(i))
plt.show()




