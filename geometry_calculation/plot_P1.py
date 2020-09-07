import numpy as np 
import matplotlib.pyplot as plt 

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

N = int(1000)
xi = np.linspace(-1, 1, 100)
P1 = np.zeros(len(xi), dtype="complex")
particular = np.zeros(len(xi), dtype="complex")
term1 = np.zeros(len(xi), dtype="complex")
term2 = np.zeros(len(xi), dtype="complex")

Sc = 1
F0 = 1
kappa = 1
pi = np.pi 
omega = 1
Nt = 10
epsilon = 0.3
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(0, 2*np.pi/kappa, 90)
full_solution = np.zeros((len(xi), len(eta), len(T)))

for i in range(len(T)):
	t = T[i]
	for n in range(N):
		n = 2*n +1
		alpha_n = n*n*pi*pi*Sc/4

		an = 4*Sc*F0*((-1)**(0.5*(n-1)))*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n + omega*omega))
		bn = -2*kappa*kappa*F0*((-1)**(0.5*(n-1)))/(kappa*kappa + n*n*pi*pi/4)

		d_n = bn - an*n*pi*kappa*kappa/2 
		c_n = an*n*n*pi*pi/2 + bn*n*pi/(kappa*kappa + n*n*pi*pi/4)

		theta_n = -2*kappa*F0*((-1)**(0.5*(n-1)))/(kappa*kappa + n*n*pi*pi/4)
		denominator = 1j*omega/Sc + kappa*kappa + n*n*pi*pi/4
		lambda_n = (1-n*n*pi*pi*0.5*(1/denominator + 1/(kappa*kappa+n*n*pi*pi/4)))


		particular += -2*kappa*F0*(-1)**(0.5*(n-1))*( xi*np.sin(n*pi*xi/2) + (n*pi*np.cos(n*pi*xi/2))/(kappa*kappa + n*n*pi*pi/4))/(kappa*kappa + n*n*pi*pi/4)
		#term1 += (-1)**(0.5*(n-1))*(theta_n*lambda_n*np.sinh(kappa) - d_n*np.cosh(kappa))/(1j*omega/Sc + kappa*kappa + n*n*pi*pi/4)
		#	term2 += (-1)**(0.5*(n-1))*(theta_n*lambda_n*np.cosh(kappa) + d_n*np.sinh(kappa))/(1j*omega/Sc + kappa*kappa + n*n*pi*pi/4)

	term1 *= np.exp(1j*omega*t)*(1j*omega*np.cosh(kappa*xi)/(Sc*kappa*np.cosh(2*kappa)))
	term2 *= np.exp(1j*omega*t)*(1j*omega*np.cosh(kappa*xi)/(Sc*kappa*np.cosh(2*kappa)))
	particular *= np.exp(1j*omega*t)

	for j in range(len(xi)):
		full_solution[j, :, i] = np.real(particular[j]*np.cos(kappa*eta) + term1[j]*np.cos(kappa*eta) + term2[j]*np.sin(kappa*eta))

	particular = np.zeros(len(xi), dtype="complex")
	term1 = np.zeros(len(xi), dtype="complex")
	term2 = np.zeros(len(xi), dtype="complex")

for i in range(len(T)):
	plt.title(r"$T=$"+str(T[i])[:4] + r"$[\omega^{-1}]$", fontsize=12)
	X, Y = np.meshgrid(eta, xi)
	Z = epsilon*full_solution[:,:,i]

	cp = plt.contourf(X, Y, Z, 40)#levels=np.linspace(-0.2, 0.2, 40))
	cb = plt.colorbar(cp)
	plt.draw()
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	plt.savefig("figures/pressure_epsilon1_nr"+str(i)+".pdf")

	plt.pause(0.2)
	plt.clf()

plt.show()
