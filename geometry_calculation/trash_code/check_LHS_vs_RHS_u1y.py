import numpy as np 
import matplotlib.pyplot as plt

N = int(2000)
epsilon = 0.10

Sc = 1
F0 = 1
pi = np.pi 
omega = 1
kappa = 1

T = np.linspace(0, 4*np.pi/omega, 8)
eta = np.linspace(0, 4*np.pi/kappa, 500)
xi = np.linspace(-1, 1, 1e4)

u_x = np.zeros(len(xi))
u_y = np.zeros(len(xi))
LHS = np.zeros((len(xi), len(eta), len(T)), dtype="complex")
RHS = np.zeros((len(xi), len(eta), len(T)), dtype="complex")

temp_RHS = np.zeros(len(xi), dtype="complex")
temp_LHS = np.zeros(len(xi), dtype="complex")


for i in range(len(T)):
	t = T[i]

	for n in range(N):
		n = 2*n + 1

		term1 = (0.5 - n*n*pi*pi*0.25/(kappa*kappa + n*n*pi*pi/4))  *  (1j*n*pi/2*np.exp(1j*n*pi*xi/2)/(n*pi/2 - omega/Sc) - 1j*n*pi/2*np.exp(-1j*n*pi*xi/2)/(n*pi/2 + omega/Sc))
		term2 = 1j*n*pi/4  *  (np.exp(1j*n*pi*xi/2)/(n*pi/2 - omega/Sc ) - np.exp(-1j*n*pi*xi/2)/(n*pi/2 + omega/Sc ) + 0.5*1j*n*pi*xi*np.exp(1j*n*pi*xi/2)/(n*pi/2 - omega/Sc ) + 0.5*1j*n*pi*xi*np.exp(-1j*n*pi*xi/2)/(n*pi/2 + omega/Sc ) )
		term3 = 0.25*n*pi  *  ( 0.5*1j*n*pi*np.exp(-1j*n*pi*xi/2)/(n*pi/2 + omega/Sc)**2 - (0.5*1j*n*pi*np.exp(1j*n*pi*xi/2))/(n*pi/2 - omega/Sc)**2 )
		TERM2 = 2*kappa*F0*np.exp(1j*omega*t)*(-1)**(0.5*(n-1))*(term1 + term2 + term3)/(kappa*kappa + n*n*pi*pi/4)

		term1 = (0.5 - n*n*pi*pi*0.25/(kappa*kappa + n*n*pi*pi/4)) * (np.exp(-1j*n*pi*xi/2)/(n*pi/2 + omega/Sc) + np.exp(1j*n*pi*xi/2)/(n*pi/2 - omega/Sc))
		term2 = 0.25*n*pi*1j*xi*(np.exp(1j*n*pi*xi/2)/(n*pi/2 - omega/Sc) - np.exp(-1j*n*pi*xi/2)/(n*pi/2 + omega/Sc))
		term3 = -0.25*n*pi* (-np.exp(1j*n*pi*xi/2)/(n*pi/2-omega/Sc)**2 + np.exp(-1j*n*pi*xi/2)/(n*pi/2 + omega/Sc)**2) 
		TERM1 = -2*1j*omega*kappa*F0*np.exp(1j*omega*t)*(-1)**(0.5*(n-1))*(term1 + term2 + term3)/(Sc*kappa*kappa + Sc*n*n*pi*pi/4)
		
		temp_RHS += TERM1 + TERM2

		temp_LHS += (np.sin(n*pi*xi/2) + 0.5*n*pi*xi*np.cos(n*pi*xi/2)-0.5*n*n*pi*pi*np.sin(n*pi*xi/2)/(kappa*kappa + n*n*pi*pi/4))*((-1)**(0.5*(n-1)))/(kappa*kappa + n*n*pi*pi/4) + 0*1j

	temp_LHS *= 2*kappa*F0*np.exp(1j*omega*t)

	#plt.plot(xi, np.real(temp_RHS))
	#plt.plot(xi, np.real(temp_LHS), "--")

	#plt.plot(xi, np.imag(temp_RHS))
	#plt.plot(xi, np.imag(temp_LHS), "--")

	plt.plot(xi, np.real(temp_RHS)/np.real(temp_LHS))

	plt.show()
	temp_RHS = np.zeros(len(xi), dtype="complex")
	temp_LHS = np.zeros(len(xi), dtype="complex")

	#for j in range(len(xi)):
		#LHS[j, :, i] = np.cos(kappa*eta)*TEMP_LHS
		#RHS[j, :, i] = np.cos(kappa*eta)*temp_RHS