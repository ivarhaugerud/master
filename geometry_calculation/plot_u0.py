import numpy as np 
import matplotlib.pyplot as plt 


N = int(1e2)
xi = np.linspace(-1, 1, 1e4)
u_x = np.zeros(len(xi))

Sc = 0.2
F0 = 1
T = np.linspace(0, 0.2*np.pi, 10)
pi = np.pi 
omega = 1

for i in range(len(T)):
	t = T[i]
	for n in range(N):
		n = 2*n +1
		alpha_n = n*n*pi*pi*Sc/4
		u_x += (-1)**(0.5*(n-1))*np.cos(n*np.pi*xi/2)*(alpha_n*np.cos(omega*t) + omega*np.sin(omega*t))/(n*(alpha_n**2 + omega**2))

	u_x *= 4/np.pi
	plt.plot(xi, u_x, label=r"$t=$"+str(t))

	u_x = np.zeros(len(xi))

plt.show()