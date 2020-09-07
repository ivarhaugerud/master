import numpy as np 
import matplotlib.pyplot as plt 

xi = np.linspace(-1, 1, int(1e6))
N = int(1e4)
sol = np.zeros(len(xi))
pi = np.pi 

for i in range(N):
	n = 2*i+1
	sol += xi*((-1)**(0.5*(n-1)))*np.sin(n*pi*xi/2)

plt.plot(xi, sol)
plt.show()
