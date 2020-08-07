import numpy as np 
import matplotlib.pyplot as plt 

N = int(1e4)
kappa = 1
epsilon = 0.2
eta = np.linspace(-15, 15, N)
y = np.linspace(-1, 1, N)
xi = 1/(1+epsilon*np.sin(kappa*eta))

plt.plot(eta, xi)
plt.plot(eta, -xi)
plt.show()