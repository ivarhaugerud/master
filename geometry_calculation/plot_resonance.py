import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

kappa = np.array([0.25, 0.5, 1.0, 2, 4, 8])
Lx = 2*np.pi/kappa
tau = np.array([0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8])
data = np.load("data_test/find_resonance_tdata.npy")
D = np.zeros((len(kappa), len(tau)))
difference = np.zeros(np.shape(D))
T = 750

for i in range(len(kappa)):
	for j in range(len(tau)):
		D[i, j] = sci.trapz(  np.trim_zeros(data[j, i, :, 8])[-T:],  np.trim_zeros(data[j, i, :, 0])[-T:] )/tau[j]
		difference[i, j] = abs(D[i, j] - sci.trapz(  np.trim_zeros(data[j, i, :, 8])[-2*T:-T],  np.trim_zeros(data[j, i, :, 0])[-2*T:-T] )/tau[j])/D[i,j]
		plt.plot(np.trim_zeros(data[j, i,   :, 0])/tau[j], np.trim_zeros(data[j, i,   :, 8]))
		plt.plot(np.trim_zeros(data[j, i, :, 0])[-T:]/tau[j], np.trim_zeros(data[j, i, :, 8])[-T:])
	plt.show()

for i in range(len(tau)):
	plt.plot(kappa, difference[:, i], "o", label=str(tau[i]))
plt.yscale("log")
plt.legend(loc="best")
plt.show()
