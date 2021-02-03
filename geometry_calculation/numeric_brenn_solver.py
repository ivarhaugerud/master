import numpy as np 
import matplotlib.pyplot as plt 
from scipy import integrate

data1 = np.load("data_test/tdatas_zeroeps_res200.npy")
data2 = np.load("data_test/tdatas_zeroeps_res300.npy")

visc = np.array([1.5, 3.0, 5.0])
data1_B2 = np.zeros(len(visc))
data2_B2 = np.zeros(len(visc))

tau = 3
dt = 0.006
timesteps = int(tau/dt)
print(np.shape(data1))
for i in range(len(visc)):
	data1_B2[i] = integrate.trapz(data1[i, -timesteps:, 8], data1[i, -timesteps:, 0])/tau
	data2_B2[i] = integrate.trapz(data2[i, -timesteps:, 8], data2[i, -timesteps:, 0])/tau

plt.plot(visc, data1_B2)
plt.plot(visc, data2_B2)
plt.show()
