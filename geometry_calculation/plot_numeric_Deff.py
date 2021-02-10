import numpy as np 
import scipy.integrate as sci
import matplotlib.pyplot as plt 

data = np.load("data_test/tdatas_large_run.npy")
epsilon = np.arange(0.10, 0.46, 0.05)
kappa  = np.array([0.2, 0.6 , 1.0, 1.4, 1.7, 2.1])
Lx = 2*np.pi/kappa
T = 3.0
dt = 0.003
datafiles = int(T/dt)

U = np.zeros((len(epsilon), len(kappa)))
D = np.zeros((len(epsilon), len(kappa)))

print(np.shape(data))
for i in range(len(epsilon)):
   for j in range(len(kappa)):
      U[i,j] = sci.trapz(data[i, j, -datafiles:, 4], data[i, j, -datafiles:, 0])/T
      D[i,j] = sci.trapz(data[i, j, -datafiles:, 8], data[i, j, -datafiles:, 0])/T

      plt.plot(data[i, j, :, 0], data[i, j, :, 8])
      #plt.plot(data[i, j, :, 0], data[i, j, :, 4])
      plt.plot(data[i, j, -datafiles:, 0], data[i, j, -datafiles:, 8])
      #plt.plot(data[i, j, :, 0], D[i,j]*np.ones(len(data[i,j,:,0])))
      plt.title(str(epsilon[i]) + " and " + str(kappa[j]))
      plt.show()
   #plt.plot(kappa, D[i, :])
   #plt.show()

for j in range(len(kappa)):
	plt.figure(1)
	plt.plot(epsilon, U[:, j], "-", label=r"$\kappa=$"+str(kappa[j])[:4])

	plt.figure(2)
	plt.plot(epsilon, D[:, j], "-", label=r"$\kappa=$"+str(kappa[j])[:4])

plt.figure(2)
plt.xlabel(r"Boundary amplitude $\epsilon$")
plt.ylabel(r"Effective dispersion $D_\parallel$")
plt.legend(loc="best")

plt.figure(1)
plt.xlabel(r"boundary amplitude $\epsilon$")
plt.ylabel(r"Average velocity squared $U^2$")
plt.legend(loc="best")
plt.show()



for i in range(len(epsilon)):
	plt.figure(1)
	plt.plot(kappa, U[i, :], "-", label=r"$\epsilon=$"+str(epsilon[i])[:4])

	plt.figure(2)
	plt.plot(kappa, D[i, :], "-", label=r"$\epsilon=$"+str(epsilon[i])[:4])

plt.figure(2)
plt.xlabel(r"Wave number $\kappa$")
plt.ylabel(r"Effective dispersion $D_\parallel$")

#omega = 2*np.pi/T 
#Sc = 1.2
#gamma = np.sqrt(1j*omega/Sc)
#plt.plot([2*np.real(gamma), 2*np.real(gamma)], [0.95, 1.07])
plt.legend(loc="best")

plt.figure(1)
plt.xlabel(r"Wave number $\kappa$")
plt.ylabel(r"Average velocity squared $U^2$")
plt.legend(loc="best")
plt.show()