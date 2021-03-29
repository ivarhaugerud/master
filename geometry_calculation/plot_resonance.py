import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

tau = np.array([2.0, 2.6, 3.69])
Lx = np.array([4.487, 4.654, 4.833, 5.026, 5.235, 5.463, 5.711, 5.983, 6.283, 6.613, 6.981, 7.391, 7.853, 8.377, 8.975])
kappa = 2*np.pi/Lx 

base = "data_test/find_resonance_try3/"
D = np.zeros((len(kappa), len(tau)))
difference = np.zeros(np.shape(D))
dt = tau/500
Ts = np.ones(len(tau), dtype="int")*int(500)

for i in range(len(kappa)):
	for j in range(len(tau)):
		T = Ts[j]
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau"+ str(round(tau[j], 3)) +"_eps0.5_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt" + str(round(dt[j], 6)) + "/tdata.dat")
		print(kappa[i], tau[j], np.shape(data), np.shape(D))
		D[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau[j]
		difference[i, j] = abs(D[i, j] - sci.trapz(  np.trim_zeros(data[:, 8])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau[j])/D[i,j]

		if D[i, j] > 2.0:
			D[i,j] = 0
			#plt.plot(np.trim_zeros(data[:, 0])/tau[j], np.trim_zeros(data[:, 8]))
			#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[j], np.trim_zeros(data[:, 8])[-T:])
			#plt.title(str(kappa[i]) + ","+ str(tau[j]))
			#plt.show()

for i in range(len(tau)):
	plt.plot(kappa, D[:, i], "o", label=str(tau[i]))

#plt.xscale("log")
#plt.yscale("log")
plt.legend(loc="best")
plt.show()
