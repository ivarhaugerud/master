import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

kappa = np.array([0.25, 0.5, 1.0, 2, 4, 8])
tau = np.array([0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8])
Lx = [25.13, 12.56, 6.283, 3.141, 1.57, 0.785]

base = "data_test/find_resonance_try2/"
D = np.zeros((len(kappa), len(tau)))
difference = np.zeros(np.shape(D))
dt = tau/500
#dt[1:5] = tau[1:5]/500
Ts = np.ones(len(tau), dtype="int")*int(500)
#Ts[1:5] = np.ones(len(tau[1:5]), dtype="int")*int(500)

for i in range(len(kappa)):
	for j in range(len(tau)):
		T = Ts[j]
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau"+ str(round(tau[j], 3)) +"_eps0.5_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt" + str(round(dt[j], 6)) + "/tdata.dat")
		print(np.shape(data), np.shape(D))
		D[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau[j]
		difference[i, j] = abs(D[i, j] - sci.trapz(  np.trim_zeros(data[:, 8])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau[j])/D[i,j]
		#plt.plot(np.trim_zeros(data[:, 0])/tau[j], np.trim_zeros(data[:, 8]))
		#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[j], np.trim_zeros(data[:, 8])[-T:])
		#plt.title(str(kappa[i]) + ","+ str(tau[j]))
		#plt.show()

for i in range(len(tau)):
	plt.plot(kappa, D[:, i], label=str(tau[i]))

plt.xscale("log")
#plt.yscale("log")
plt.legend(loc="best")
plt.show()
