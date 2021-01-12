import numpy as np 
import matplotlib.pyplot as plt 

T = 5
periods = 4000
datafiles = periods*100
half_way = int(datafiles/2)
t = np.linspace(0, T*periods, datafiles)

kappas = np.array([0.1, 0.4, 0.7, 1.0, 1.3, 1.6])
Lxs    = 2*np.pi/kappas
files = ["62_8", "15_7", "8_97", "6_28", "4_83", "3_92"]

var    = np.zeros((len(kappas), datafiles))
D_para = np.zeros((len(kappas), datafiles-1))
D      = np.zeros((len(kappas), 2))

for i in range(len(kappas)):
	var[i, :] = np.load("data/Lx"+files[i]+"_var.npy")

var[0, :] /= 4

for i in range(len(kappas)):
	plt.plot(t, var[i, :], label=files[i])

plt.legend(loc="best")	
plt.show()
for i in range(len(kappas)):
	D_para[i, :] = var[i, 1:]/t[1:]
	plt.plot(t[1:], D_para[i, :], label=files[i])

plt.legend(loc="best")
plt.show()


for i in range(len(kappas)):
	D[i, 0] = np.mean(D_para[i, half_way:])
	D[i, 1] = np.std(D_para[i, half_way:])

D[0, 0] = np.mean(D_para[0, int(3*half_way/2):])
D[0, 1] =  np.std(D_para[0, int(3*half_way/2):])

plt.errorbar(kappas, D[:,0], yerr=D[:,1], fmt="o")
plt.show()
np.save("data/D_eff_vs_kappa", D)