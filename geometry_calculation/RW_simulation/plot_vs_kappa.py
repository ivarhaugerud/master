import numpy as np 
import matplotlib.pyplot as plt 

T = 5
periods = 4000
datafiles = periods*100
half_way = int(datafiles/2)
t = np.linspace(0, T*periods, datafiles)

kappas = np.array([0.4, 0.7, 1.0, 1.6])
Lxs    = 2*np.pi/kappas
files = ["15_7", "8_97", "6_28", "3_92"]

var    = np.zeros((len(kappas), datafiles))
D_para = np.zeros((len(kappas), datafiles-1))
D      = np.zeros((len(kappas), 2))

for i in range(len(kappas)):
	var[i, :] = np.load("data/Lx"+files[i]+"_var.npy")

	plt.plot(t, var[i, :])
plt.show()
for i in range(len(kappas)):
	D_para[i, :] = var[i, 1:]/t[1:]
	plt.plot(t[1:], D_para[i, :])
plt.show()


for i in range(len(kappas)):
	D[i, 0] = np.mean(D_para[i, half_way:])
	D[i, 1] = np.std(D_para[i, half_way:])

plt.plot(kappas, D[:,0], "o")
plt.show()