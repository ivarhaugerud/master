import numpy as np
import matplotlib.pyplot as plt

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 20000
datafiles = periods*10
skip = int(periods*timesteps/datafiles)

#geometry parameters
visc = np.array([1.5, 3.0, 5.0])

dirr = []

for i in range(len(visc)):
    dirr.append("flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

var  = np.zeros(datafiles)
mean = np.zeros(datafiles)

N = 1000
x = np.zeros(N)
y = np.zeros(N)

print(np.sqrt(2*0.006/10))
for l in range(len(visc)):
	for i in range(datafiles):
		plt.clf()
		pos = np.load(dirr[l]+"pos2/RW_positions_"+str(int(i*skip))+".npy")
		#var[i]  = np.square(np.std(pos[0, :]))
		mean[i] = np.mean(pos[0, :])
		plt.scatter(pos[0, :], pos[1, :])
		plt.pause(0.01)
		#x[i] = pos[0, 5]
		#y[i] = pos[1, 5]

	plt.plot(x, y)
	plt.show()

	plt.plot(mean)
	plt.show()
	#np.save(dirr[l]+"var",   var)
	#np.save(dirr[l]+"mean", mean)
	#print("saved ", Lx[l])
