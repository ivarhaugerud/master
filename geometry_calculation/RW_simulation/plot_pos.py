import numpy as np
import matplotlib.pyplot as plt

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 5000
datafiles = periods*20
period    = tau
datafiles = periods*20
skip = int(periods*timesteps/datafiles)

#geometry parameters
Lx =  np.array([1.05, 2.09, 6.28, 9.42, 12.56, 15.71, 25.13])
visc = np.array([1.5, 3.0, 5.0])

dirr = []

for i in range(len(visc)):
    dirr.append("flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

var  = np.zeros(datafiles)
mean = np.zeros(datafiles)

x = np.zeros(datafiles)
y = np.zeros(datafiles)

datafiles = 5000

x = np.zeros(datafiles)
y = np.zeros(datafiles)
prev_pos = np.zeros((2, 1000))
for l in range(len(visc)):
	for i in range(10):
		plt.clf()
		pos = np.load(dirr[l]+"pos/RW_positions_"+str(int(i*skip))+".npy")
		var[i]  = np.square(np.std(pos[0, :]))
		mean[i] = np.mean(pos[0, :])
		x[i] = np.mean(pos[0, :]-prev_pos[0, :])
		#y[i] = np.mean(pos[1, :]-prev_pos[1, :])
		prev_pos = np.copy(pos)
		plt.scatter(pos[0,:], pos[1, :])
		#plt.axis([-5, 5, -1, 1])
		plt.axis("equal")
		plt.title(str(dt*i*skip/tau))
		plt.pause(0.51)
	plt.plot(x)
	plt.show()

	#np.save(dirr[l]+"var",   var)
	#np.save(dirr[l]+"mean", mean)
	#print("saved ", Lx[l])