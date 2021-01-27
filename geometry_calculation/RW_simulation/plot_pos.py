import numpy as np
import matplotlib.pyplot as plt

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 10000
datafiles = periods*20
period    = tau
periods = 5000
datafiles = periods*20
skip = int(periods*timesteps/datafiles)

#geometry parameters
Lx =  np.array([1.05, 2.09, 6.28, 9.42, 12.56, 15.71, 25.13])
dirr = []

for i in range(len(Lx)):
    dirr.append("flow_fields/non_zero_eps/Lx"+str(Lx[i])+"_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

var  = np.zeros(datafiles)
mean = np.zeros(datafiles)

t = np.linspace(0, period*periods, datafiles)
x = np.zeros(datafiles)
y = np.zeros(datafiles)

for l in range(len(Lx)):
	for i in range(datafiles):
		pos = np.load(dirr[l]+"pos/RW_positions_"+str(int(i*skip))+".npy")
		var[i]  = np.square(np.std(pos[0, :]))
		mean[i] = np.mean(pos[0, :])

	np.save(dirr[l]+"var2",   var)
	np.save(dirr[l]+"mean2", mean)
	print("saved ", Lx[l])