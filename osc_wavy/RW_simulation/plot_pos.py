import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

dt = 0.006
tau = 3.0 
timesteps = int(tau/(dt))
periods   = 6000
datafiles = 25000
skip = int(periods*timesteps/datafiles)
t = np.linspace(0, periods*tau, datafiles)

#geometry parameters
visc = np.array([1.5, 3.0, 5.0])
var  = np.zeros(datafiles)
mean = np.zeros(datafiles)
dirr = []

U = np.zeros(len(visc))

for i in range(len(visc)):
    dirr.append("flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006/")
    tdat = np.loadtxt(dirr[i] +"/tdata.dat")
    U[i]  = np.sqrt(integrate.trapz(tdat[-timesteps:, 4], tdat[-timesteps:, 0])/tau)


Pe = 10
D  = U/Pe

for l in range(len(visc)):
	for i in range(datafiles):
		pos = np.load(dirr[l]+"pos/RW_positions_"+str(int(i*skip))+".npy")
		var[i]  = np.var(pos[0, :])
		mean[i] = np.mean(pos[0, :])

	var /= (2*D[l]*t)
	plt.plot(var)
	np.save(dirr[l]+"pos/var",   var)
	np.save(dirr[l]+"pos/mean", mean)
	print("saved ", visc[l])
plt.show()