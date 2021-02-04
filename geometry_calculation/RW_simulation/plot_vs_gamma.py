import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py

Pe           = 1
dt           = 0.006
tau          = 3.0 
timesteps    = int(tau/(dt))
#geometry parameters
epsilon = 0.0
visc = np.array([1.5, 3.0, 5.0])
Lx = 12.56
omega = 2*np.pi/tau
xi = np.linspace(-1, 1, int(1e5))

#flow parameters
periods_of_flow = 8/3
exp_u2 = np.zeros(len(visc))
exp_D  = np.zeros(len(visc))
D_ana  = np.zeros(len(visc))
num_D_para  = np.zeros(len(visc))
RW_D_para   = np.zeros((2, len(visc)))
dirr = []

for i in range(len(visc)):
    dirr.append("flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006")


for i in range(len(visc)):
    Sc = visc[i]
    F0 = 12/visc[i]
    tdat = np.loadtxt(dirr[i] +"/tdata.dat")
    time = tdat[:,0]
    u2   = tdat[:,4]
    exp_u2[i] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)
    gamma   = np.sqrt(1j*omega/Sc)
    gamma_c = np.conj(gamma)
    a       = np.sqrt(1j*omega)
    a_c     = np.conj(a)
    rho = a 
    rho_c = a_c
    Pe = 1 
    num_D_para[i] = integrate.trapz(np.loadtxt(dirr[i] +"/tdata.dat")[-timesteps:, 8], np.loadtxt(dirr[i]+"/tdata.dat")[-timesteps:, 0])/tau
    D_ana[i] = 1 + np.real(Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*omega*omega*omega*(Sc*Sc-1))*(1/(gamma*np.tanh(gamma)) + 1j/(gamma*np.tanh(gamma_c)) - 1/(rho*np.tanh(rho)) - 1j/(rho*np.tanh(rho_c))))
    print(D_ana[i], num_D_para[i])

plt.figure(4)
plt.plot(visc, num_D_para, "ro", label="numerisk brenner")
plt.plot(visc, D_ana, "ko", label="analytisk")
plt.plot(visc, num_D_para, "r-")
plt.plot(visc, D_ana, "k-")
#plt.show()


#simulation parameters
periods      = 20000
datafiles    = periods*10
t = np.linspace(0, tau*periods, datafiles)
half_way     = int(2*datafiles/4)
var    = np.zeros((len(visc), datafiles))


for l in range(len(visc)):
	var = np.load(dirr[l] + "/var.npy")
	mean = np.load(dirr[l]+ "/mean.npy")
	Dm = 1#np.sqrt(exp_u2[l])/Pe
	plt.figure(2)
	plt.plot(t/tau, var/Dm, label="visc="+str(visc[l]))

	plt.figure(3)
	plt.plot(t[1:]/tau, var[1:]/(Dm*t[1:]), label="visc="+str(visc[l]))
	plt.plot(t[half_way:]/tau, var[half_way:]/(Dm*t[half_way:]))

	RW_D_para[0,l] = np.mean(var[half_way:]/(2*Dm*t[half_way:]))
	RW_D_para[1,l] = np.std( var[half_way:]/(2*Dm*t[half_way:]))

plt.legend(loc="best")

plt.figure(3)
plt.xlabel("time [period]")
plt.ylabel("MSD/t")

plt.figure(2)
plt.xlabel("time [period]")
plt.ylabel("MSD")

plt.figure(4)
plt.errorbar(visc, (RW_D_para[0, :]), yerr=RW_D_para[1, :], label="RW")
plt.legend(loc="best")
plt.xlabel("viskositet")
plt.ylabel("D_eff")
plt.show()