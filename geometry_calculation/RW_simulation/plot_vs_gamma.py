import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py

#simulation parameters
dt           = 0.006
tau          = 3.0 
T            = tau
timesteps    = int(tau/(dt))
periods      = 5000
datafiles    = periods*20
N            = int(1e3)
Pe           = 1

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
var    = np.zeros((len(visc), datafiles))
dirr = []

for i in range(len(visc)):
    dirr.append("flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006")


for i in range(len(visc)):
	Sc = visc[i]#/Dm
	F0 = 12/visc[i]
	gamma   = np.sqrt(1j*omega/Sc)
	gamma_c = np.conj(gamma)
	a       = np.sqrt(1j*omega)
	a_c     = np.conj(a)
	rho = a 
	rho_c = a_c

	factor  = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
	D_ana[i] = (1 + np.real(factor * 0.5 * integrate.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)))
	

for i in range(len(visc)):
    tdat = np.loadtxt(dirr[i] +"/tdata.dat")
    time = tdat[:,0]
    u2   = tdat[:,4]
    exp_u2[i] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)
    exp_D[i] = integrate.trapz(tdat[-timesteps:, 8], time[-timesteps:])/(tau)
    Dm = np.sqrt(exp_u2[i])/Pe
    var[i, :] = np.load(dirr[i] +"/var.npy")/Dm
    #plt.plot(np.loadtxt("flow_fields/zero_eps/mu_"+str(nus[i]) +"/tdata.dat")[:, 8])
    num_D_para[i] = integrate.trapz(np.loadtxt(dirr[i] +"/tdata.dat")[-3000:, 8], np.loadtxt(dirr[i]+"/tdata.dat")[-3000:, 0])/T

plt.plot(visc, num_D_para, "o", label="numerisk")
plt.plot(visc, D_ana, "o", label="analytisk")
plt.legend(loc="best")
plt.xlabel("viskositet")
plt.ylabel("D_eff")
plt.show()

for i in range(len(nus)):
	D_para[i, :] = var[i, 1:]/t[1:]
	plt.plot(t[1:], D_para[i, :], label=str(nus[i]))
	plt.plot(t, D_ana[i]*np.ones(len(t)), label=str(nus[i]))
plt.legend(loc="best")
plt.show()


for i in range(len(nus)):
	D[i, 0] = np.mean(D_para[i, half_way:])
	D[i, 1] = np.std(D_para[i, half_way:])

D[0, 0] = np.mean(D_para[0, int(half_way):])
D[0, 1] =  np.std(D_para[0, int(half_way):])

plt.plot(nus, D_ana)
plt.errorbar(nus, D[:,0], yerr=D[:,1], fmt="o")
plt.show()
np.save("data/D_eff_vs_gamma", D)