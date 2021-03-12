import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py

Pe           = 10
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

from numpy import *
for i in range(len(visc)):
    Sc = visc[i]
    F0 = 12/visc[i]
    tdat = np.loadtxt(dirr[i] +"/tdata.dat")
    time = tdat[:,0]
    u2   = tdat[:,4]
    D = 1.0
    Pe = 1/D
    exp_u2[i] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/(tau)
    gamma   = np.sqrt(1j*omega/Sc)
    gamma_c = np.conj(gamma)
    a       = np.sqrt(1j*omega/D)
    a_c     = np.conj(a)
    rho = a 
    rho_c = a_c
    num_D_para[i] = integrate.trapz(np.loadtxt(dirr[i] +"/tdata.dat")[-timesteps:, 8], np.loadtxt(dirr[i]+"/tdata.dat")[-timesteps:, 0])/tau
    D_ana[i] = 1 + np.real(Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*omega*omega*omega*(Sc*Sc-1))*(1/(gamma*np.tanh(gamma)) + 1j/(gamma*np.tanh(gamma_c)) - 1/(rho*np.tanh(rho)) - 1j/(rho*np.tanh(rho_c))))
    D2 = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
    D_ana[i] = D2
    plt.plot(np.loadtxt(dirr[i] +"/tdata.dat")[-timesteps:, 0], np.loadtxt(dirr[i] +"/tdata.dat")[-timesteps:, 8]-D2)


    amp = ( (1/(rho*tanh(rho)) - 1/(sinh(rho)**2) + 1/(gamma*tanh(gamma)) - 1/(sinh(gamma)**2) - 4*(rho/tanh(rho) - gamma/tanh(gamma))/(rho*rho-gamma*gamma))/2)
    amp *= Pe*Pe*F0*F0*tanh(gamma)*tanh(gamma)/(4*gamma*gamma*(rho*rho-gamma*gamma)**2)
    plt.plot(np.loadtxt(dirr[i] +"/tdata.dat")[-timesteps:, 0], np.ones(timesteps)*np.sqrt(4*np.real(amp)**2 + 4*np.imag(amp)**2))
    #print(D_ana[i], num_D_para[i], D2)
plt.show()
"""
plt.figure(4)
plt.plot(visc, num_D_para, "ro", label="numerisk brenner")
plt.plot(visc, D_ana, "ko", label="analytisk")
plt.plot(visc, num_D_para, "r-")
plt.plot(visc, D_ana, "k-")
#plt.show()


#simulation parameters
periods      = 1000
datafiles    = periods*25
t = np.linspace(0, tau*periods, datafiles)
half_way     = int(2*datafiles/4)
var    = np.zeros((len(visc), datafiles))
"""
"""
for l in range(len(visc)):
    var = np.load(dirr[l] + "/pos/var.npy")
    mean = np.load(dirr[l]+ "/pos/mean.npy")
    Dm = 1#np.sqrt(exp_u2[l])/Pe
    plt.figure(2)
    plt.plot(t/tau, var/Dm, label="visc="+str(visc[l]))

    plt.figure(3)
    plt.plot(t[1:]/tau, var[1:]/(Dm*t[1:]), label="visc="+str(visc[l]))
    plt.plot(t[half_way:]/tau, var[half_way:]/(Dm*t[half_way:]))

    RW_D_para[0,l] = np.mean(var[half_way:]/(2*Dm*t[half_way:]))
    RW_D_para[1,l] = np.std( var[half_way:]/(2*Dm*t[half_way:]))

"""
"""
data1 = np.load(dirr[0]+"/var_over_t.npy")#/visc_1.5_var_over_t.npy")
data2 = np.load(dirr[1]+"/var_over_t.npy")#/visc_3.0_var_over_t.npy")
data3 = np.load(dirr[2]+"/var_over_t.npy")#/visc_5.0_var_over_t.npy")
cutoff = int(0.8*len(data1))

plt.figure(5)
plt.plot(data1, label="nu=1.5")
plt.plot(data2, label="nu=3.0")
plt.plot(data3, label="nu=5.0")
plt.legend(loc="best")
plt.xlabel("data points (time)")
plt.ylabel("Effective D")
RW_D_para[0, 0] = (np.mean(data1[cutoff:]))
RW_D_para[0, 1] = (np.mean(data2[cutoff:]))
RW_D_para[0, 2] = (np.mean(data3[cutoff:]))

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

#simulation parameters
periods      = 1000
datafiles    = periods*25
t = np.linspace(0, tau*periods, datafiles)
half_way     = int(2*datafiles/4)
var    = np.zeros((len(visc), datafiles))


for l in range(len(visc)):
    var = np.load(dirr[l] + "/pos/var.npy")
    mean = np.load(dirr[l]+ "/pos/mean.npy")
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
"""