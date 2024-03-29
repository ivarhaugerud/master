import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci
import os
import seaborn as sns
from scipy import integrate
import h5py
from scipy.interpolate import griddata

def velocity(xi, t, omega, gamma):
    #print(gamma, omega)
    return np.real(F0/(gamma*gamma)*(1-np.cosh(gamma*xi)/np.cosh(gamma))*np.exp(1j*omega*t))

#simulation parameters
dt           = 0.006
tau          = 3.0 
periods      = 5000
timesteps    = int(tau/dt)
datafiles    = periods*25
N            = int(4*1e3)
RW_timesteps = 30
Pe           = 1
#geometry parameters
epsilon = 0.0
omega = 2*np.pi/tau 
visc = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
exp_u2 = np.zeros(len(visc))
l = 6.28
"""
dirr = []
for i in range(len(visc)):
    dirr.append("flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

    tdat = np.loadtxt(dirr[i] +"tdata.dat")
    exp_u2[i] = integrate.trapz(tdat[-timesteps:, 4], tdat[-timesteps:, 0])/(tau)
"""

for i in range(len(visc)):
    nu = visc[i]
    F0 = 12/nu 
    Sc = nu 
    gamma = np.sqrt(1j*omega/Sc)

    prev_pos = np.zeros((2, N))
    pos      = np.zeros((2, N))
    pos[1,:]      = np.random.uniform(-0.99, 0.99, N)
    prev_pos[1,:] = np.copy(pos[1,:])

    U = np.sqrt(exp_u2[i])
    D  = 0.1 #U/Pe
    alpha = np.sqrt(2*D*dt/RW_timesteps)

    var = np.zeros(int(periods*timesteps))
    t = np.linspace(0, tau-dt, timesteps)

    for k in range(int(periods*timesteps)):
        pos[:, :] +=  dt * velocity(pos[1, :], t[(k+timesteps)%timesteps], omega, gamma)

        for l in range(RW_timesteps):
            pos[:, :] += alpha*np.random.normal(loc=0, scale=1, size=(2,N))
            pos[:, np.where( abs(pos[1, :]) > 1)] = prev_pos[:, np.where( abs(pos[1, :]) >  1)]
            prev_pos = np.copy(pos)

        var[k] = np.var(pos[0, :])

    t = np.linspace(0, periods*tau, periods*timesteps)
    np.save("data/var_over_2Dt_nu"+str(nu)+"_D0.1", var[1:]/(2*D*t[1:]))
    #plt.plot(t[1:], var[1:]/(2*D*t[1:]))
    #plt.show()
    print("DONE WITH RUN FOR NU: ", visc[i], np.mean(var[-int(periods*timesteps/2):]/(2*D*t[-int(periods*timesteps/2):])))