import numpy as np
import os
import sys
import matplotlib.pyplot as plt 

def velocity(xi, eta, t):
    xi = xi/(1+epsilon*np.sin(kappa*eta))
    uy = (-P_1*kappa**2*xi*np.cosh(kappa*xi)/(2*gamma**2) + P_1*kappa*kappa_p*xi*np.sinh(kappa)*np.cosh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/gamma**2) + Ay*np.sinh(kappa_pp*xi)
    ux = F0*xi**2*np.cosh(gamma*xi)/(4*np.cosh(gamma)) - P_1*kappa**2*xi*np.sinh(kappa*xi)/(2*gamma**2) + P_1*kappa_p*kappa_p*np.sinh(kappa)*xi*np.sinh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/gamma**2 + Ax*np.cosh(kappa_pp*xi)
    ux_no_eta = P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma)

    ux_full = np.real(np.exp(1j*omega*t)*F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma))+ epsilon*np.real(np.exp(1j*omega*t)*np.sin(kappa*eta)*(  (P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_prime*xi)/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi)/np.cosh(kappa_prime) - xi*np.sinh(gamma*xi)/np.sinh(gamma))  ))+ epsilon*epsilon*np.real(np.exp(1j*omega*t)*ux*np.cos(2*kappa*eta)) + np.real(np.exp(1j*omega*t)*ux_no_eta)
    uy_full = epsilon*np.real((np.exp(1j*omega*t)*np.cos(kappa*eta)*kappa*P1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_prime*xi)/np.sinh(kappa_prime) - np.sinh(kappa*xi)/np.sinh(kappa))) + epsilon*epsilon*np.real(np.exp(1j*omega*t)*uy*np.sin(2*kappa*eta))
    return np.array([ux_full, uy_full])

#simulation parameters
dt           = 0.01
tau          = 3.0 
periods      = 50
timesteps    = int(tau/dt)
datafiles    = 1000
N            = int(15*1e3)
RW_timesteps = 20
D            = 0#1e-6
alpha        = np.sqrt(2*D*dt/RW_timesteps)
var          = np.zeros(int(periods*timesteps))
t            = np.linspace(0, tau-dt, timesteps)+tau/2
pos_saves    = np.linspace(0, timesteps*periods, datafiles, dtype="int")
saves_counter= 0
print("Diffusive transport:", alpha)
print("Advective transport:", dt)

epsilons = np.arange(0.1, 0.31, 0.10)
omega = 2*np.pi/tau
kappa = float(sys.argv[1])
l = 2*np.pi/kappa
nu = 1.2
F0 = 12/nu 
Sc = nu 
from cmath import *
gamma = np.sqrt(1j*omega/Sc)
kappa_prime = np.sqrt(gamma*gamma + kappa*kappa)
P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
P_1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))
Ax = sqrt(gamma**2 + 4*kappa**2)*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(4*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
Ay = kappa*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(2*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
psi_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

epsilon = 0.25
test_y = np.linspace(-4*np.pi/kappa, 5*np.pi/kappa, 1000)
test_x = np.linspace(-1-epsilon, 1+epsilon, 500)
u = np.zeros((len(test_x), len(test_y)))

X, Y = np.meshgrid(test_y,test_x)

for i in range(len(epsilons)):
    epsilon = epsilons[i]
    dirr = "data/analytic_D0_eps"+str(epsilon)+"_kappa"+str(kappa)
    os.system("mkdir " + dirr)

    prev_pos = np.zeros((2, N))
    pos      = np.zeros((2, N))

    pos[1,:]      = np.linspace(-1-epsilon+0.01, 1+epsilon-0.01, N)#np.random.uniform(-1-epsilon+ 0.01, 1+epsilon-0.01, N)
    pos[0, :]     = l/4
    prev_pos[1,:] = np.copy(pos[1,:])

    for k in range(int(periods*timesteps)):
        pos[:, :] +=  dt * velocity(pos[1, :], pos[0, :], t[(k+timesteps)%timesteps])

        pos[:, np.where( pos[1, :] >   1+epsilon*np.sin(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.sin(kappa*pos[0,:]))] #checks if y-coordinate outside
        pos[:, np.where( pos[1, :] <  -1-epsilon*np.sin(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.sin(kappa*pos[0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)

        for l in range(RW_timesteps):
            pos[:, :] += alpha*np.random.normal(loc=0, scale=1, size=(2,N))
            pos[:, np.where( pos[1, :] >   1+epsilon*np.sin(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] >  1+epsilon*np.sin(kappa*pos[0,:]))] #checks if y-coordinate outside
            pos[:, np.where( pos[1, :] <  -1-epsilon*np.sin(kappa*pos[0,:]))] = prev_pos[:, np.where( pos[1, :] < -1-epsilon*np.sin(kappa*pos[0,:]))] #checks if y-coordinate outside
            prev_pos = np.copy(pos)
        #plt.clf()
        #if k % 5 == 0:
        #    plt.scatter(pos[0, :], pos[1, :], s=0.5)
        #    plt.plot(test_y, 1+epsilon*np.sin(kappa*test_y), "k")
        #    plt.plot(test_y,-1-epsilon*np.sin(kappa*test_y), "k")
        #    plt.pause(0.01)
        var[k] = np.var(pos[0, :])
        if k % pos_saves[saves_counter] == 0:
            print("timestep:", k, " out of ", int(periods*timesteps), " proportion: ", k/int(periods*timesteps))
            np.save(dirr+"/pos_"+str(pos_saves[saves_counter]), pos)
            saves_counter += 1

    t = np.linspace(0, periods*tau, periods*timesteps)
    np.save(dirr + "/var_over_2Dt", var[1:]/(2*D*t[1:]))
    print("DONE WITH: ", dirr, np.mean(var[-int(periods*timesteps/2):]/(2*D*t[-int(periods*timesteps/2):])))
    plt.scatter(pos[0, :], pos[1, :])
    plt.show()