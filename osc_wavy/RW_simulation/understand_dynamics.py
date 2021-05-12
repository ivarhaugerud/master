import numpy as np
import os
import sys
import matplotlib.pyplot as plt 
import scipy.integrate as sci

def velocity(xi, eta, t):
    xi = xi/(1+epsilon*np.sin(kappa*eta))
    uy = (-P_1*kappa**2*xi*np.cosh(kappa*xi)/(2*gamma**2) + P_1*kappa*kappa_p*xi*np.sinh(kappa)*np.cosh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/gamma**2) + Ay*np.sinh(kappa_pp*xi)
    ux = F0*xi**2*np.cosh(gamma*xi)/(4*np.cosh(gamma)) - P_1*kappa**2*xi*np.sinh(kappa*xi)/(2*gamma**2) + P_1*kappa_p*kappa_p*np.sinh(kappa)*xi*np.sinh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/gamma**2 + Ax*np.cosh(kappa_pp*xi)
    ux_no_eta = P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma)

    ux_full = np.real(np.exp(1j*omega*t)*F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma))+ epsilon*np.real(np.exp(1j*omega*t)*np.sin(kappa*eta)*(  (P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_prime*xi)/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi)/np.cosh(kappa_prime) - xi*np.sinh(gamma*xi)/np.sinh(gamma))  ))+ epsilon*epsilon*np.real(np.exp(1j*omega*t)*ux*np.cos(2*kappa*eta)) + np.real(np.exp(1j*omega*t)*ux_no_eta)
    uy_full = epsilon*np.real((np.exp(1j*omega*t)*np.cos(kappa*eta)*kappa*P1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_prime*xi)/np.sinh(kappa_prime) - np.sinh(kappa*xi)/np.sinh(kappa))) + epsilon*epsilon*np.real(np.exp(1j*omega*t)*uy*np.sin(2*kappa*eta))
    return np.array([ux_full, uy_full])

#simulation parameters
dt           = 0.5
tau          = 40.0
omega        = 2*np.pi/tau
kappa        = 2*np.pi/12.5
periods      = 5
timesteps    = int(tau/dt)
datafiles    = 5
N            = int(2*1e3)
RW_timesteps = 10
D            = 0.1
alpha        = np.sqrt(2*D*dt/RW_timesteps)
var          = np.zeros(int(periods*timesteps))
t            = np.linspace(0, tau-dt, timesteps)+tau/2
pos_saves    = np.linspace(0, timesteps*periods, datafiles, dtype="int")
saves_counter= 0
print("Diffusive transport:", alpha)
print("Advective transport:", dt)

epsilon = 0.3
l = 2*np.pi/kappa
nu = 1000
F0 = 2.5*1000/nu 
from cmath import *
gamma = np.sqrt(1j*omega/nu)
kappa_prime = np.sqrt(gamma*gamma + kappa*kappa)
P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
P_1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))
Ax = sqrt(gamma**2 + 4*kappa**2)*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(4*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
Ay = kappa*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(2*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
psi_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

test_y = np.linspace(-7*np.pi/kappa, 7*np.pi/kappa, 1000)
test_x = np.linspace(-1-epsilon, 1+epsilon, 500)
u = np.zeros((len(test_x), len(test_y)))

X, Y = np.meshgrid(test_y,test_x)
xi = np.linspace(-1, 1, int(1e4))

dirr = "data/understand_dynamics_analytic_again_eps"+str(epsilon)+"_kappa"+str(kappa)[:5]
os.system("mkdir " + dirr)

Nt = int(periods*timesteps)
prev_pos = np.zeros((Nt, 2, N))
pos      = np.zeros((Nt, 2, N))
U        = np.zeros(Nt)
U[0] = sci.trapz((velocity(xi, 0, t[0])[0]),xi)/2

pos[0, 1, :]     = np.linspace(-1-epsilon+0.01, 1+epsilon-0.01, N)
pos[0, 0, :]     = l/4

for k in range(1, int(periods*timesteps)):
    pos[k, :, :] +=  pos[k-1, :, :] + dt * velocity(pos[k-1, 1, :], pos[k-1, 0, :], t[(k+timesteps)%timesteps])
    U[k] = sci.trapz((velocity(xi, 0, t[(k+timesteps)%timesteps])[0]),xi)/2
    pos[k, :, np.where( pos[k, 1, :] >   1+epsilon*np.sin(kappa*pos[k, 0,:]))] = pos[k-1, :, np.where( pos[k, 1, :] >  1+epsilon*np.sin(kappa*pos[k, 0,:]))] #checks if y-coordinate outside
    pos[k, :, np.where( pos[k, 1, :] <  -1-epsilon*np.sin(kappa*pos[k, 0,:]))] = pos[k-1, :, np.where( pos[k, 1, :] < -1-epsilon*np.sin(kappa*pos[k, 0,:]))] #checks if y-coordinate outside

    for l in range(RW_timesteps):
        pos[k, :, :] += alpha*np.random.normal(loc=0, scale=1, size=(2,N))
        pos[k, :, np.where( pos[k, 1, :] >   1+epsilon*np.sin(kappa*pos[k, 0,:]))] = pos[k-1, :, np.where( pos[k, 1, :] >  1+epsilon*np.sin(kappa*pos[k, 0,:]))] #checks if y-coordinate outside
        pos[k, :, np.where( pos[k, 1, :] <  -1-epsilon*np.sin(kappa*pos[k, 0,:]))] = pos[k-1, :, np.where( pos[k, 1, :] < -1-epsilon*np.sin(kappa*pos[k, 0,:]))] #checks if y-coordinate outside
        prev_pos = np.copy(pos)
    plt.clf()
    
    if k % 2 == 0:
        plt.scatter(pos[k, 0, :], pos[k, 1, :], s=0.5)
        plt.plot(test_y, 1+epsilon*np.sin(kappa*test_y), "k")
        plt.plot(test_y,-1-epsilon*np.sin(kappa*test_y), "k")
        #plt.axis("equal")
        plt.pause(0.01)
    var[k] = np.var(pos[0, :])
    
    
    if k % pos_saves[saves_counter] == 0:
        print("timestep:", k, " out of ", int(periods*timesteps), " proportion: ", k/int(periods*timesteps))
        np.save(dirr+"/pos_"+str(pos_saves[saves_counter]), pos)
        saves_counter += 1
    

#plt.plot(U)
#plt.show()
"""
t = np.linspace(0, periods*tau, Nt)
UMI = np.argmax(U[int(periods*timesteps/2):])+int(periods*timesteps/2)
PI  = int(timesteps/4)
plt.plot(t, U)
plt.plot(t[UMI-PI:UMI], U[UMI-PI:UMI])
plt.show()
Dx_avg = 2*sci.trapz(U[UMI-PI:UMI], t[UMI-PI:UMI])
print("max delta X in a period:", Dx_avg)
print("Wavelength:", 2*np.pi/kappa)
print("Diffusion length in that time: ", np.sqrt(D*tau))

"""
plt.plot(var)
plt.show()
n = Nt
FPI = np.argmax(abs(pos[-1, 0, :]))
U -= np.min(U)
U /= np.max(U)
test_y = np.linspace(0, np.max(pos[:, 0, FPI]), 1000)

fig = plt.figure()
ax = fig.add_subplot(111)
s = 10 # Segment length

for i in range(0,n-s,s):
    ax.plot(pos[i:i+s+1, 0, FPI],pos[i:i+s+1, 1, FPI],color=(U[i],0.0,U[i]))

#plt.scatter(pos[:, 0, FPI], pos[:, 1, FPI],  c=U/np.max(U), cmap="hsv", s=1)
plt.plot(test_y, 1+epsilon*np.sin(kappa*test_y), "k")
plt.plot(test_y,-1-epsilon*np.sin(kappa*test_y), "k")
#plt.axis("equal")

plt.savefig("figure.png")
plt.show()
