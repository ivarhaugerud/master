import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata
import matplotlib.ticker as tick 
from scipy import integrate
import seaborn as sns
import matplotlib
import os
import scipy.integrate as sci


plt.style.use("bmh")
sns.color_palette("hls", 1)

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

dirr = "results_oscwavychannel/"
base = "data_test/vary_geometry/"

#simulation paramters
#dt = 0.02
#tau = 5.0
#timesteps = int((tau/0.02))

epsilon = np.arange(0.1, 0.601, 0.10)
Lx = np.array([12.0])#, 9.0, 6.0, 3.0])
kappa = 2*np.pi/Lx 
kappas = kappa
tau = 3.0
omega = 2*np.pi/tau 

U = np.zeros((len(kappa), len(epsilon)))
difference = np.zeros(np.shape(U))
T = 750 

plt.figure(1)
kappa_cont = np.linspace(min(kappa), max(kappa), int(1e4))

for i in range(len(kappa)):
	for j in range(len(epsilon)):
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau3.0_eps"+str(epsilon[j])[:3]+"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt0.004/tdata.dat")
		U[i, j] = sci.trapz(  data[:, 4][-T:],  data[:, 0][-T:] )/tau
		difference[i, j] = abs(U[i, j] - sci.trapz(  np.trim_zeros(data[:, 4])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau)/U[i,j]

		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 4])[-T:])
		plt.xlabel(r" Time [periods]", fontsize=8)
		plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.show()

F0 = 12/1.2
omega = 2*np.pi/tau
nu = 1.2
D = 1#0.3
Sc = nu
gamma = np.sqrt(1j*omega/Sc)

#epsilon = np.logspace(-2, np.log10(0.6), 8)
#epsilon = np.arange(0, 0.6, 0.05)#[0.01, 0.01794823, 0.03221389, 0.05781823, 0.10377349, 0.185]
#print(epsilon)
kappa = kappa[0]
Lx = 2*np.pi/kappa


###
#ANALYTIC RESULTS
###

pi = np.pi 
Nt = 120
T = np.linspace(0, 2*np.pi/omega, Nt)
eta = np.linspace(0, 2*np.pi/kappa, 100)
xi = np.linspace(-1, 1, 350)

u_x = np.zeros((len(T), len(xi), len(eta), 3))
u_y = np.zeros((len(T), len(xi), len(eta), 3))

kappa_prime = np.sqrt(1j*omega/Sc + kappa*kappa)
P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
P_1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))

from cmath import *
Ax = sqrt(gamma**2 + 4*kappa**2)*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(4*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
Ay = kappa*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(2*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
psi_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

#zeroth order
xi = np.linspace(-1, 1, int(1e5))
u_x0 = (F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma))

#first order
u_x1 = ((  (P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_prime*xi)/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi)/np.cosh(kappa_prime) - xi*np.sinh(gamma*xi)/np.sinh(gamma))  ))
u_y1 = (kappa*P1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_prime*xi)/np.sinh(kappa_prime) - np.sinh(kappa*xi)/np.sinh(kappa))

#second order
u_x2 = P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma)


U_ana = np.zeros(len(epsilon))
for e in range(len(epsilon)):
	eps = epsilon[e]
	kin = u_x0*np.conjugate(u_x0) + eps*eps*(np.conjugate(u_x0)*u_x1 + 2*np.conjugate(u_x0)*u_x2 + np.conjugate(u_y1)*u_y1/2)
	U_ana[e] = sci.trapz( kin, xi)/2

plt.plot(epsilon, U_ana, "-", label="Analytic")
for i in range(len(kappas)):
	plt.plot(epsilon, U[i,:], "ko")
#plt.plot(epsilon[0], exp_u2[0], "ko", label="Numerical")
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=14)
plt.ylabel(r"Kinetic energy of fluid $\langle u^2 \rangle/2$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.show()

for i in range(len(kappas)):
	plt.plot(epsilon, abs((U_ana-U[i,:])), "ko")
plt.plot(epsilon[:12], epsilon[:12]**4)
plt.yscale("log")
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=14)
plt.ylabel(r"Relative difference analytic and numerical kinetic energy $\langle u^2 \rangle/2$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.show()
