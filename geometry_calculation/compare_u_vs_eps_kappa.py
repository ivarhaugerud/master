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

plt.style.use("bmh")
sns.color_palette("hls", 1)

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

dirr = "results_oscwavychannel/run_12_11/"

#simulation paramters
dt = 0.01
tau = 5.0 

epsilon = np.arange(0.0, 0.61, 0.1)
kappas  = np.array([0.1, 1.0, 1.5, 3, 5])
kappas = np.array([0.1,  0.7, 1.5, 3, 5])
Lx = 2*np.pi/kappas

omega = 2*np.pi/tau
nu = 1.2
D = 1
f1 = 3
F0 = f1/nu
Sc = nu
gamma = np.sqrt(1j*omega/Sc)

exp_u2 = np.zeros((len(epsilon), len(kappas)))
periods = 3

for i in range(len(epsilon)):
	eps = epsilon[i]
	for j in range(len(kappas)):
		res = int(100*(1+float(eps)))
		filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(eps)[:3])+"_nu1.2_D0.3_fzero0.0_fone3.0_res"+str(res)+"_dt0.01/tdata.dat"
		#try:
		tdat = np.loadtxt(dirr + filename)
		#except:
		#	print("no file for kappa="+str(kappas[j]) + " epsilon =" + str(eps), filename)
		t = tdat[:,0]
		u2 = tdat[:,4]
		start_index = np.argmin(abs(t -(periods-2.25)*tau))
		end_index   = np.argmin(abs(t -(periods-0.25)*tau))
		#plt.plot(t, u2)
		#plt.plot(t[start_index:end_index], u2[start_index:end_index], "--")
		exp_u2[i, j] = integrate.trapz(u2[start_index:end_index], t[start_index:end_index])/(2*tau)
	#plt.show()

###
#ANALYTIC RESULTS
###

pi = np.pi 
Nt = 200
T = np.linspace(0, 2*np.pi/omega, Nt)
N_eta = 220
xi = np.linspace(-1, 1, 240)

u_x = np.zeros((len(T), len(xi), N_eta, 3))
u_y = np.zeros((len(T), len(xi), N_eta, 3))
new_eps = epsilon#	np.linspace(0, max(epsilon), 20)
u_squared_ana = np.zeros((len(new_eps),len(kappas)))

for j in range(len(kappas)):
	kappa = kappas[j]
	eta = np.linspace(0, 2*np.pi/kappa, N_eta)
	kappa_prime = np.sqrt(1j*omega/Sc + kappa*kappa)
	P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
	kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
	kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
	P_1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))

	from cmath import *
	Ax = sqrt(gamma**2 + 4*kappa**2)*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(4*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
	Ay = kappa*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(2*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
	psi_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

	uy = (-P_1*kappa**2*xi*np.cosh(kappa*xi)/(2*gamma**2) + P_1*kappa*kappa_p*xi*np.sinh(kappa)*np.cosh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/gamma**2) + Ay*np.sinh(kappa_pp*xi)
	ux = F0*xi**2*np.cosh(gamma*xi)/(4*np.cosh(gamma)) - P_1*kappa**2*xi*np.sinh(kappa*xi)/(2*gamma**2) + P_1*kappa_p*kappa_p*np.sinh(kappa)*xi*np.sinh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/gamma**2 + Ax*np.cosh(kappa_pp*xi)
	ux_no_eta = P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma)


	for t in range(len(T)):
		for x in range(len(xi)):
			#zeroth order
			u_x[t,x,:, 0] = np.real(np.exp(1j*omega*T[t])*F0*(1-np.cosh(gamma*xi[x])/np.cosh(gamma))/(gamma*gamma))

			#first order
			u_x[t,x,:, 1] = np.real(np.exp(1j*omega*T[t])*np.sin(kappa*eta)*(  (P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi[x])/np.cosh(kappa) - np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi[x])/np.cosh(kappa_prime) - xi[x]*np.sinh(gamma*xi[x])/np.sinh(gamma))  ))
			u_y[t,x,:, 1] = np.real((np.exp(1j*omega*T[t])*np.cos(kappa*eta)*kappa*P1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_prime*xi[x])/np.sinh(kappa_prime) - np.sinh(kappa*xi[x])/np.sinh(kappa)))

			#second order
			u_x[t,x,:, 2]  = np.real(np.exp(1j*omega*T[t])*ux[x]*np.cos(2*kappa*eta)) + np.real(np.exp(1j*omega*T[t])*ux_no_eta[x])
			u_y[t,x,:, 2]  = np.real(np.exp(1j*omega*T[t])*uy[x]*np.sin(2*kappa*eta))


	u = np.zeros((len(T), len(xi), len(eta), 3))
	first_int  = np.zeros((len(T), len(eta)))
	second_int = np.zeros(len(T))

	after_xi_integral = np.zeros((len(T),len(eta)))
	after_xi_eta_integral = np.zeros(len(T))

	for e in range(len(new_eps)):
		eps = new_eps[e]
		for i in range(len(T)):
			#plt.clf()
			u[i, :, :, 0] = u_x[i, :, :, 0] + eps*u_x[i, :, :, 1] + eps*eps*u_x[i, :, :, 2]
			u[i, :, :, 1] =                   eps*u_y[i, :, :, 1] + eps*eps*u_y[i, :, :, 2]

			u[i, :, :, 2] = u[i,:,:,0]*u[i,:,:,0] + u[i,:,:,1]*u[i,:,:,1]

			for k in range(len(eta)):
				after_xi_integral[i,k] = integrate.trapz(u[i,:,k,2]*(1+eps*np.sin(kappa*eta[k])), xi)/2

			after_xi_eta_integral[i] = integrate.trapz(after_xi_integral[i,:], eta)/(2*np.pi/kappa)

		u_squared_ana[e, j] = integrate.trapz(after_xi_eta_integral, T)/(2*pi/omega)

for j in range(len(kappas)):
	plt.plot(new_eps, u_squared_ana[:,j], "-", label="$\kappa=$"+str(kappas[j])[:5], color=sns.color_palette()[j])
	plt.plot(epsilon, exp_u2[:, j], "o", color=sns.color_palette()[j])

plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=14)
plt.ylabel(r"Kinetic energy of fluid $\langle u^2 \rangle/2$", fontsize=14)
plt.legend(loc="best", fontsize=12)

filename = "figures/comparison_numeric_analytic.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

for i in range(len(kappas)):
	plt.plot(epsilon, abs((u_squared_ana[:,i]-exp_u2[:,i])/exp_u2[:,i]), "o", label="$\kappa=$"+str(kappas[i])[:5], color=sns.color_palette()[i])
	plt.plot(epsilon, abs((u_squared_ana[:,i]-exp_u2[:,i])/exp_u2[:,i]), "--", linewidth=1, color=sns.color_palette()[i])

epsilon = np.linspace(0.01, max(epsilon), int(1e3))
plt.plot(epsilon, epsilon**4/(1-epsilon))
plt.yscale("log")
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=14)
plt.ylabel(r"Relative difference analytic and numerical kinetic energy $\langle u^2 \rangle/2$", fontsize=14)
plt.legend(loc="best", fontsize=12)

filename = "figures/comparison_numeric_analytic_difference.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
plt.show()
