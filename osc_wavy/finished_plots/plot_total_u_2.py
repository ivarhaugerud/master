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

root = "../../../master_latex/results/"
plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	

dirr = "../results_oscwavychannel/run_12_11/"

#simulation paramters
dt = 0.01
tau = 5.0 
T = int(tau/dt)
epsilon = np.arange(0.0, 0.51, 0.05)
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
periods = 2
#epsilon = np.
for i in range(len(epsilon)):
	eps = epsilon[i]
	for j in range(len(kappas)):
		res = int(100*(1+float(eps)))
		try:
			filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(eps)[:4])+"_nu1.2_D0.3_fzero0.0_fone3.0_res"+str(res)+"_dt0.01/tdata.dat"
			tdat = np.loadtxt(dirr + filename)
		except:
			filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(eps)[:3])+"_nu1.2_D0.3_fzero0.0_fone3.0_res"+str(res)+"_dt0.01/tdata.dat"
			tdat = np.loadtxt(dirr + filename)
		t = tdat[:,0]
		u2 = tdat[:,4]
		exp_u2[i, j] = integrate.trapz(u2[-T:], t[-T:])/(tau)

###
#ANALYTIC RESULTS
###

pi = np.pi 
xi = np.linspace(-1, 1, int(1e5))
new_eps = epsilon
u_squared_ana = np.zeros((len(new_eps),len(kappas)))

for i in range(len(new_eps)):
	epsilon = new_eps[i]
	for j in range(len(kappas)):
		kappa = kappas[j]
		kappa_prime = np.sqrt(1j*omega/Sc + kappa*kappa)
		P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
		kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
		kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
		gamma   = np.sqrt(1j*omega/nu)
		rho     = np.sqrt(1j*omega/D)
		kappa_p = np.sqrt(gamma*gamma + kappa*kappa)
		P_1     = ((F0*gamma*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p))) )

		from cmath import *
		Ax = sqrt(gamma**2 + 4*kappa**2)*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(4*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
		Ay = kappa*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(2*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
		psi_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

		uy2 = 0.5*((-P_1*kappa**2*xi*np.cosh(kappa*xi)/(2*gamma**2) + P_1*kappa*kappa_p*xi*np.sinh(kappa)*np.cosh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/gamma**2) + Ay*np.sinh(kappa_pp*xi) )
		ux2_sin = 0.5*(F0*xi**2*np.cosh(gamma*xi)/(4*np.cosh(gamma)) - P_1*kappa**2*xi*np.sinh(kappa*xi)/(2*gamma**2) + P_1*kappa_p*kappa_p*np.sinh(kappa)*xi*np.sinh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/gamma**2 + Ax*np.cosh(kappa_pp*xi))

		#analytic results
		ux0 = F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma)
		#ux1 = 0.5*((  (P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p*xi)/np.cosh(kappa_p))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma))  ))
		#uy1 = 0.5*(kappa*P_1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))
		#ux2 = 0.5*(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))
		ux2 = 0.5*(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))
		uy1 = 0.5*(kappa*P_1*np.sinh(kappa)/(gamma*gamma))*( np.sinh(kappa_p*xi)/np.sinh(kappa_p) - np.sinh(kappa*xi)/np.sinh(kappa))
		ux1 = 0.5*((P_1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_p*xi)/np.cosh(kappa_p))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_p*xi)/np.cosh(kappa_p) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))
		u = ux0*np.conjugate(ux0) + epsilon*epsilon*ux0*np.conjugate(ux1) + 2*epsilon*epsilon*ux0*np.conjugate(ux2) + epsilon*epsilon*uy1*np.conjugate(uy1) + 0.5*epsilon**4*ux2*np.conjugate(ux2) + epsilon**3*ux1*np.conjugate(ux2)
		u += ux2*ux2*0.5*epsilon**4 + 0.5*epsilon**4*uy2*uy2 + 0.5*epsilon**3*uy1*uy2
		u_squared_ana[i, j] = integrate.trapz(u, xi)

print(np.shape(exp_u2))
print(len(kappas))
print(len(kappas))
epsilon = new_eps

for j in range(len(kappas)):
	plt.plot(new_eps, u_squared_ana[:,j], "o", markersize=3, label=r"$\kappa=$"+str(kappas[j])[:5], color=sns.color_palette()[j])
	plt.plot(epsilon, exp_u2[:, j], "-", color=sns.color_palette()[j])

plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Kinetic energy of fluid $\langle u^2 \rangle$", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.axis([-0.02, 0.52, 0.02, 0.37])
#filename = root + "figures/comparison_numeric_analytic.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

for i in range(len(kappas)):
	plt.plot(epsilon, abs((u_squared_ana[:,i]-exp_u2[:,i])/exp_u2[:,i]), "o", markersize=3, label=r"$\kappa=$"+str(kappas[i])[:5], color=sns.color_palette()[i])
	plt.plot(epsilon, abs((u_squared_ana[:,i]-exp_u2[:,i])/exp_u2[:,i]), "-", linewidth=1, color=sns.color_palette()[i])

epsilon = np.linspace(0.03, max(epsilon), int(1e3))
plt.plot(epsilon, epsilon**4/(1-epsilon), "k")
plt.yscale("log")
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Relative difference kinetic energy $\langle u^2 \rangle$", fontsize=8)
plt.legend(loc="best", fontsize=8, ncol=2)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)

#filename = root+"figures/comparison_numeric_analytic_difference2.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
