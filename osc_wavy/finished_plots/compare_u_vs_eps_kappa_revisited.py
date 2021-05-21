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
from cmath import *

root = "../../../master_latex/results/"
plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	


Lx = np.array([17.95, 31.41, 41.89])
Lx2 = np.array([17.95, 41.88])
kappas  = 2*np.pi/Lx
kappas2 = 2*np.pi/Lx2
epsilon = np.array([0.0, 0.0178, 0.0316, 0.056, 0.1, 0.177, 0.245, 0.31])
base = "../data_test/benchmark_3/"

U           = np.zeros((len(epsilon), len(kappas)))
U2          = np.zeros(np.shape(U))

T = int(6/(0.003))
tau = 6
omega = 2*np.pi/tau

for i in range(len(epsilon)):
	for j in range(len(kappas)):
		res = int(170)
		data = np.loadtxt(base+"Lx"+str(Lx[j])+"_tau6.0_eps"+str(epsilon[i])+"_nu5.0_D25.0_fzero0.0_fone250.0_res"+str(res)+"_dt0.003/tdata.dat")
		U[i, j] = integrate.trapz(  data[-T:, 4],  data[-T:, 0] )/tau

#for i in range(len(epsilon)):#
#	for j in range(len(kappa2)):
#		data = np.loadtxt(base+"Lx"+str(Lx2[j])+"_tau6.0_eps"+str(epsilon[i])+"_nu5.0_D4.0_fzero0.0_fone250.0_res"+str(res)+"_dt0.003/tdata.dat")
#		D2[i, j] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau


tau = 6.0
nu  = 5
F0  = 250/nu 
omega = 2*np.pi/tau
gamma = np.sqrt(1j*omega/nu)
pi = np.pi 

###
#ANALYTIC RESULTS
###

xi = np.linspace(-1, 1, int(1e5))
u  = np.zeros((len(epsilon), len(kappas)))

for j in range(len(kappas)):
	kappa = kappas[j]
	kappa_p = np.sqrt(1j*omega/nu + kappa*kappa)
	kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
	"""
	P1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_prime*np.tanh(kappa)/(kappa*np.tanh(kappa_prime)))
	kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
	kappa_pp = np.sqrt(4*kappa*kappa+gamma*gamma)
	P_1 = (gamma*F0*np.tanh(gamma)/(kappa*np.cosh(kappa)))/(1-kappa_p*np.tanh(kappa)/(kappa*np.tanh(kappa_p)))

	
	Ax = sqrt(gamma**2 + 4*kappa**2)*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(4*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
	Ay = kappa*(F0*gamma**2*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*gamma**2*sinh(kappa)*sinh(2*kappa)*tanh(kappa_p) + 2*P_1*kappa**2*cosh(kappa)*cosh(2*kappa)*tanh(kappa_p) - 2*P_1*kappa*kappa_p*sinh(kappa)*cosh(2*kappa))/(2*gamma**2*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
	psi_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

	ux0 = 0.5*F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)
	ux1 = 0.5*((P1*kappa*np.cosh(kappa)/(gamma*gamma))*(np.cosh(kappa*xi)/np.cosh(kappa) - np.cosh(kappa_prime*xi)/np.cosh(kappa_prime))   + (F0*np.tanh(gamma)/(gamma))*(np.cosh(kappa_prime*xi)/np.cosh(kappa_prime) - xi*np.sinh(gamma*xi)/np.sinh(gamma)))
	uy1 = 0.5*(kappa*P1*np.sinh(kappa)/(gamma*gamma)*( np.sinh(kappa_prime*xi)/np.sinh(kappa_prime) - np.sinh(kappa*xi)/np.sinh(kappa)))
	uy2 = 0.5*((-P_1*kappa**2*xi*np.cosh(kappa*xi)/(2*gamma**2) + P_1*kappa*kappa_p*xi*np.sinh(kappa)*np.cosh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.sinh(2*kappa*xi)/gamma**2) + Ay*np.sinh(kappa_pp*xi))
	ux2 = 0.5*(F0*xi**2*np.cosh(gamma*xi)/(4*np.cosh(gamma)) - P_1*kappa**2*xi*np.sinh(kappa*xi)/(2*gamma**2) + P_1*kappa_p*kappa_p*np.sinh(kappa)*xi*np.sinh(kappa_p*xi)/(2*gamma**2*np.sinh(kappa_p)) - 2*kappa*psi_2*np.cosh(2*kappa*xi)/gamma**2 + Ax*np.cosh(kappa_pp*xi))
	ux2_no_eta = 0.5*(P_1*kappa*kappa*np.sinh(kappa)*(xi*np.sinh(kappa*xi)/np.sinh(kappa) - np.cosh(gamma*xi)/np.cosh(gamma))/(2*gamma*gamma) + F0*np.cosh(gamma*xi)*(1-xi*xi)/(4*np.cosh(gamma)) + P_1*kappa_p*kappa_p*np.sinh(kappa)*(np.cosh(gamma*xi)/np.cosh(gamma)-xi*np.sinh(kappa_p*xi)/np.sinh(kappa_p))/(2*gamma*gamma))
	"""
	from numpy import *
	ux0 = 0.5*F0*(1-np.cosh(gamma*xi)/np.cosh(gamma))/(gamma*gamma)
	P_1 = (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p)))
	P_2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))

	uy1 = 0.5*((kappa*P_1*sinh(kappa)/(gamma*gamma))*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa)))/2
	ux1 = 0.5*((P_1*kappa*cosh(kappa)/(gamma*gamma))*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + (F0*tanh(gamma)/gamma)*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))/2
	uy1 = 0.5*((kappa*P_1*sinh(kappa)/(gamma*gamma))*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa)))/2
	ux2_no_eta = 0.5*(P_1*sinh(kappa)*(kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa) - kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p) + gamma*gamma*cosh(gamma*xi)/cosh(gamma))/(2*gamma*gamma) + cosh(gamma*xi)*F0*(1-xi*xi)/(4*cosh(gamma)))
	ux2 = 0.5* (P_1*sinh(kappa)*(kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p) -kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa) - gamma*gamma*cosh(kappa_pp*xi)/cosh(kappa_pp))/(2*gamma*gamma) + F0*(xi*xi*cosh(gamma*xi)/cosh(gamma) - cosh(kappa_pp*xi)/cosh(kappa_pp))/4 - P_2*cosh(2*kappa)*(cosh(2*kappa*xi)/cosh(2*kappa) - cosh(kappa_pp*xi)/cosh(kappa_pp)) )
	uy2 = 0.5* ((P_1*kappa*sinh(kappa)/(2*gamma*gamma*tanh(kappa_p)))*(kappa_p*xi*cosh(kappa_p*xi)/cosh(kappa_p) - kappa_p*sinh(kappa_pp*xi)/sinh(kappa_pp) -kappa*(tanh(kappa_p)/tanh(kappa))*(xi*cosh(kappa*xi)/cosh(kappa) - sinh(kappa_pp*xi)/sinh(kappa_pp))) - P_2*sinh(2*kappa)*(sinh(2*kappa*xi)/sinh(2*kappa) - sinh(kappa_pp*xi)/sinh(kappa_pp)) )

	u[:, j] =  2*epsilon*epsilon*0.5*integrate.trapz( 0.5*2*ux0*np.conjugate(ux1) + 0.5*ux1*np.conjugate(ux1) + 2*ux2_no_eta*np.conjugate(ux0) + 0.5*uy1*np.conjugate(uy1) , xi)#  + integrate.trapz(ux0*np.conjugate(ux0), xi)

for i in range(len(kappas)):
	plt.plot(epsilon, u[:, i], color="C"+str(i))
	plt.plot(epsilon, U[:,i]-U[0,0], "o", markersize=3, color="C"+str(i))

plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Kinetic energy of fluid $\langle u^2 \rangle^{(2)}$", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#plt.axis([-0.05, 0.52, -0.1, 0.45])
plt.show()

for i in range(len(kappas)):
	plt.plot(epsilon, (abs(u[:,i]-(U[:,i]-U[0,0])))/U[:,i], "o", markersize=3, label=r"$\kappa=$"+str(kappas[i])[:5], color=sns.color_palette()[i])

#epsilon = np.linspace(0.06, max(epsilon), int(1e3))
plt.plot(epsilon, epsilon**4/(1-epsilon), "k", label=r"$\epsilon^4$")
plt.yscale("log")
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=8)
plt.ylabel(r"Difference kinetic energy $\langle u^2 \rangle^{(2)}$", fontsize=8)
plt.legend(loc="best", fontsize=8, ncol=2)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xscale("log")
#filename = root+"figures/comparison_numeric_analytic_difference4.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
