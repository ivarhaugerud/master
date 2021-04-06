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

dirr = "results_oscwavychannel/run_06_12/"

#simulation paramters
dt = 0.01
tau = 5.0 
epsilon = np.arange(0.0, 0.31, 0.05)
kappas  = np.array([1.2,  1.4, 1.6, 1.8, 2.0])
Lx = 2*np.pi/kappas

omega = 2*np.pi/tau
nu = 1
D =  1
f1 = 10
F0 = f1/nu
Sc = 1/nu
Pe = F0/(nu*D) 
gamma = np.sqrt(1j*omega/Sc)
gamma_c = np.conj(gamma)
a = np.sqrt(1j*omega)
a_c = np.sqrt(-1j*omega)
xi = np.linspace(-1, 1, int(1e5))

factor = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
D_para0 = 1+0.5*sci.trapz(np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)
print(np.sqrt(2*omega/Sc))

exp_u2 = np.zeros((len(epsilon), len(kappas)))
periods = 4

for i in range(len(epsilon)):
	eps = epsilon[i]
	for j in range(len(kappas)):
		try:
			data = np.loadtxt(dirr + "Lx"+str(Lx[j])[:4]+"_tau5.0_eps"+str(eps)[:4]+"_nu1.0_D1.0_fzero0.0_fone10.0_res"+str(100*(1+eps))[:3]+"_dt0.01/tdata.dat")
			Bren = data[:, 8]
			time = data[:, 0]
		except:
			Bren = np.zeros(len(Bren))

		start_index = np.argmin(abs(time -(periods-2.35)*tau))
		end_index   = np.argmin(abs(time -(periods-0.35)*tau))
		plt.plot(time, Bren, label=r"$\kappa$="+str(kappas[j]))[:4]
		plt.plot(time[start_index:end_index], Bren[start_index:end_index], "--")
		exp_u2[i, j] = integrate.trapz(Bren[start_index:end_index], time[start_index:end_index])/(2*tau)
	plt.legend(loc="best")
	plt.show()
plt.show()

for j in range(len(epsilon)):
	plt.plot(kappas, exp_u2[j, :], "o", label=r"$\epsilon$="+str(epsilon[j])[:4])
	#plt.plot(kappas, 105*(exp_u2[j, :]-1)/(2), "-")

plt.xlabel(r"Boundary wave number $\kappa$", fontsize=14)
plt.ylabel(r"Total $D_\parallel $", fontsize=14)
plt.legend(loc="best", fontsize=12)

#filename = "figures/comparison_numeric_analytic_D.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


"""
for i in range(len(kappas)):
	plt.plot(epsilon, abs((u_squared_ana[:,i]-exp_u2[:,i])/exp_u2[:,i]), "o", label="$\kappa=$"+str(kappas[i])[:5], color=sns.color_palette()[i])
	plt.plot(epsilon, abs((u_squared_ana[:,i]-exp_u2[:,i])/exp_u2[:,i]), "--", linewidth=1, color=sns.color_palette()[i])

epsilon = np.linspace(0.01, max(epsilon), int(1e3))
plt.plot(epsilon, epsilon**4/(1-epsilon))
plt.yscale("log")
plt.xlabel(r"Boundary amplitude $\epsilon$", fontsize=14)
plt.ylabel(r"Relative difference analytic and numerical kinetic energy $\langle u^2 \rangle/2$", fontsize=14)
plt.legend(loc="best", fontsize=12)

filename = "figures/comparison_numeric_analytic_difference_D.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
plt.show()
"""