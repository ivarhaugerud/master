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
Lx = np.array([1.25, 2.09, 4.18, 62.8, 8.97])
kappaa = 2*np.pi/Lx 
epsilon = np.arange(0.0, 0.31, 0.05)
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

for i in range(len(epsilon)):
	eps = epsilon[i]
	for j in range(len(kappas)):
		res = int(100*(1+float(eps)))
		try:
			filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(eps)[:4])+"_nu1.2_D0.3_fzero0.0_fone3.0_res"+str(res)+"_dt0.01/tdata.dat"
			tdat = np.loadtxt(dirr + filename)
		#try:
		except:
			filename = "Lx"+str(Lx[j])[:4]+"_tau"+str(tau)+"_eps"+str(str(eps)[:3])+"_nu1.2_D0.3_fzero0.0_fone3.0_res"+str(res)+"_dt0.01/tdata.dat"
			tdat = np.loadtxt(dirr + filename)
		#except:
		#	print("no file for kappa="+str(kappas[j]) + " epsilon =" + str(eps), filename)
		t = tdat[:,0]
		u2 = tdat[:,6]
		start_index = np.argmin(abs(t -(periods-2.25)*tau))
		end_index   = np.argmin(abs(t -(periods-0.25)*tau))
		plt.plot(t, u2)
		plt.plot(t[start_index:end_index], u2[start_index:end_index], "--")
		exp_u2[i, j] = integrate.trapz(u2[start_index:end_index], t[start_index:end_index])/(2*tau)
	#plt.show()
plt.show()
###
#ANALYTIC RESULTS
###


for j in range(len(epsilon)):
	plt.plot(kappas, exp_u2[j, :], "o", color=sns.color_palette()[j])
	plt.plot(kappas, exp_u2[j, :], "-", color=sns.color_palette()[j])

plt.xlabel(r"Boundary wave number $\kappa$", fontsize=14)
plt.ylabel(r"Total $D_\parallel $", fontsize=14)
plt.legend(loc="best", fontsize=12)

filename = "figures/comparison_numeric_analytic_D.pdf"
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

filename = "figures/comparison_numeric_analytic_difference_D.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
plt.show()
