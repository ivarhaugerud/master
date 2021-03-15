"""
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import os 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

tau = np.logspace(-1, 1, 7)
omega = 2*np.pi/tau
T = 500
dt  = tau/T
D   = np.zeros((len(tau)))
difference = np.zeros(np.shape(D))

for i in range(len(tau)):
	print(i, tau[i])
	try:
		filename = "data_test/vary_tau_Lx4.48/Lx4.487_tau"+ str(round(tau[i], 3)) +"_eps0.5_nu1.2_D0.1_fzero0.0_fone12.0_res150_dt" + str(round(dt[i], 6)) + "/tdata.dat"
		tdat = np.loadtxt(filename)
	except:
		filename = "data_test/vary_tau_Lx4.48/Lx4.487_tau"+ str((tau[i]))[:5] +"_eps0.5_nu1.2_D0.1_fzero0.0_fone12.0_res150_dt" + str(round(dt[i], 6)) + "/tdata.dat"
		tdat = np.loadtxt(filename)
	t = tdat[:,0]
	B = tdat[:,8]
	D[i] = sci.trapz(B[-T:], t[-T:])/(tau[i])
	difference[i] = abs(D[i]-sci.trapz(B[-2*T:-T], t[-2*T:-T])/(tau[i]))/D[i]

plt.plot(tau, difference, "o")
plt.yscale("log")
plt.legend(loc="best")
plt.xscale("log")
plt.show()

plt.plot(tau, D, "o", markersize=3)
plt.xscale("log")
plt.xlabel(r"Viscosity $[\nu]$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $[D_\parallel]$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = "figures/vary_visc.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
"""


import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import os 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

tau = np.logspace(-2, 2, 10)
Lx = np.array([2.855, 3.49])
kappa = 2*np.pi/Lx
omega = 2*np.pi/tau
T = 750
dt  = tau/T
D   = np.zeros((len(tau), len(Lx)))
difference = np.zeros(np.shape(D))

for j in range(len(Lx)):
	for i in range(len(tau)):
		print(i, tau[i])
		try:
			filename = "data_test/vary_tau/Lx"+str(Lx[j])+"_tau"+ str(round(tau[i], 3)) +"_eps0.3_nu1.3_D0.7_fzero0.0_fone12.0_res150_dt" + str(round(dt[i], 6)) + "/tdata.dat"
			tdat = np.loadtxt(filename)
		except:
			try:
				filename = "data_test/vary_tau/Lx"+str(Lx[j])+"_tau"+ str((tau[i]))[:5] +"_eps0.3_nu1.3_D0.7_fzero0.0_fone12.0_res150_dt" + str(round(dt[i], 6)) + "/tdata.dat"
				tdat = np.loadtxt(filename)
			except:
				print("no file for: ", Lx[j], tau[i])
		t = tdat[:,0]
		B = tdat[:,8]
		D[i,j] = sci.trapz(B[-T:], t[-T:])/(tau[i])
		plt.plot(t/tau[i], B)
		plt.plot(t[-T:]/tau[i], B[-T:])
		difference[i,j] = abs(D[i,j]-sci.trapz(B[-2*T:-T], t[-2*T:-T])/(tau[i]))/D[i, j]
		plt.show()

plt.plot(tau, difference, "o")
plt.yscale("log")
plt.legend(loc="best")
plt.xscale("log")
plt.show()

for i in range(len(Lx)):
	plt.plot(tau, D[:,i], "o", markersize=3, label=r"$\kappa=%3.2f$ " % (kappa[i]))

plt.xscale("log")
plt.xlabel(r"Viscosity $[\nu]$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $[D_\parallel]$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8)
plt.axis([0.007, 125, 0.88, 1.15])
#filename = "figures/vary_visc.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
