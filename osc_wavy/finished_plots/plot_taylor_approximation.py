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
root = "../../../master_latex/results/"

base = "../data_test/vary_nu/"
base = "../data_test/taylor_approximation_try2/taylor_approximation_try2/"
tau = np.logspace(-1, 2, 5)

T    = 500
kappa = 0.2
dt   = tau/T
D  = np.zeros((len(tau)))
U  = np.zeros(np.shape(D))
difference = np.zeros(np.shape(D))
D_ana = np.zeros(len(tau))

nu = 8000
epsilon = 0.2
Dm = 0.1
F0 = 96000/nu 
omega = 2*np.pi/tau 
RHO = np.sqrt(1j*omega/Dm)
Pe = 1/Dm
from numpy import *
D0 = np.zeros(len(tau))

for i in range(len(tau)):
	rho = RHO[i]
	gamma = np.sqrt(1j*omega[i]/nu)
	D_ana[i] = 0.5*F0**2*Pe**2*(0.583333333333333*rho**8*conjugate(rho)**2 + 0.5*rho**8 + 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho)) + 1.75*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 1.75*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 0.5*rho**7*conjugate(rho)**4/tanh(rho) + 0.5*rho**7*conjugate(rho)**4/tanh(rho)**3 - 2.75*rho**6*conjugate(rho)**4 - 1.0*rho**6*conjugate(rho)**2 - 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4.0*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4.0*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 2.25*rho**6*conjugate(rho)**4/tanh(rho)**2 + 1.0*rho**5*conjugate(rho)**6/tanh(rho) - 0.25*rho**5*conjugate(rho)**4/tanh(rho) - 1.0*rho**5*conjugate(rho)**6/tanh(rho)**3 + 2.75*rho**4*conjugate(rho)**6 + 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 0.25*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 2.25*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4.0*rho**4*conjugate(rho)**6/tanh(rho)**2 - 0.5*rho**3*conjugate(rho)**8/tanh(rho) + 4.0*rho**3*conjugate(rho)**6/tanh(rho) + 0.5*rho**3*conjugate(rho)**8/tanh(rho)**3 - 0.583333333333333*rho**2*conjugate(rho)**8 + 1.0*rho**2*conjugate(rho)**6 + 1.75*rho**2*conjugate(rho)**8/tanh(rho)**2 - 1.75*rho*conjugate(rho)**8/tanh(rho) - 0.5*conjugate(rho)**8)/(rho**4*(1.0*rho**6 - 3.0*rho**4*conjugate(rho)**2 + 3.0*rho**2*conjugate(rho)**4 - 1.0*conjugate(rho)**6)*conjugate(rho)**4)
	D_ana[i] = 1/2*F0**2*Pe**2*(7/12*rho**8*conjugate(rho)**2 + rho**8/2 + rho**8*conjugate(rho)**3/(2*tanh(conjugate(rho))) + 7/4*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 7/4*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1/2*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 1/2*rho**7*conjugate(rho)**4/tanh(rho) + 1/2*rho**7*conjugate(rho)**4/tanh(rho)**3 - 11/4*rho**6*conjugate(rho)**4 - rho**6*conjugate(rho)**2 - rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 9/4*rho**6*conjugate(rho)**4/tanh(rho)**2 + rho**5*conjugate(rho)**6/tanh(rho) - 1/4*rho**5*conjugate(rho)**4/tanh(rho) - rho**5*conjugate(rho)**6/tanh(rho)**3 + 11/4*rho**4*conjugate(rho)**6 + 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 1/4*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 9/4*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4*rho**4*conjugate(rho)**6/tanh(rho)**2 - 1/2*rho**3*conjugate(rho)**8/tanh(rho) + 4*rho**3*conjugate(rho)**6/tanh(rho) + 1/2*rho**3*conjugate(rho)**8/tanh(rho)**3 - 7/12*rho**2*conjugate(rho)**8 + rho**2*conjugate(rho)**6 + 7/4*rho**2*conjugate(rho)**8/tanh(rho)**2 - 7/4*rho*conjugate(rho)**8/tanh(rho) - 1/2*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4)

	D0[i]    = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))
	#D_ana[i] -= -1/(kappa*kappa) + 1/(kappa*tanh(kappa)) + 1
	try:
		data = np.loadtxt(base+"Lx31.41_tau" + str(tau[i])[:4] + "_eps0.2_nu8000.0_D0.1_fzero0.0_fone96000.0_res150_dt"+str(dt[i])[:6] +"/tdata.dat")
		D[i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau[i]
		U[i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau[i]
		difference[i] = abs(D[i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau[i])/D[i]
		#plt.plot(np.trim_zeros(data[:, 0])/tau[i], np.trim_zeros(data[:, 8]))
		#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[i], np.trim_zeros(data[:, 8])[-T:])
		#plt.show()
	
	except:
		data = np.loadtxt(base+"Lx31.41_tau" + str(tau[i])[:5] + "_eps0.2_nu8000.0_D0.1_fzero0.0_fone96000.0_res150_dt"+str(dt[i])[:6] +"/tdata.dat")
		D[i] = sci.trapz(  data[-T:, 8],  data[-T:, 0] )/tau[i]
		U[i] = sci.trapz(  data[-T:, 4],  data[-T:, 0] )/tau[i]
		difference[i] = abs(D[i] - sci.trapz(  data[-2*T:-T, 8],  data[-2*T:-T, 0] )/tau[i])/D[i]
		#plt.plot(np.trim_zeros(data[:, 0])/tau[i], np.trim_zeros(data[:, 8]))
		#plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau[i], np.trim_zeros(data[:, 8])[-T:])
		#plt.show()	


omega_contin = np.logspace(np.log10(min(2*np.pi/tau)), np.log10(max(2*np.pi/tau)))
D0_contin    = np.zeros(len(omega_contin))
D_ana_contin = np.zeros(len(omega_contin))

plt.plot(np.sqrt(2*np.pi/(tau*nu)), difference)
plt.xscale("log")
plt.yscale("log")
plt.show()
for i in range(len(omega_contin)):
	rho = np.sqrt(1j*omega_contin[i]/Dm)
	gamma = np.sqrt(1j*omega_contin[i]/nu)
	D_ana_contin[i] = 0.5*F0**2*Pe**2*(0.583333333333333*rho**8*conjugate(rho)**2 + 0.5*rho**8 + 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho)) + 1.75*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 1.75*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 0.5*rho**7*conjugate(rho)**4/tanh(rho) + 0.5*rho**7*conjugate(rho)**4/tanh(rho)**3 - 2.75*rho**6*conjugate(rho)**4 - 1.0*rho**6*conjugate(rho)**2 - 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4.0*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4.0*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 2.25*rho**6*conjugate(rho)**4/tanh(rho)**2 + 1.0*rho**5*conjugate(rho)**6/tanh(rho) - 0.25*rho**5*conjugate(rho)**4/tanh(rho) - 1.0*rho**5*conjugate(rho)**6/tanh(rho)**3 + 2.75*rho**4*conjugate(rho)**6 + 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 0.25*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 2.25*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4.0*rho**4*conjugate(rho)**6/tanh(rho)**2 - 0.5*rho**3*conjugate(rho)**8/tanh(rho) + 4.0*rho**3*conjugate(rho)**6/tanh(rho) + 0.5*rho**3*conjugate(rho)**8/tanh(rho)**3 - 0.583333333333333*rho**2*conjugate(rho)**8 + 1.0*rho**2*conjugate(rho)**6 + 1.75*rho**2*conjugate(rho)**8/tanh(rho)**2 - 1.75*rho*conjugate(rho)**8/tanh(rho) - 0.5*conjugate(rho)**8)/(rho**4*(1.0*rho**6 - 3.0*rho**4*conjugate(rho)**2 + 3.0*rho**2*conjugate(rho)**4 - 1.0*conjugate(rho)**6)*conjugate(rho)**4)
	D0_contin[i]    = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))

D_tot_ana = D0 + epsilon*epsilon*D_ana
D_tot_ana_contin = D0_contin + epsilon*epsilon*D_ana_contin


tau_num   = np.logspace(-1, 2, 5)
omega_num = 2*np.pi/tau
D_num = np.load("../finite_element/data/vary_omega.npy")

plt.figure(2)
plt.plot((omega/nu), D, "ko", markersize=3, label="Full numerical")
plt.plot((omega/nu),  D0 + epsilon*epsilon*(D_num), "ro", markersize=3, label="Approximate numerical")
plt.plot((omega_contin/nu), D_tot_ana_contin , label="Analytical")
plt.xscale("log")
plt.xlabel(r" Womersley number $\frac{\omega a^2}{\nu}$", fontsize=8)
plt.ylabel(r" Effective Diffusion coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xscale("log")
filename = root + "figures/D_eff_taylor_approx.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()



plt.figure(2)
plt.ylabel(r"$D_{eff}$")
plt.xlabel(r"viscosity $\nu$")
plt.plot((omega/nu), np.ones(len(tau))*kappa**2,   label=r"$\mathcal{O}(\kappa^2)$")
plt.plot((omega/nu), np.ones(len(tau))*epsilon**4, label=r"$\mathcal{O}(\epsilon^4)$")

plt.plot((omega/nu), abs(D-(D0 + epsilon*epsilon*D_num))/D, "o",color="C2", markersize=3, label="Numerical")
plt.plot((omega/nu), abs(D-(D0 + epsilon*epsilon*D_num))/D, "-", color="C2", linewidth=0.5)
plt.plot((omega/nu), abs(D_tot_ana-D)/D, "ko", markersize=3, color="C3", label="Analytic")
plt.plot((omega/nu), abs(D_tot_ana-D)/D, "-", linewidth=0.5, color="C3")

plt.xscale("log")
plt.xlabel(r" Womersley number $\frac{\omega a^2}{\nu}$", fontsize=8)
plt.ylabel(r" Releative difference $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xscale("log")
filename = root + "figures/D_eff_taylor_reldiff.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()




#D_ana[i] = 1/2*F0**2*Pe**2*(7/12*rho**8*conjugate(rho)**2 + rho**8/2 + rho**8*conjugate(rho)**3/(2*tanh(conjugate(rho))) + 7/4*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 7/4*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1/2*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 1/2*rho**7*conjugate(rho)**4/tanh(rho) + 1/2*rho**7*conjugate(rho)**4/tanh(rho)**3 - 11/4*rho**6*conjugate(rho)**4 - rho**6*conjugate(rho)**2 - rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 9/4*rho**6*conjugate(rho)**4/tanh(rho)**2 + rho**5*conjugate(rho)**6/tanh(rho) - 1/4*rho**5*conjugate(rho)**4/tanh(rho) - rho**5*conjugate(rho)**6/tanh(rho)**3 + 11/4*rho**4*conjugate(rho)**6 + 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 1/4*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 9/4*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4*rho**4*conjugate(rho)**6/tanh(rho)**2 - 1/2*rho**3*conjugate(rho)**8/tanh(rho) + 4*rho**3*conjugate(rho)**6/tanh(rho) + 1/2*rho**3*conjugate(rho)**8/tanh(rho)**3 - 7/12*rho**2*conjugate(rho)**8 + rho**2*conjugate(rho)**6 + 7/4*rho**2*conjugate(rho)**8/tanh(rho)**2 - 7/4*rho*conjugate(rho)**8/tanh(rho) - 1/2*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4)
