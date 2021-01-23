import numpy as np 
import matplotlib.pyplot as plt 
from scipy import integrate
from scipy import integrate as sci

tau = 3
dt = 0.006
timesteps = int(tau/(dt))
periods = 4000
datafiles = periods*25
half_way = int(datafiles/2)
t = np.linspace(0, tau*periods, datafiles)
periods_of_flow = 8/3
xi = np.linspace(-1, 1, int(1e5))
Pe = 6.0
Sc = 3.6
F0 = 12.0/Sc
omega = 2*np.pi/tau
gamma   = np.sqrt(1j*omega/Sc)
gamma_c = np.conj(gamma)
a       = np.sqrt(1j*omega)
a_c     = np.conj(a)
rho = a 
rho_c = a_c

factor  = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
D_ana = np.ones(len(t))*(1 + np.real(factor * 0.5 * sci.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)))


kappas = np.array([0.4, 0.5, 0.66, 1])
exp_u2 = np.zeros(len(kappas))
Lxs    = 2*np.pi/kappas
files = ["flow_fields/non_zero_eps/Lx15.71_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/",
		 "flow_fields/non_zero_eps/Lx12.56_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/",
		 "flow_fields/non_zero_eps/Lx9.42_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/",
		 "flow_fields/non_zero_eps/Lx6.28_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/"]

var    = np.zeros((len(kappas), datafiles))
D_para = np.zeros((len(kappas), datafiles-1))
D      = np.zeros((len(kappas), 2))

for i in range(len(kappas)):
	var[i, :] = np.load(files[i]+"pos/var.npy")
	tdat = np.loadtxt(files[i] +"tdata.dat")
	time = tdat[:,0]
	u2   = tdat[:,4]
	exp_u2[i] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/tau


for i in range(len(kappas)):
	plt.plot(t, var[i, :], label=files[i])

plt.legend(loc="best")	
plt.show()
Pe = 6

for i in range(len(kappas)):
	Dm = exp_u2[i]/Pe
	D_para[i, :] = var[i, 1:]/(Dm*t[1:])
	plt.plot(t[1:], D_para[i, :], label=str(kappas[i]))
plt.plot(t, D_ana)

plt.legend(loc="best")
plt.show()

print(np.real(np.sqrt(2*1j*omega/Sc)))
for i in range(len(kappas)):
	D[i, 0] = np.mean(D_para[i, half_way:])
	D[i, 1] = np.std(D_para[i, half_way:])

D[0, 0] = np.mean(D_para[0, int(3*half_way/2):])
D[0, 1] =  np.std(D_para[0, int(3*half_way/2):])

plt.errorbar(kappas, D[:,0], yerr=D[:,1], fmt="o")
plt.show()
np.save("data/D_eff_vs_kappa", D)

eps = 0.25
plt.errorbar(kappas, (D[:,0]-D_ana[0])/(eps*eps), yerr=D[:,1]/(eps*eps), fmt="o")
plt.show()
np.save("data/D_eff_vs_kappa", D)