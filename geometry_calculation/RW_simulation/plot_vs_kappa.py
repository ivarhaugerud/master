import numpy as np 
import matplotlib.pyplot as plt 
from scipy import integrate
from scipy import integrate as sci

tau = 3
dt = 0.006
timesteps = int(tau/(dt))
periods = 10000
datafiles = periods*20
half_way = int(datafiles/2)
skip = int(periods*timesteps/datafiles)

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

#geometry parameters
Lx =  np.array([1.05, 2.09, 6.28, 9.42, 12.56, 15.71, 25.13])
dirr = []

for i in range(len(Lx)):
    dirr.append("flow_fields/non_zero_eps/Lx"+str(Lx[i])+"_tau3.0_eps0.25_nu3.6_D1.0_fzero0.0_fone12.0_res100_dt0.006/")

kappas = 2*np.pi/Lx
exp_u2 = np.zeros(len(kappas))

var    = np.zeros((len(kappas), datafiles-1))
D_para = np.zeros((len(kappas), datafiles-1))
D      = np.zeros((len(kappas), 2))
t = t[1:]

Pe = 1

for i in range(len(kappas)):
	print(np.shape(np.load(dirr[i]+"var.npy")))
	var[i, :] = np.load(dirr[i]+"var.npy")[1:]
	var[i, :] -= var[i, 0]
	tdat = np.loadtxt(dirr[i] +"tdata.dat")[1:]
	time = tdat[:,0]
	u2   = tdat[:,4]
	exp_u2[i] = integrate.trapz(u2[-timesteps:], time[-timesteps:])/tau
	Dm = np.sqrt(exp_u2[i])/Pe
	D_para[i, :] = var[i, :]/(Dm*t)
	plt.plot(t, D_para[i, :], label=str(kappas[i])[:4])

plt.legend(loc="best")
plt.show()

for i in range(len(kappas)):
	plt.plot(t, var[i, :], label=str(Lx[i]))

plt.legend(loc="best")	
plt.show()

print(np.real(np.sqrt(2*1j*omega/Sc)))
for i in range(len(kappas)):
	D[i, 0] = np.mean(D_para[i, half_way:])
	D[i, 1] = np.std(D_para[i, half_way:])

plt.errorbar(kappas, D[:,0], yerr=D[:,1], fmt="o")
plt.show()
np.save("data/D_eff_vs_kappa", D)

eps = 0.25
plt.errorbar(kappas, (D[:,0]-D_ana[0]), yerr=D[:,1], fmt="o")
plt.show()
np.save("data/D_eff_vs_kappa", D)


### PLOT POSITIONS
datafiles = 2000
x = np.zeros(datafiles)
y = np.zeros(datafiles)

for l in range(len(Lx)):
	for i in range(datafiles):
		pos = np.load(dirr[l]+"pos/RW_positions_"+str(int(i*skip))+".npy")
		x[i] = pos[0, 7]
		y[i] = pos[1, 7]
	plt.title("Lx = "+str(Lx[l]))
	plt.plot(x, y)
	plt.show()

