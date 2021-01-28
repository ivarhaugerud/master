import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

T = 3
omega = 2*np.pi/T
epsilon = 0.0
Pe      = 10
xi      = np.linspace(-1, 1, int(1e5))

periods = 5000
datafiles = periods*25
half_way = int(datafiles/2)
t = np.linspace(0, T*periods, datafiles)

nus = np.array([0.5, 4.0])
U = np.zeros(len(nus))
U[0] =  0.6204847210770329
U[1] = 0.03588638734002204

var    = np.zeros((len(nus), datafiles))
D_para = np.zeros((len(nus), datafiles-1))
D      = np.zeros((len(nus), 2))
D_ana = np.zeros(len(nus))

x = np.linspace(1e-4, 1, int(1e3))
gamma = np.sqrt(1j*x)

for i in range(len(nus)):
	Sc = nus[i]#/Dm
	F0 = 3/nus[i]
	gamma   = np.sqrt(1j*omega/Sc)
	gamma_c = np.conj(gamma)
	a       = np.sqrt(1j*omega)
	a_c     = np.conj(a)
	rho = a 
	rho_c = a_c

	factor  = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
	D_ana[i] = (1 + np.real(factor * 0.5 * sci.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)))
	
num_D_para = np.zeros(len(nus))

for i in range(len(nus)):
	Dm = U[i]/Pe
	var[i, :] = np.load("flow_fields/zero_eps/mu_"+str(nus[i]) +"/pos_2/var.npy")/Dm
	#plt.plot(np.loadtxt("flow_fields/zero_eps/mu_"+str(nus[i]) +"/tdata.dat")[:, 8])

	num_D_para[i] = sci.trapz(np.loadtxt("flow_fields/zero_eps/mu_"+str(nus[i]) +"/tdata.dat")[-3000:, 8], np.loadtxt("flow_fields/zero_eps/mu_"+str(nus[i]) +"/tdata.dat")[-3000:, 0])/T

plt.plot(nus, num_D_para, "o", label="numerisk")
plt.plot(nus, D_ana, "o", label="analytisk")
plt.legend(loc="best")
plt.xlabel("viskositet")
plt.ylabel("D_eff")
plt.show()

var[-1, :] /= 10
for i in range(len(nus)):
	D_para[i, :] = var[i, 1:]/t[1:]
	plt.plot(t[1:], D_para[i, :], label=str(nus[i]))
	plt.plot(t, D_ana[i]*np.ones(len(t)), label=str(nus[i]))
plt.legend(loc="best")
plt.show()


for i in range(len(nus)):
	D[i, 0] = np.mean(D_para[i, half_way:])
	D[i, 1] = np.std(D_para[i, half_way:])

D[0, 0] = np.mean(D_para[0, int(half_way):])
D[0, 1] =  np.std(D_para[0, int(half_way):])

plt.plot(nus, D_ana)
plt.errorbar(nus, D[:,0], yerr=D[:,1], fmt="o")
plt.show()
np.save("data/D_eff_vs_gamma", D)