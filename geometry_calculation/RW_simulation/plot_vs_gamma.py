import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci

T = 3
omega = 2*np.pi/T
epsilon = 0.0
U_scale = 5
Pe = 20 
Dm = U_scale/Pe

periods = 5000
datafiles = periods*25
half_way = int(datafiles/2)
t = np.linspace(0, T*periods, datafiles)

nus = np.array([0.5, 4.0])

var    = np.zeros((len(nus), datafiles))
D_para = np.zeros((len(nus), datafiles-1))
D      = np.zeros((len(nus), 2))
D_ana = np.zeros(len(nus))

x = np.linspace(1e-4, 1, int(1e3))
gamma = np.sqrt(1j*x)


plt.plot(x, abs(np.real(np.tanh(np.conj(gamma))/gamma + 1j*np.tanh(gamma)/gamma)))
plt.plot(x, x*(2/3-34*x*x/315))
plt.plot(x, x*x*x*x*x*x, "--")
#plt.yscale("log")
#plt.xscale("log")
plt.show()

for i in range(len(nus)):
	Sc = nus[i]#/Dm
	F0 = 3/nus[i]
	gamma   = np.sqrt(1j*omega/Sc)
	gamma_c = np.conj(gamma)
	a       = np.sqrt(1j*omega)
	a_c     = np.conj(a)
	rho = a 
	rho_c = a_c
	xi      = np.linspace(-1, 1, int(1e5))

	factor  = Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(omega*omega*omega*(Sc-1)*(Sc-1))
	#print("num: ", 0.5 * sci.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi))
	#print("ana: ",  (1/(rho*rho+gamma*gamma))*(-rho/np.tanh(rho)+rho_c/np.tanh(rho_c)+gamma_c/np.tanh(gamma_c) - gamma/np.tanh(gamma))  + (rho/np.tanh(rho) - rho_c/np.tanh(rho_c))/(2*rho*rho) + (gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c))/(2*gamma*gamma) )
	#print(  ((1-Sc)/(2*rho*rho*(1+Sc)))*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)) - ((1-Sc)/(2*gamma*gamma*(1+Sc)))*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c))   )
	D_ana[i] = (1 + np.real(factor * 0.5 * sci.trapz( np.sinh(a*xi)*np.sinh(a_c*xi)/(np.sinh(a)*np.sinh(a_c)) + np.sinh(gamma*xi)*np.sinh(gamma_c*xi)/(np.sinh(gamma)*np.sinh(gamma_c)) - np.sinh(a*xi)*np.sinh(gamma_c*xi)/(np.sinh(a)*np.sinh(gamma_c)) - np.sinh(gamma*xi)*np.sinh(a_c*xi)/(np.sinh(gamma)*np.sinh(a_c)), xi)))
	print(rho_c/rho)
	print(D_ana[i])
	print(1+Sc*Sc*Sc*Pe*Pe*F0*F0/(2*omega*omega*omega*(Sc*Sc-1))*(np.tanh(gamma_c)/gamma + 1j*np.tanh(gamma)/gamma - 1j*np.tanh(gamma)*np.tanh(gamma_c)/(rho*np.tanh(rho_c)) - np.tanh(gamma)*np.tanh(gamma_c)/(rho*np.tanh(rho))))
for i in range(len(nus)):
	var[i, :] = np.load("flow_fields/zero_eps/mu_"+str(nus[i]) +"/pos_2/var.npy")/Dm
"""

for i in range(len(nus)):
	plt.plot(t, var[i, :], label=str(nus[i]))

plt.legend(loc="best")	
plt.show()
"""
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