import numpy as np 
from scipy import integrate
import scipy.integrate as sci
import matplotlib.pyplot as plt 

data1 = np.load("data_test/newest_visc_1.5_var_over_t.npy")
data2 = np.load("data_test/newest_visc_3.0_var_over_t.npy")
data3 = np.load("data_test/newest_visc_5.0_var_over_t.npy")

visc = np.array([1.5, 3.0, 5.0])
analytic = np.zeros(len(visc))
U        = np.zeros(len(visc))
U2       = np.zeros(len(visc))
tau = 3.0
omega = 2*np.pi/tau
Pe_sim = 10
dt           = 0.006
tau          = 3.0 
timesteps    = int(tau/(dt))
dirr = []
for i in range(len(visc)):
    dirr.append("RW_simulation/flow_fields/zero_eps/Lx12.56_tau3.0_eps0.0_nu"+str(visc[i])+"_D1.0_fzero0.0_fone12.0_res100_dt0.006")

colors = ["r", "b", "g"]
for i in range(len(visc)):
    Sc = visc[i]
    F0 = 12/visc[i]

    tdat = np.loadtxt(dirr[i] +"/tdata.dat")
    U[i]  = integrate.trapz(tdat[-timesteps:, 4], tdat[-timesteps:, 0])/tau
    plt.plot(tdat[:, 0], tdat[:, 4])
    plt.plot(tdat[-timesteps:, 0], tdat[-timesteps:, 4])
    plt.show()
    D = np.sqrt(U[i])/Pe_sim
    gamma   = np.sqrt(1j*omega/Sc)
    gamma_c = np.conj(gamma)
    rho = np.sqrt(1j*omega/D)
    rho_c = np.conj(rho)
    Pe = 1/D
    #print(Pe, D, U)
    analytic[i] = 1 + np.real(Sc*Sc*Sc*Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*omega*omega*omega*(Sc*Sc-1))*(1/(gamma*np.tanh(gamma)) + 1j/(gamma*np.tanh(gamma_c)) - 1/(rho*np.tanh(rho)) - 1j/(rho*np.tanh(rho_c))))
    D2          = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
    analytic[i] = D2
    print(analytic[i], D2)

plt.plot(data1, color=colors[0], label=r"$\nu=1.5$")
plt.plot(data2, color=colors[1], label=r"$\nu=3.0$")
plt.plot(data3, color=colors[2], label=r"$\nu=5.0$")

plt.plot(np.ones(len(data1))*analytic[0], "--", color=colors[0])
plt.plot(np.ones(len(data1))*analytic[1], "--", color=colors[1])
plt.plot(np.ones(len(data1))*analytic[2], "--", color=colors[2])
plt.legend(loc="best")
plt.xlabel(r"datapoints")
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$")
plt.savefig("figures/comparison_analytic_and_numeric.png")
plt.show()



D = np.zeros(len(visc))
cutoff = int(0.8*len(data1))
D[0] = (np.mean(data1[cutoff:]))
D[1] = (np.mean(data2[cutoff:]))
D[2] = (np.mean(data3[cutoff:]))

plt.plot(visc, D, "o")
plt.plot(visc, analytic, "o")
plt.show()

plt.plot(visc, ((D)/(analytic)) , "o")
plt.show()