import os
import numpy as np 
from scipy import integrate
import scipy.integrate as sci
import matplotlib.pyplot as plt 

plt.style.use(['science','no-latex', 'grid'])
#sns.color_palette("hls", 1)
import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

visc = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
data = {}

for i in range(len(visc)):
    data[str(i)] = np.load("data_test/D0.25__visc_"+str(visc[i])+"_var_over_t.npy")


analytic = np.zeros(len(visc))
tau = 3.0
omega = 2*np.pi/tau
tau          = 3.0 
D = 0.25
periods = 1000
t = np.linspace(0, periods, int(len(data["1"])))

for i in range(len(visc)):
    Sc = visc[i]
    F0 = 12/visc[i]
    gamma   = np.sqrt(1j*omega/Sc)
    gamma_c = np.conj(gamma)
    rho = np.sqrt(1j*omega/D)
    rho_c = np.conj(rho)
    Pe = 1/D
    analytic[i] = np.real(1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c))) )

D = np.zeros(len(visc))
error = np.zeros(len(visc))
cutoff = int(0.5*len(data["1"]))
plot_skip = 100

for i in range(len(visc)):
    D[i] = (np.mean(data[str(i)][cutoff:]))
    error = (np.std(data[str(i)][cutoff:]))

for i in range(len(visc)-4):
    fig = plt.figure(1)
    plt.plot(t[::plot_skip], data[str(i)][::plot_skip], "C"+str(i), label=r"$\nu=$"+str(visc[i]))
    plt.plot(t[::plot_skip], np.ones(len(t[::plot_skip]))*analytic[i], "--", color="C"+str(i))

plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Time [periods $T$]", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.savefig("figures/comparison_analytic_and_numeric.png")
name = "figures/D_eff_vs_t.pdf"
fig.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))

fig = plt.figure(2)
plt.errorbar(visc, D, yerr=error, fmt="o", markersize=3 , label="Random Walk")


visc = np.linspace(0.4, 5, 60)
analytic = np.zeros(len(visc))

for i in range(len(visc)):
    D = 0.25
    Sc = visc[i]
    F0 = 12/visc[i]
    gamma   = np.sqrt(1j*omega/Sc)
    gamma_c = np.conj(gamma)
    rho = np.sqrt(1j*omega/D)
    rho_c = np.conj(rho)
    Pe = 1/D
    analytic[i] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))

plt.plot(visc, analytic, label="analytic")
plt.legend(loc="best", fontsize=8)
plt.xlabel(r"Viscosity $\nu$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $D_\parallel$", fontsize=8)
plt.savefig("figures/comparison_analytic_and_numeric.png")
plt.tick_params(axis='both', which='major', labelsize=8)
name = "figures/D_eff_vs_visc.pdf"
fig.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()
