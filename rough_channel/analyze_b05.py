import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import os
import scipy.interpolate as sci 

plt.style.use("bmh")
sns.color_palette("hls", 1)

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def linear_regresion(x, y):
    n = float(len(x))
    D = float(np.sum(np.square(x)) - (np.sum(x)**2)/n)
    E = float(np.sum(x*y) - np.sum(x)*np.sum(y)/n)
    F = float(np.sum(np.square(y)) - (np.sum(y)**2)/n)

    delta_m = np.sqrt((1/(n-2))*(D*F-E**2)/(D**2))
    delta_c = np.sqrt(1/(n-2)*(D/n+np.mean(x)**2)*(D*F-E**2)/(D**2))
    m = E/D
    c = np.mean(y)-m*np.mean(x)

    return m, c, delta_m, delta_c
    #using linear regression from Squires, with uncertainty to find slope and constant term

Brenner = np.loadtxt("data_square/Re0.0_b0.5_res300.dat")
Pe = Brenner[:, 0]
D_eff_B = Brenner[:, 1]

interpool = sci.interp1d(Pe, D_eff_B)

data_strings = ["data_square/RandomWalkers/Re0.000000_b0.500000_Dm1.000000_U1.000000_dt0.000015_Nrw3000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.500000_Dm0.250000_U1.000000_dt0.000060_Nrw2000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.500000_Dm0.015625_U1.000000_dt0.000960_Nrw2000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.500000_Dm0.003906_U1.000000_dt0.003840_Nrw2000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.500000_Dm0.000977_U1.000000_dt0.015360_Nrw2000/tdata.dat"]

Dm     = [1, 0.25, 0.0625, 0.015625, 0.003906, 0.000977]
D_eff = np.zeros((len(Dm), 2))
t_cuts = [6.4, 32.3, 800, 460, 363]
Pe2 = 1/np.array(Dm)

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    t = data[:, 0]
    variance = data[:, 2]
    t_cut = np.argmin(abs(t-t_cuts[i]))

    Ds = variance[1:]/(Dm[i]*2*t[1:])
    t_prime = t[1:]

    plt.plot(t_prime, Ds)
    plt.plot(t_prime, interpool(Pe2[i])*np.ones(len(t_prime)))
    plt.plot(t_prime[t_cut:], interpool(Pe2[i])*np.ones(len(t_prime[t_cut:])))

    plt.title(r"Diffusion $D_m=$" + str(Dm[i]))
    plt.ylabel(r"Effective diffusion coefficient $D_{\parallel}$ $[D_m]$", fontsize=14)
    plt.xlabel(r"Timesteps $N \times 10^{6}$", fontsize=14)
    plt.legend(loc="best", fontsize=12)
    #plt.savefig("../../oppgave/figures/b_0_Dm"+str(Dm[i])+".pdf")

    plt.show()
	
    D_eff[i, 0] = np.mean(variance[t_cut:]/(Dm[i]*2*t[t_cut:]))
    #D_eff[i, 1] = np.max(abs(variance[t_cut:]/(2*t[t_cut:])-D_eff[i,0]))
    D_eff[i, 1] = np.std((variance[t_cut:]/(2*t[t_cut:])))

Pe = np.linspace(1, 1e3, 1e5)
analytic = interpool(Pe)

plt.errorbar(Pe2, D_eff[:, 0], yerr=D_eff[:,1], fmt="o", label="Random Walks")
plt.plot(Pe, analytic, label=r"Brenner")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"Peclet number Pe", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/b05.pdf")
plt.show()

analytic = interpool(Pe2)

plt.errorbar(Pe2, (D_eff[:, 0]-analytic)/analytic, yerr=D_eff[:,1]/analytic, fmt='o')
plt.plot(Pe2, np.zeros(len(Pe2)), "k")
plt.ylabel(r"Relative difference Brenner and simulated $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/Pe0_differenece.pdf")
plt.show()


data_strings = ["data_square/RandomWalkers/Re0.000000_b0.000000_Dm1.000000_U1.000000_dt0.000003_Nrw1000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.200000_U1.000000_dt0.000100_Nrw3000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.040000_U1.000000_dt0.000500_Nrw3000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.008000_U1.000000_dt0.002500_Nrw3000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.001600_U1.000000_dt0.012500_Nrw3000/tdata.dat"]

Dm     = [1, 0.2, 0.04, 0.008, 0.0016]
Pe = 1/np.array(Dm)
D_para = 1+2*Pe*Pe/105

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    t = data[:, 0]
    variance = data[:, 2]

    Ds = variance[1:]/(2*t[1:]*Dm[i])
    t_prime = t[1:]

    plt.plot(t_prime*(Dm[i]), (Ds-D_para[i])/D_para[i], label=r"$Dm=$"+str(Dm[i]))

plt.xlabel(r"Time $t$ [$1/(D_m b^2)$]")
plt.ylabel(r"Relative difference RW-Brenner $D_{\parallel}$ [$D_m b$]")
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/data_collapse.pdf")
plt.show()
