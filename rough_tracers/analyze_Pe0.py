import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import os
import scipy.interpolate as sci 

plt.style.use("bmh")
sns.color_palette("hls", 1)

def func2(x, a, b):
    return a*(1 + b/np.sqrt(x))

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

Brenner = np.loadtxt("data_square/Pe0_300.dat")
b = Brenner[:, 0]
D_eff_B = Brenner[:, 1]

interpool = sci.interp1d(b, D_eff_B)

data_strings = ["data_square/RandomWalkers/Re0.000000_b0.000000_Dm1.000000_U0.000000_dt0.000010_Nrw5000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.200000_Dm1.000000_U0.000000_dt0.000001_Nrw10000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.400000_Dm1.000000_U0.000000_dt0.000003_Nrw5000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.600000_Dm1.000000_U0.000000_dt0.000003_Nrw10000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.800000_Dm1.000000_U0.000000_dt0.000003_Nrw8000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.000000_Dm1.000000_U0.000000_dt0.000003_Nrw5000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.200000_Dm1.000000_U0.000000_dt0.000002_Nrw500/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.400000_Dm1.000000_U0.000000_dt0.000002_Nrw500/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.600000_Dm1.000000_U0.000000_dt0.000002_Nrw700/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.800000_Dm1.000000_U0.000000_dt0.0000005_Nrw500/tdata.dat"]

B     = np.linspace(0, 1.8, 10)
B[0] += 1e-5
D_eff = np.zeros((len(B), 2))
t_cuts = [20, 2, 7.5, 12.5, 12.5, 19.4, 36.43, 47.12, 60, 76.5]

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    t = data[:, 0]
    variance = data[:, 2]
    t_cut = np.argmin(abs(t-t_cuts[i]))

    Ds = variance[1:]/(2*t[1:])
    t_prime = t[1:]
    """
    plt.plot(t_prime, Ds)
    plt.plot(t_prime, interpool(B[i])*np.ones(len(t_prime)))

    plt.title(r"Roughness $b=$" + str(B[i]))
    plt.ylabel(r"Effective diffusion coefficient $D_{\parallel}$ $[D_m]$", fontsize=14)
    plt.xlabel(r"Timesteps $N \times 10^{6}$", fontsize=14)
    plt.legend(loc="best", fontsize=12)
    plt.savefig("../../oppgave/figures/Pe0_r"+str(B[i])+".pdf")

    plt.show()
    """
	
    D_eff[i, 0] = np.mean(variance[t_cut:]/(2*t[t_cut:]))
    #D_eff[i, 1] = np.max(abs(variance[t_cut:]/(2*t[t_cut:])-D_eff[i,0]))
    D_eff[i, 1] = np.std((variance[t_cut:]/(2*t[t_cut:])))

plt.errorbar(B, D_eff[:, 0], yerr=D_eff[:,1], fmt="o", label="Random Walks")
plt.plot(b, D_eff_B, label=r"Brenner")
b_test = np.linspace(0, 2-1e-5, 1e4)
prop_open = 2*b_test/(2+b_test)
plt.plot(b_test, (1-prop_open)*(1+prop_open), label="Reguera")
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/Pe0.pdf")
plt.show()

plt.errorbar(B, (D_eff[:, 0]-interpool(B))/interpool(B), yerr=D_eff[:,1]/interpool(B), fmt='o')
plt.plot(B, np.zeros(len(B)), "k")
plt.ylabel(r"Relative difference Brenner and simulated $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/Pe0_differenece.pdf")
plt.show()



data_strings = ["data_square/RandomWalkers/Re0.000000_b0.600000_Dm1.000000_U0.000000_dt0.000003_Nrw10000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.800000_Dm1.000000_U0.000000_dt0.000003_Nrw8000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.000000_Dm1.000000_U0.000000_dt0.000003_Nrw5000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.200000_Dm1.000000_U0.000000_dt0.000002_Nrw500/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.400000_Dm1.000000_U0.000000_dt0.000002_Nrw500/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.600000_Dm1.000000_U0.000000_dt0.000002_Nrw700/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.800000_Dm1.000000_U0.000000_dt0.000002_Nrw1000/tdata.dat"]

B     = np.linspace(0.6, 1.8, 7)
D_eff = np.zeros((len(B), 2))

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    t = data[:, 0]
    variance = data[:, 2]

    Ds = variance[1:]/(2*t[1:])
    t_prime = t[1:]

    plt.plot(t_prime/(B[i]**2), (Ds-interpool(B[i]))/B[i], label=r"$b=$"+str(B[i]))

plt.xlabel(r"Time $t$ [$1/(D_m b^2)$]")
plt.ylabel(r"Relative difference RW-Brenner $D_{\parallel}$ [$D_m b$]")
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/data_collapse.pdf")
plt.show()