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

Brenner = np.loadtxt("data_square/Re0.0_b1.5_res200.dat")
Pe = Brenner[:, 0]
D_eff_B = Brenner[:, 1]

interpool = sci.interp1d(Pe, D_eff_B)

data_strings = ["data_square/RandomWalkers/Re0.000000_b1.500000_Dm1.000000_U1.000000_dt0.000002_Nrw1000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.500000_Dm0.250000_U1.000000_dt0.000008_Nrw1000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.500000_Dm0.062500_U1.000000_dt0.000032_Nrw1000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.500000_Dm0.015625_U1.000000_dt0.000128_Nrw1000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.500000_Dm0.004000_U1.000000_dt0.000512_Nrw1000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b1.500000_Dm0.001000_U1.000000_dt0.000128_Nrw1000/tdata.dat"]

Dm     = [1, 0.25, 0.0625, 0.015625, 0.004, 0.001]
D_eff = np.zeros((len(Dm), 2))
t_cuts = [203.4, 736, 800, 460, 363, 1881]
Pe2 = 1/np.array(Dm)

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    t = data[:, 0]
    variance = data[:, 2]
    t_cut = np.argmin(abs(t-t_cuts[i]))

    Ds = variance[1:]/(Dm[i]*2*t[1:])
    t_prime = t[1:]

    D_eff[i, 0] = np.mean(variance[t_cut:]/(Dm[i]*2*t[t_cut:]))
    D_eff[i, 1] = np.std((variance[t_cut:]/(2*t[t_cut:])))

Pe = np.linspace(1, 1e3, 1e5)
analytic = interpool(Pe)

c=sns.color_palette()[2]
plt.errorbar(Pe2, D_eff[:, 0], yerr=D_eff[:,1], fmt="o", label=r"RW - $b=1.5$", color=c)
c=sns.color_palette()[0]
plt.plot(Pe, analytic, label=r"Brenner - $b=1.5$", color=c)

Brenner = np.loadtxt("data_square/Pe0_300.dat")
b = Brenner[:, 0]
D_eff_B = Brenner[:, 1]

interpool = sci.interp1d(b, D_eff_B)

data_strings = ["data_square/RandomWalkers/Re0.000000_b0.000000_Dm1.000000_U1.000000_dt0.000003_Nrw1000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.200000_U1.000000_dt0.000100_Nrw3000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.040000_U1.000000_dt0.000500_Nrw3000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.008000_U1.000000_dt0.002500_Nrw3000/tdata.dat",
                "data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.001600_U1.000000_dt0.012500_Nrw3000/tdata.dat"]

Dm     = [1, 0.2, 0.04, 0.008, 0.0016]
D_eff = np.zeros((len(Dm), 2))
t_cuts = [34.1, 271, 673, 4608, 21912]
Pe2 = 1/np.array(Dm)

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    t = data[:, 0]
    variance = data[:, 2]
    t_cut = np.argmin(abs(t-t_cuts[i]))

    Ds = variance[1:]/(Dm[i]*2*t[1:])
    t_prime = t[1:]
    
    D_eff[i, 0] = np.mean(variance[t_cut:]/(Dm[i]*2*t[t_cut:]))
    D_eff[i, 1] = np.std((variance[t_cut:]/(2*t[t_cut:])))

Pe = np.linspace(1, 1e3, 1e5)
D_para = 1+2*Pe*Pe/105

c=sns.color_palette()[3]
plt.errorbar(Pe2, D_eff[:, 0], yerr=D_eff[:,1], fmt="o", label=r"RW - $b=0$", color=c)
c=sns.color_palette()[1]
plt.plot(Pe, D_para, label=r"Aris - $b=0$", color=c)

plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"Peclet number Pe", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/comparison.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("../../oppgave/figures/comparison.pdf", "../../oppgave/figures/comparison.pdf"))
plt.show()
