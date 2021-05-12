import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import os
import scipy.interpolate as sci 
root = "../master_latex/results/figures/rough/"

plt.style.reload_library()
plt.style.use(['science','no-latex', 'grid'])
#sns.color_palette("hls", 1)
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

Pe = np.linspace(1, 1e3, int(1e5))
analytic = interpool(Pe)

plt.figure(1)
c= sns.color_palette()[0]
plt.scatter(Pe2, D_eff[:, 0], s=8.5, label=r"RW: $b=1.5$", color="C1")
plt.plot(Pe, analytic, label=r"Brenner $b=1.5$", color="C1")

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
counter = 0

Brenner = np.loadtxt("data_square/Re0.0_b0.0_res1000.dat")
Pe_temp = Brenner[:, 0]
D_eff_B = Brenner[:, 1]
plt.figure(2)

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    t = data[:, 0]
    variance = data[:, 2]
    t_cut = np.argmin(abs(t-t_cuts[i]))


    Ds = variance[1:]/(Dm[i]*2*t[1:])
    t_prime = t[1:]
    indx = int(len(t_prime)*2/3)

    a = np.argmin(abs(Pe_temp-1/Dm[i]))
    
    D_eff[i, 0] = np.mean(variance[indx:]/(Dm[i]*2*t[indx:]))
    D_eff[i, 1] = np.std((variance[indx:]/(2*t[indx:])))

    plt.plot(np.linspace(0, 1, len(t_prime)), Ds*Dm[i], label=r"Pe$=%3.0f$" % (1/Dm[i]), color="C"+str(counter))
    plt.plot(np.linspace(0, 1, len(t_prime)), D_eff_B[a]*np.ones(len(t_prime))*Dm[i], color="C"+str(counter))
    plt.fill_between(np.linspace(0, 1, len(t_prime)), Ds*Dm[i],  D_eff_B[a]*np.ones(len(t_prime))*Dm[i], color="C"+str(counter), alpha=0.5)
    counter += 1

plt.yscale("log")
plt.ylabel(r"Effective diffusion coefficient $D_{\parallel}$ [$D_m$/Pe]", fontsize=8)
plt.xlabel(r"Time $[T_{max}]$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=3)
name = root+"vary_Pe_RW.pdf"
plt.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))

plt.figure(1)
Pe = np.linspace(1, 1e3, int(1e5))
D_para = 1+2*Pe*Pe/105

Brenner = np.loadtxt("data_square/Re0.0_b0.0_res1000.dat")
Pe_temp = Brenner[:, 0]
D_eff_B = Brenner[:, 1]
interpool = sci.interp1d(Pe_temp, D_eff_B)

c=sns.color_palette()[2]
plt.scatter(Pe2, D_eff[:, 0], s=8.5, label=r"RW: $b=0$", color="C2")
plt.plot(Pe, interpool(Pe), label=r"Brenner $b=0$", color="C2")

plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"Peclet number Pe", fontsize=8)
plt.ylabel(r"Effective diffusion $D_{\parallel}$ $[D_m]$", fontsize=8)
plt.plot(Pe, D_para, "--", label=r"Taylor-Aris", color="k") 
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
name = root+"comparison.pdf"
plt.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()
