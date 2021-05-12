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

Brenner = np.loadtxt("data_square/Pe0_300.dat")
b_original  = Brenner[:, 0]
D_eff       = Brenner[:, 1]

b_original2 = np.zeros(len(b_original)+1)
b_original2[:-1] = b_original
b_original2[-1]  = 2.0

D_eff_n0 = np.zeros(len(b_original)+1)
D_eff_n0[:-1] = D_eff
D_eff_n0[-1]  = 0

interpool = sci.interp1d(b_original2, D_eff_n0)
b_test = np.linspace(0, 2-1e-5, int(1e4))

n0 = interpool(b_test)


b_N = 5
data_strings = ["data_square/Pe0_frac_res800_bmax0.25.dat",
                "data_square/Pe0_frac_res600_bmax0.5625.dat",
                "data_square/Pe0_frac_res300_bmax1.125.dat",
                "data_square/Pe0_frac_res300_bmax1.6875.dat",
                "data_square/Pe0_frac_res300_bmax1.985.dat"]

total_Bs = int(5*len(data_strings)+1)
bs = np.zeros(total_Bs)
diff = np.zeros(total_Bs)

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    b = data[:, 0]
    diffusion = data[:, 1]

    bs[i*5:(i+1)*5] = b
    diff[i*5:(i+1)*5] = diffusion

bs[-1]   = 2.0
diff[-1] = 0

interpool = sci.interp1d(bs, diff)
n1 = interpool(b_test)


b_N = 5
data_strings = ["data_square/Pe0_frac_2_res800_bmax0.25.dat",
                "data_square/Pe0_frac_2_res400_bmax0.5625.dat",
                "data_square/Pe0_frac_2_res400_bmax0.875.dat",
                "data_square/Pe0_frac_2_res400_bmax0.9375.dat",
                "data_square/Pe0_frac_2_res350_bmax1.25.dat",
                "data_square/Pe0_frac_2_res300_bmax1.5625.dat",
                "data_square/Pe0_frac_2_res325_bmax1.875.dat",
                "data_square/Pe0_frac_2_res250_bmax1.98.dat"]

total_Bs = int(5*len(data_strings)+1)
bs = np.zeros(total_Bs)
diff = np.zeros(total_Bs)

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    b = data[:, 0]
    diffusion = data[:, 1]

    bs[i*5:(i+1)*5] = b
    diff[i*5:(i+1)*5] = diffusion

bs[-1]   = 2.0
diff[-1] = 0

interpool = sci.interp1d(bs, diff)
n2 = interpool(b_test)

plt.plot(b_test, n0, label=r"$n=1$")
plt.plot(b_test, n1, label=r"$n=2$")
plt.plot(b_test, n2, label=r"$n=3$")

"""
prop_open = 2*b_test/(2+b_test)
plt.plot(b_test, (1-prop_open)*(1+prop_open), label="Reguera")
"""
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("figures/Pe0_compare.pdf")
plt.show()

plt.plot(b_test, n1-n0, label=r"$n_1-n_0$")
plt.plot(b_test, n2-n0, label=r"$n_2-n_0$")
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.ylabel(r"Rel. diff. effective diffusion $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.legend(loc="best", fontsize=12)
#plt.yscale("log")
plt.savefig("figures/Pe0_frac_difference_n0.pdf")
plt.show()

plt.plot(b_test, (n1-n0)/n0, label=r"$(n_1-n_0)/n_0$")
plt.plot(b_test, (n2-n1)/n1, label=r"$(n_2-n_1)/n_1$")
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.ylabel(r"Rel. diff. effective diffusion $D_{\parallel}$ $[D_m]$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.yscale("log")
plt.savefig("figures/Pe0_frac_rel_difference.pdf")
plt.show()