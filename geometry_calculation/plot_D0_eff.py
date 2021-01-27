import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd
import matplotlib
import os

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

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

D_eff = np.load("data_test/D_eff_2.npy")[:,:,0]*105/2

omega   = np.logspace(-3, 3, 30)
schmidt = np.logspace(-3, 3, 60)
print(schmidt)


fig = plt.figure(1)
x_, y_ = np.meshgrid(omega, schmidt)

ax1 = plt.contourf(x_,y_, np.transpose((D_eff)) )#, levels=np.linspace(0, 1, 101))
cbar = fig.colorbar(ax1)
#cbar.ax.locator_params(nbins=10)
cbar.ax.set_ylabel(r'Geometric factor $\tilde{g}$', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r"Diffusion coefficient $D$", fontsize=14)
plt.xlabel(r"Kinematic viscosity $\nu$", fontsize=14)
plt.savefig("figures/D_0_eff.pdf", bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%("figures/D_0_eff.pdf", "figures/D_0_eff.pdf"))
#plt.savefig("../figures/g_tilde_D0.pdf")
plt.show()

plt.style.use("bmh")
sns.color_palette("hls", 1)


for i in range(len(schmidt[:])):
    if int(i)%10 == 0:
        plt.plot(omega, D_eff[:, i], label="$\log_{10}$(Sc) = " + str(np.log10(schmidt[i])))

#plt.ylabel(r"Geometric factor $\tilde{g}$", fontsize=14)
plt.xlabel(r"Driving frequency $\omega$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.xscale("log")
plt.yscale("log")
plt.show()

for i in range(len(omega[:])):
    if int(i)%10 == 0:
        plt.plot(schmidt, D_eff[i, :], label="$\log_{10}(\omega) = $" + str(np.log10(omega[i])))

#plt.ylabel(r"Geometric factor $\tilde{g}$", fontsize=14)
plt.xlabel(r"Schmidt number Sc", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.xscale("log")
plt.yscale("log")
plt.show()




