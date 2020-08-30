import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd

plt.style.use("bmh")
sns.color_palette("hls", 1)
import matplotlib

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

Bruus = np.load("data/D_eff_Bruus.npy") 
Ivaar = np.load("data/D_eff_Ivar.npy")
Omega = np.logspace(-4, 4, 70)

Ivaar *= 105/2
Bruus *= 105/2

start_index = np.argmin(abs(Omega-10))
m, c, delta_m, delta_c = linear_regresion(np.log10(Omega[start_index:]), np.log10(Ivaar[start_index:]))
print(m, delta_m)

plt.plot(Omega, Bruus, label="Bruus")
plt.plot(Omega, Ivaar, "--", label="Our result")
plt.plot(Omega[start_index:], np.power(10, c+m*np.log10(Omega[start_index:])), "k", label="Slope -3.502(2)")

#plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"Driving frequency $\omega$", fontsize=14)
plt.ylabel(r"Geometrical factor $\tilde{g}$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("figures/D0.pdf")
plt.show()