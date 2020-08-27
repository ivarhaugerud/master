import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd

plt.style.use("bmh")
sns.color_palette("hls", 1)
import matplotlib

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

Bruus = np.load("D_eff_Bruus.npy") 
Ivaar = np.load("D_eff_Ivar.npy")
Omega = np.logspace(-3, 3, 40)

plt.plot(Omega, Bruus*105/2, label="Bruus")
plt.plot(Omega, Ivaar*105/2, "--", label="Our result")
plt.xscale("log")
plt.xlabel(r"Driving frequency $\omega$", fontsize=14)
plt.ylabel(r"Geometrical factor $\tilde{g}$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("figures/D0.pdf")
plt.show()