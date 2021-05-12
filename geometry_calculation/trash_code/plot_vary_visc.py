import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import os 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

visc = np.logspace(-2, 2, 10)
data = np.load("data_test/tdata_tau3_eps05_Lx4.448_D0.1_dt0.004_vary_viscF0.npy")
D = np.zeros((len(visc)))
difference = np.zeros(np.shape(D))
T = int(3.0/0.004)
tau = 3.0

for i in range(len(visc)):
	D[i] = sci.trapz(  np.trim_zeros(data[i, :, 8])[-T:],  np.trim_zeros(data[i, :, 0])[-T:] )/tau
	difference[i] = abs(D[i] - sci.trapz(  np.trim_zeros(data[i, :, 8])[-2*T:-T],  np.trim_zeros(data[i, :, 0])[-2*T:-T] )/tau)/D[i]
	#plt.plot(np.trim_zeros(data[i,   :, 0])/tau, np.trim_zeros(data[i,   :, 8]))
	#plt.plot(np.trim_zeros(data[i, :, 0])[-T:]/tau, np.trim_zeros(data[i, :, 8])[-T:])
	#plt.show()

plt.plot(visc, difference, "o")
plt.yscale("log")
plt.legend(loc="best")
plt.xscale("log")
plt.show()

plt.plot(visc, D, "o", markersize=3)
plt.xscale("log")
plt.xlabel(r"Viscosity $[\nu]$", fontsize=8)
plt.ylabel(r"Effective diffusion coefficient $[D_\parallel]$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = "figures/vary_visc.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
