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
root = "../../../master_latex/results/figures/"

fig,ax = plt.subplots(figsize=(6,6))
sizes = np.array([5, 10, 20, 40, 80, 160, 320, 640, 1280])
f = 1e-8
tau = 2 
visc = (2-1/2)/3
U = sizes*sizes*f/(3*visc*4)
error = np.zeros(len(sizes))

for i in range(len(sizes)):
	ux = np.loadtxt("../data/benchmark_testing_pousielle"+str(sizes[i])+"_ux.txt")
	y = np.linspace(-1, 1, sizes[i])
	U_ana = 3*U[i]*(1-y*y)/2
	U_num = np.mean(ux, axis=0)

	plt.plot(y, U_num/U[i], label=r"$N_y$=" + str(sizes[i]))
	error[i] = (sci.trapz(abs(U_ana-U_num), y)/sizes[i])/U[i]

plt.plot(y, U_ana/U[i], color="k", label="analytic")
plt.xlabel(r" Vertical position [$a$]", fontsize=8)
plt.ylabel(r" Horizontal velocity [$U$]",  fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8)
filename = root + "rough/change_U.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

plt.plot(sizes, error, "o", markersize=3, label="Measured value")

plt.plot(sizes[1:], 8*np.power(sizes[1:], -2.0), label=r"$N_y^{-2}$")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Number of vertical lattice sites $N_y$", fontsize=8)
plt.ylabel(r"Relative difference in U",  fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8)
filename = root + "rough/relative_change_U.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()