import os
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import scipy.interpolate as scp 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	
root = "../../../master_latex/results/"
base = "../data_test/vary_geometry/"
Map = matplotlib.cm.get_cmap('Spectral_r')

epsilon = np.arange(0.1, 0.601, 0.10)
Lx = np.array([12.0, 9.0, 6.0, 3.0])
kappa = 2*np.pi/Lx 
tau = 3.0

D = np.zeros((len(kappa), len(epsilon)))
difference = np.zeros(np.shape(D))
T = 750 

plt.figure(1)
kappa_cont = np.linspace(min(kappa), max(kappa), int(1e4))

for i in range(len(kappa)):
	for j in range(len(epsilon)):
		data = np.loadtxt(base+"Lx" +  str(Lx[i]) +"_tau3.0_eps"+str(epsilon[j])[:3]+"_nu1.2_D1.0_fzero0.0_fone12.0_res150_dt0.004/tdata.dat")
		D[i, j] = sci.trapz(  data[:, 8][-T:],  data[:, 0][-T:] )/tau
		difference[i, j] = abs(D[i, j] - sci.trapz(  np.trim_zeros(data[:, 8])[-2*T:-T],  np.trim_zeros(data[:, 0])[-2*T:-T] )/tau)/D[i,j]

		plt.plot(np.trim_zeros(data[:, 0])[-T:]/tau, np.trim_zeros(data[:, 8])[-T:])
		plt.xlabel(r" Time [periods]", fontsize=8)
		plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.show()

plt.figure(2)
for i in range(len(epsilon)):
	plt.plot(kappa, difference[:, i], "o", markersize=3, label=r"$\epsilon=%3.2f$" % epsilon[i])
plt.yscale("log")
plt.legend(loc="best")
plt.show()


plt.figure(3)
for i in range(len(epsilon)):
	plt.plot(kappa, D[:, i], "o", color="C"+str(i), markersize=3, label=r"$\epsilon=%3.2f$" % epsilon[i])
	plt.plot(kappa, D[:, i], color="C"+str(i))

plt.xlabel(r" Wave number $\kappa$", fontsize=8)
plt.ylabel(r" Effective Diffusion Coefficient $ D_\parallel $",  fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/vary_geo.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()


fig = plt.figure(4)
x_, y_ = np.meshgrid(kappa, epsilon)
ax1 = plt.contourf(x_,y_, np.transpose(D), cmap=Map, levels=18)
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Effective diffusion coefficient $D_\parallel$', fontsize=8)
plt.ylabel(r"Boundary amplitude $\epsilon$",    fontsize=8)
plt.xlabel(r"Wave number $\kappa$",       fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
filename = root + "figures/vary_geo_contour.pdf"
plt.savefig(filename, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()