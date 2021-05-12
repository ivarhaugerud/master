import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import seaborn as sns
import os

plt.style.use(['science','no-latex', "grid"])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
root = "../master_latex/results/"

Reynolds = [0, 1, 3.16, 10, 31.6, 100]
b_s = np.linspace(0.0, 0.7, 8)
RZ = np.zeros((len(Reynolds), len(b_s)))

prop_RZs = np.load("analyzed_data/RZ_prop_try2.npy")
print(np.shape(prop_RZs))
bs = np.linspace(0, 1.9, 20)
colors_b = plt.cm.plasma(np.linspace(0,1,int(len(Reynolds))))
index_RZs = [0, 2, 3, 4, 5, 6]
exps = [0, 0, 0.5, 1.0, 1.5, 2]

for i in range(len(Reynolds)):
	N_b = len(np.trim_zeros(prop_RZs[index_RZs[i],:], "b"))
	b = np.linspace(0, (N_b-1)/10, N_b)
	
	if i != 0:
		plt.scatter(b, np.trim_zeros(prop_RZs[index_RZs[i],:], "b"), s=3, color=sns.color_palette()[i], label="Re = $10^{%2.1f}$" % exps[i] )
	else:
		plt.scatter(b, np.trim_zeros(prop_RZs[index_RZs[i],:], "b"), s=3, color=sns.color_palette()[i], label="Re = 0")
	plt.plot(b, np.trim_zeros(prop_RZs[index_RZs[i],:], "b"), "-", linewidth=1, color=sns.color_palette()[i])

plt.plot(bs, bs/2, "k", label="Cavity area")
plt.ylabel(r"Relative reciruclation zone area $\frac{A_{RZ}}{A}$", fontsize=8)
plt.xlabel(r"Roughness $b$", fontsize=8)
plt.legend(loc="best", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
name = root+"figures/reciruclation_area.pdf"
plt.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()
"""
#plt.xscale("log")
#plt.yscale("log")

prop_RZs[np.where(prop_RZs == 0)] = -0.02

fig = plt.figure(i)
x_, y_ = np.meshgrid(bs, Reynolds[1:])
ax1 = plt.contourf(x_,y_, prop_RZs[1:, :], np.linspace(0, np.max(np.max(prop_RZs)), 30))
cbar = fig.colorbar(ax1)
cbar.ax.set_ylabel(r'Effective diffusion coeff $D$', fontsize=14)
plt.yscale('log')
#plt.savefig("../figures/relative_diff.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/relative_diff.pdf", "../figures/relative_diff.pdf"))
plt.ylabel(r"Reynolds number Re", fontsize=14)
plt.xlabel(r"Roughness $b$", fontsize=14)
#plt.savefig("figures/g_tilde_contour.pdf")
plt.show()
"""