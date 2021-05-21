import os
import numpy as np 
import matplotlib.pyplot as plt 

plt.style.use(['science','no-latex', "grid"])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
root = "../../master_latex/results/figures/rough/"

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


b = np.array([0.5, 0.7, 1.0, 1.3, 1.5])
Pe = np.array([1, 10, 100])
slopes = np.zeros((len(Pe), len(b), 2))

data = {}

for i in range(len(Pe)):
    for j in range(len(b)):
        data[str(i)+str(j)+"o"] = np.load("data_analyzed/prop_in_RZ_b"+str(b[j])+"_Pe" + str(Pe[i]) + "Re31.npy")
        data[str(i)+str(j)+"t"] = np.load("data_analyzed/occupation_time_main_channel_b"+str(b[j])+"_Pe" + str(Pe[i]) + "Re31.npy")


Nbins = 100

plt.figure(5)
for i in range(len(Pe)):
    for j in range(len(b)):
        data[str(i)+str(j)+"n"], data[str(i)+str(j)+"b"], patches = plt.hist(data[str(i)+str(j)+"t"], bins=Nbins, align="mid", density=True, alpha=0.5, histtype='bar', ec='black', label=r"$b=$"+str(b[j]))



for i in range(len(Pe)):
    plt.figure(1+i)
    for j in range(len(b)):
        x = data[str(i)+str(j)+"b"]
        y = data[str(i)+str(j)+"n"]

        x = x[np.where(abs(y) > 1e-12)]
        y = y[np.where(abs(y) > 1e-12)]
        slopes[i, j, 0], c, slopes[i, j, 1], delta_c = linear_regresion(x, np.log(y))
        plt.plot(data[str(i)+str(j)+"b"][1:], data[str(i)+str(j)+"n"], "o", markersize=3, color="C"+str(j), label=r"$b=$"+str(b[j]))
        plt.plot(data[str(i)+str(j)+"b"][1:], np.exp(c+slopes[i,j,0]*data[str(i)+str(j)+"b"][1:]), color="C"+str(j))

    plt.xlabel("Occupation time", fontsize=8)
    plt.ylabel('Normalized occurance density', fontsize=8)
    plt.legend(loc="best", fontsize=8)
    plt.yscale("log")
    name = (root+"MC_slope_Re31_Pe"+str(Pe[i])+".pdf")
    plt.savefig(name, bbox_inches="tight")
    os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()

prop_in_RZ = np.zeros((len(Pe), len(b), 2))

for i in range(len(Pe)):
    for j in range(len(b)):
        prop_in_RZ[i,j,0] = np.mean(data[str(i)+str(j)+"o"])
        prop_in_RZ[i,j,1] =  np.std(data[str(i)+str(j)+"o"])


for i in range(len(Pe)):
    plt.errorbar(b, 1-prop_in_RZ[i,:,0], fmt="o", markersize=3, label="Pe="+str(Pe[i]), color="C"+str(i))
    plt.errorbar(b, 1-prop_in_RZ[i,:,0], linewidth=0.5, color="C"+str(i))
plt.xlabel(r"Roughness $b$", fontsize=8)
plt.ylabel(r"Proportion of particles in MC", fontsize=8)
plt.legend(loc="best", fontsize=8)
name = (root+"occupation_prop_Re31_MC.pdf")
plt.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()

for i in range(len(Pe)):
    plt.errorbar(b, -1/(slopes[i,:,0]), yerr=slopes[i,:,1]/(slopes[i,:,0]**2), fmt="o", markersize=3, label="Pe="+str(Pe[i]), color="C"+str(i))
    plt.errorbar(b, -1/(slopes[i,:,0]), linewidth=0.5, color="C"+str(i))

plt.xlabel(r"Roughness $b$", fontsize=8)
plt.ylabel(r"Characteristic CC occupation time $\tau$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, framealpha=0.5)
plt.yscale("log")
name = root + "slopes_vs_b_Re31_MC.pdf"
plt.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()