import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import os
import scipy.interpolate as sci 

plt.style.use("bmh")
sns.color_palette("hls", 1)

def func2(x, a, b):
    return a*np.exp(-b*x)

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

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


data_strings = ["data_square/Re0.0_b0.0_res700.dat",
                "data_square/Re0.0_b0.1_res700.dat",
                "data_square/Re0.0_b0.2_res400.dat",
                "data_square/Re0.0_b0.3_res400.dat",
                "data_square/Re0.0_b0.4_res320.dat",
                "data_square/Re0.0_b0.5_res300.dat",
                "data_square/Re0.0_b0.6_res300.dat",
                "data_square/Re0.0_b0.7_res270.dat",
                "data_square/Re0.0_b0.8_res250.dat",
                "data_square/Re0.0_b0.9_res250.dat",
                "data_square/Re0.0_b1.0_res235.dat",
                "data_square/Re0.0_b1.1_res200.dat",
                "data_square/Re0.0_b1.2_res220.dat",
                "data_square/Re0.0_b1.3_res210.dat",
                "data_square/Re0.0_b1.4_res210.dat",
                "data_square/Re0.0_b1.5_res200.dat",
                "data_square/Re0.0_b1.6_res195.dat",
                "data_square/Re0.0_b1.7_res190.dat",
                "data_square/Re0.0_b1.8_res185.dat",
                "data_square/Re0.0_b1.9_res180.dat"]

Pe = np.logspace(0, 3, 40)
D_eff = np.zeros((len(data_strings), len(Pe), 3))
b = np.linspace(0, 1.9, len(data_strings))
max_index = np.zeros(len(data_strings)-1, dtype='int')

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    Pe = data[:, 0]
    D_eff[i, :, 0] = data[:, 1]
    D_eff[i, :, 1] = data[:, 2]
    D_eff[i, :, 2] = 105*(data[:, 1]-1)/(2*Pe*Pe)
    if i > 0:
    	max_index[i-1] = np.argmax(D_eff[i, :, 2])

    if i % 3 == 0:
	    plt.figure(1)
	    plt.plot(Pe, D_eff[i, :, 0], label=r"$b=$"+str(b[i]))

	    plt.figure(2)
	    plt.plot(Pe, D_eff[i, :, 2], label=r"$b=$"+str(b[i]))

plt.figure(5)
for i in range(len(Pe)):
	if i%7 == 0:
		plt.plot(b, D_eff[:, i, 2], label=r"Pe=" + str(Pe[i])[:5])


plt.figure(1)
plt.xlabel(r"Peclet number Pe")
plt.ylabel(r"Effective diffusion coefficient")
plt.legend(loc="best", fontsize=12)
plt.xscale("log")
plt.yscale("log")
plt.savefig("../../oppgave/figures/R0.pdf")

plt.figure(2)
plt.xlabel(r"Peclet number Pe")
plt.ylabel(r"Geometric factor $\tilde{g}$")
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/g_tilde.pdf")

plt.figure(3)
m, c, delta_m, delta_c = linear_regresion( np.log(b[1:-5]), np.log(Pe[max_index[:-5]]) )
plt.plot(b[1:-5], np.exp(c+m*np.log(b[1:-5])), label=r"Slope: "+str(m)[0:5] + "(" + str(delta_m*100)[:1] + ")")
print(m, delta_m)
m, c, delta_m, delta_c = linear_regresion( np.log(b[-5:]), np.log(Pe[max_index[-5:]]) )
print(m, delta_m)

plt.plot(b[-5:], np.exp(c+m*np.log(b[-5:])), label=r"Slope: "+str(m)[0:5] + "(" + str(delta_m*100)[:1] + ")")
plt.plot(b[1:], Pe[max_index], "o")
plt.yscale("log")
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
plt.ylabel(r"Peclet number Pe", fontsize=14)
plt.xlabel(r"Roughness $b$", fontsize=14)

fig = plt.figure(4)
x_, y_ = np.meshgrid(Pe, b)
plt.ylabel(r"Roughness $b$", fontsize=14)
plt.xlabel(r"Peclet number Pe", fontsize=14)
ax1 = plt.contourf(x_,y_, D_eff[:, :, 2], levels=np.linspace(-13, 45, 15))
cbar = fig.colorbar(ax1)
cbar.ax.set_ylabel(r'Geometric factor $\tilde{g}$', fontsize=14)
plt.xscale('log')
#plt.savefig("../figures/relative_diff.pdf", bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%("../figures/relative_diff.pdf", "../figures/relative_diff.pdf"))
plt.xlabel(r"Peclet number Pe", fontsize=14)
plt.ylabel(r"Roughness $b$", fontsize=14)

plt.figure(5)
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.ylabel(r"Geometric factor $\tilde{g}$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.savefig("../../oppgave/figures/R0_g_vs_b.pdf")
plt.show()
