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


data_strings = ["data_square/Re0.0_b0.003_res2200.dat",
                "data_square/Re0.0_b0.006_res2200.dat",
                "data_square/Re0.0_b0.012_res2000.dat",
                "data_square/Re0.0_b0.025_res1100.dat",
                "data_square/Re0.0_b0.05_res1000.dat",
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

b = np.array([0.003, 0.006, 0.012, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9])
max_index = np.zeros(len(data_strings), dtype='int')
max_Pe = np.zeros(len(data_strings))

for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    Pe = data[:, 0]
    g_tilde = 105*(data[:, 1]-1)/(2*Pe*Pe)

    f = sci.interp1d(Pe, g_tilde)
    fine_Pe = np.linspace(min(Pe), max(Pe), int(1e4))
    fine_g_tilde = f(fine_Pe)

    max_index[i] = np.argmax(fine_g_tilde)
    max_Pe[i]    = fine_Pe[max_index[i]]

    if i < 6:
        plt.plot(fine_Pe, f(fine_Pe), label=r"$b=$" + str(b[i]))
        plt.plot(fine_Pe[max_index[i]], fine_g_tilde[max_index[i]], "ko")

plt.legend(loc="best")
plt.xscale("log")
#plt.yscale("log")
plt.show()

plt.figure(3)
m, c, delta_m, delta_c = linear_regresion( np.log(b[:-5]), np.log(max_Pe[:-5]) )
plt.plot(b[:-5], np.exp(c+m*np.log(b[:-5])), label=r"Slope: "+str(m)[0:5] + "(" + str(delta_m*100)[:1] + ")")
print(m, delta_m)
m, c, delta_m, delta_c = linear_regresion( np.log(b[-5:]), np.log(max_Pe[-5:]) )
print(m, delta_m)

plt.plot(b[-5:], np.exp(c+m*np.log(b[-5:])), label=r"Slope: "+str(m)[0:5] + "(" + str(delta_m*100)[:1] + ")")
plt.plot(b, max_Pe, "o")
plt.yscale("log")
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
plt.ylabel(r"Peclet number Pe", fontsize=14)
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.show()
