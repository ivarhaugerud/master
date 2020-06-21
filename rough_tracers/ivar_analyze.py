import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import os

plt.style.use("bmh")
sns.color_palette("hls", 1)

def func(x, a, b, c):
    return a*(1+b/(x)**c)

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
"""
data = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.400000_Dm1.000000_U0.000000_dt0.000100_Nrw1000/tdata.dat")
t = data[:, 0]
mean = data[:, 1]
variance = data[:, 2]
t_cut = np.argmin(abs(t-500))

plt.figure(1)
plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
c=sns.color_palette()[0]
plt.plot(t, variance/(2*t), label=r"$r=0.4$")
plt.plot(t, 8.333670647214470595e-01*t/t)
plt.legend(loc="best", fontsize=12)
plt.show()
plt.figure(10)
plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.plot(t, variance/(2*t), color=c, label=r"$r=0.4$")
plt.plot(t, 8.333670647214470595e-01*t/t, "--", color=c)
"""
data  = np.loadtxt("data_square/RandomWalkers/Re0.000000_b1.500000_Dm0.062500_U1.000000_dt0.000032_Nrw1000/tdata.dat")

t = data[:, 0]
mean = data[:, 1]
variance = data[:, 2]
t_cut = np.argmin(abs(t-733.0))
Dm = 0.25

start_index = np.argmax(variance[1:]/t[1:])+1
D_eff = variance[start_index:]/(2*t[start_index:])

plt.figure(2)
analytic = 9.195160475665391
c=sns.color_palette()[1]
plt.plot(t[1:], variance[1:]/(Dm*2*t[1:]), color=c)
plt.plot(t[1:], analytic*t[1:]/t[1:])
print("Average effective diffusion coefficient: ", np.mean(variance[t_cut:]/(Dm*2*t[t_cut:])))
print("Standard deviation of effective diffusion coefficient: ", np.std(variance[t_cut:]/(Dm*2*t[t_cut:])))
print("Analytic value: ", analytic)
print("Relative error in percent: ", abs(100*(analytic-np.mean(variance[t_cut:]/(Dm*2*t[t_cut:])))/analytic))

plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.legend(loc="best", fontsize=12)
m,c ,delta_m, delta_c = linear_regresion(t, mean)
plt.show()
"""
data  = np.loadtxt("data_square/RandomWalkers/Re0.000000_b1.500000_Dm1.000000_U1.000000_dt0.000002_Nrw1000/tdata.dat")

t = data[:, 0]
mean = data[:, 1]
variance = data[:, 2]
t_cut = np.argmin(abs(t-204.0))
Dm = 1

start_index = np.argmax(variance[1:]/t[1:])+1
D_eff = variance[start_index:]/(2*t[start_index:])

plt.figure(2)
analytic = 1.246554654279746
c=sns.color_palette()[1]
plt.plot(t[1:], variance[1:]/(Dm*2*t[1:]), color=c)
plt.plot(t[1:], analytic*t[1:]/t[1:])
print("Average effective diffusion coefficient: ", np.mean(variance[t_cut:]/(Dm*2*t[t_cut:])))
print("Standard deviation of effective diffusion coefficient: ", np.std(variance[t_cut:]/(Dm*2*t[t_cut:])))
print("Analytic value: ", analytic)
print("Relative error in percent: ", abs(100*(analytic-np.mean(variance[t_cut:]/(Dm*2*t[t_cut:])))/analytic))

plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.legend(loc="best", fontsize=12)
m,c ,delta_m, delta_c = linear_regresion(t, mean)
plt.show()
"""
"""
data = np.loadtxt("data_square/Re0.0_b0.5.dat")
Pe = data[:, 0]
D_eff = data[:, 1]

plt.plot(Pe, D_eff, "o")
plt.yscale("log")
plt.xscale("log")
plt.show()
plt.plot(Pe, 105/2*(D_eff -1)/(Pe*Pe), "o", label=r"$b=0.5$")

data = np.loadtxt("data_square/Re0.0_b0.2.dat")
Pe = data[:, 0]
D_eff = data[:, 1]

plt.plot(Pe, 105/2*(D_eff -1)/(Pe*Pe), "o", label=r"$b=0.2$")

data = np.loadtxt("data_square/Re0.0_b0.0.dat")
Pe = data[:, 0]
D_eff = data[:, 1]

plt.plot(Pe, 105/2*(D_eff -1)/(Pe*Pe), "o", label=r"$b=0.0$")
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
#plt.yscale("log")
plt.show()
"""
"""
data  = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.000000_Dm0.001600_U1.000000_dt0.012500_Nrw3000/tdata.dat")

t = data[:, 0]
mean = data[:, 1]
variance = data[:, 2]
t_cut = np.argmin(abs(t-1500.0))
Dm = 0.0016 

start_index = np.argmax(variance[1:]/t[1:])+1
D_eff = variance[start_index:]/(2*t[start_index:])

plt.figure(2)
analytic = 1+2/(Dm*Dm*105)
c=sns.color_palette()[1]
plt.plot(t[1:], variance[1:]/(Dm*2*t[1:]), color=c)
plt.plot(t[1:], analytic*t[1:]/t[1:])
print("Average effective diffusion coefficient: ", np.mean(variance[t_cut:]/(Dm*2*t[t_cut:])))
print("Standard deviation of effective diffusion coefficient: ", np.std(variance[t_cut:]/(Dm*2*t[t_cut:])))
print("Analytic value: ", analytic)
print("Relative error in percent: ", abs(100*(analytic-np.mean(variance[t_cut:]/(Dm*2*t[t_cut:])))/analytic))

plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.legend(loc="best", fontsize=12)
m,c ,delta_m, delta_c = linear_regresion(t, mean)
plt.show()
"""
"""
plt.figure(1)
#plt.plot(t, func(t, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
c=sns.color_palette()[1]
plt.plot(t, mean, color=c)
print("Average effective diffusion coefficient: ", np.mean(variance[t_cut:]/(2*t[t_cut:])))
print("Standard deviation of effective diffusion coefficient: ", np.std(variance[t_cut:]/(2*t[t_cut:])))
print("Analytic value: ", analytic)
print("Relative error in percent: ", abs(100*(analytic-np.mean(variance[t_cut:]/(2*t[t_cut:])))/analytic))
print("Mean velocity: ", m, " pm ", delta_m)
print("Analytic mean: ", 1)
plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.show()
"""

data  = np.loadtxt("data_square/RandomWalkers/Re0.000000_b1.500000_Dm1.000000_U1.000000_dt0.000002_Nrw1000/Positions/xy_t15.000000.pos")
x = data[:, 1]
y = data[:, 2]

plt.scatter(x, y, s=2)
plt.show()

"""
data  = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.800000_Dm1.000000_U0.000000_dt0.000003_Nrw8000/tdata.dat")

t = data[:, 0]
mean = data[:, 1]
variance = data[:, 2]
t_cut = np.argmin(abs(t-2.0))


start_index = np.argmax(variance[1:]/t[1:])+1
print("start index: ", start_index)

D_eff = variance[start_index:]/(2*t[start_index:])
popt, pcov = curve_fit(func, t[start_index:], D_eff)

print(popt, pcov)
plt.plot(t[start_index:], D_eff)
plt.plot(t[start_index:], func(t[start_index:], *popt), 'r-', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
plt.show()
"""
"""
plt.show()
plt.figure(2)
analytic = 6.661579885472934670e-01
c=sns.color_palette()[1]
plt.plot(t[1:], variance[1:]/(2*t[1:]), color=c)
plt.plot(t[1:], analytic*t[1:]/t[1:])
print("Average effective diffusion coefficient: ", np.mean(variance[t_cut:]/(2*t[t_cut:])))
print("Standard deviation of effective diffusion coefficient: ", np.std(variance[t_cut:]/(2*t[t_cut:])))
print("Analytic value: ", analytic)
print("Relative error in percent: ", abs(100*(analytic-np.mean(variance[t_cut:]/(2*t[t_cut:])))/analytic))

plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.legend(loc="best", fontsize=12)
plt.show()
"""
"""
plt.figure(1)
c=sns.color_palette()[1]
plt.plot(t, mean, color=c, label=r"$r=0.2$")
plt.plot(t[1:], 0.2*t[1:]/t[1:])
plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.legend(loc="best", fontsize=12)

plt.show()
"""
"""
data = np.loadtxt("data_square/RandomWalkers/Re0.000000_b1.500000_Dm1.000000_U1.000000_dt0.000002_Nrw1000/declinedpos.dat")
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]
#declined = data[:, 3]

plt.scatter(x, y, s=1, color="k")
#plt.axis("equal")
plt.show()
"""
"""
"""
#plt.plot(t[1:], 100*declined[1:]/(t[1:]*5000/0.000010))
#lt.show()
"""
data = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.250000_Dm1.000000_U1.000000_dt0.000020_Nrw2000/declinedpos.dat")
x = data[:, 1]
y = data[:, 2]
plt.scatter(x, y)
plt.show()
"""
"""
data = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.400000_Dm1.000000_U0.000000_dt0.000010_Nrw5000/Positions/xy_t9.999990.pos")
x = data[:, 1]
y = data[:, 2]
plt.scatter(x, y)

data = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.400000_Dm1.000000_U0.000000_dt0.000010_Nrw5000/Positions/xy_t19.999980.pos")
x = data[:, 1]
y = data[:, 2]
plt.scatter(x, y)
#data = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.400000_Dm1.000000_U0.000000_dt0.000010_Nrw5000/Positions/xy_t1000.000000.pos")
#x = data[:, 1]
#y = data[:, 2]
#plt.scatter(x, y)
#plt.axis("equal")
plt.show()
"""
"""
plt.figure(10)
c=sns.color_palette()[1]
plt.plot(t, variance/(2*t), color=c, label=r"$r=0.2$")
plt.plot(t, 9.168006067034336626e-01*t/t, "--", color=c)

data = np.loadtxt("data_square/RandomWalkers/Re0.000000_b0.000000_Dm1.000000_U0.000000_dt0.000100_Nrw1000/tdata.dat")
t = data[:, 0]
mean = data[:, 1]
variance = data[:, 2]
t_cut = np.argmin(abs(t-500))

plt.figure(3)
c=sns.color_palette()[3]
plt.plot(t, variance/(2*t), label=r"$r=0.0$")
plt.plot(t, t/t, color=c)
plt.xlabel(r"Time $t$", fontsize=14)
plt.ylabel(r"Effective diffusion $D_{\parallel}/2Dt$", fontsize=14)
plt.legend(loc="best", fontsize=12)

plt.figure(10)
c=sns.color_palette()[3]
plt.plot(t, variance/(2*t), color=c, label=r"$r=0.0$")
plt.plot(t, t/t, "--", color=c)
plt.legend(loc="best", fontsize=12)

plt.show()
"""