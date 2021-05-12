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











import numpy as np 
import matplotlib.pyplot as plt 
import h5py

b = 1.5
Pe = 100
dirr  = "data_square/RandomWalkers/Re31.000000_b%7.6f" % b + "_Dm0.010000_U1.000000_dt0.010000_Nrw1000/Positions/"
dirr2 = "data_square/RandomWalkers/Re0.000000_b%7.6f" % b + "_Dm0.010000_U1.000000_dt0.010000_Nrw1000/Positions/"
dt = 0.1
T = 100
"""
Pe = 1
dirr  = "data_square/RandomWalkers/Re31.000000_b%7.6f" % b + "_Dm1.000000_U1.000000_dt0.001000_Nrw1000/Positions/"
dirr2 = "data_square/RandomWalkers/Re0.000000_b%7.6f" % b + "_Dm1.000000_U1.000000_dt0.001000_Nrw1000/Positions/"
dt = 0.1
T = 50
"""
N = int(T/dt)-1
Nrw = 1000 
ts = np.arange(0, T-1e-5, dt)



file = "data_square/U_Re0.0_b"+str(b)+".h5"
f = h5py.File(file, 'r')
vector = np.array(list(f["VisualisationVector"]["0"]))
geometry = np.array(list(f["Mesh"]["mesh"]["geometry"]))
non_RZ = []
RZ = []



x = np.zeros((N, Nrw))
y = np.zeros((N, Nrw))
x2 = np.zeros((N, Nrw))
y2 = np.zeros((N, Nrw))

for i in range(len(ts)-1):
    try:
        pos  = np.loadtxt(dirr+"xy_t%7.6f" % ts[i+1] +".pos")
        pos2 = np.loadtxt(dirr2+"xy_t%7.6f" % ts[i+1] +".pos")
    except:
        try:
            pos  = np.loadtxt(dirr+"xy_t%8.7f" % ts[i+1] +".pos")
            pos2 = np.loadtxt(dirr2+"xy_t%8.7f" % ts[i+1] +".pos")
        except:
            pos  = np.loadtxt(dirr+"xy_t%6.5f" % ts[i+1] +".pos")
            pos2 = np.loadtxt(dirr2+"xy_t%6.5f" % ts[i+1] +".pos")

    x[i, :] = pos[:, 1]
    y[i, :] = pos[:, 2]

    x2[i, :] = pos2[:, 1]
    y2[i, :] = pos2[:, 2]

for i in range(int(N/4)-1):
    
    plt.figure(1)
    plt.clf()
    e = np.linspace(min(x[i*4,:]), max(x[i*4, :]), int(1e5))
    plt.scatter(x[i*4, :], y[i*4, :],   s=2, color="b")
    plt.scatter(x2[i*4, :], y2[i*4, :], s=2, color="r")
    xi = np.ones(len(e))*(1-b/2)
    xi[np.where( (e - b/2 + 2*b)%(int(2*b)) > b )] = 1+b/2
    plt.plot(e, xi, "k")
    plt.plot(e, -xi, "k")
    plt.pause(0.01)
    
    plt.figure(2)
    plt.clf()
    plt.hist(x[i*4,:],  bins=80, align="mid", color="b", density=True, alpha=0.5, histtype='bar', ec='black', label=r"Re=31")
    plt.hist(x2[i*4,:], bins=80, align="mid", color="r", density=True, alpha=0.5, histtype='bar', ec='black', label=r"Re=0")
    plt.legend(loc="best")
    plt.pause(0.01)

    plt.figure(3)
    denom = (2*ts[i*4]/Pe)
    plt.plot(ts[i*4], np.var(x[i*4, :])/denom, "bo")
    plt.plot(ts[i*4], np.var(x2[i*4, :])/denom, "ro")
    plt.pause(0.01)
plt.show()
