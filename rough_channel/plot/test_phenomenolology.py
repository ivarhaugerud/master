import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import seaborn as sns
import os

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

Reynolds = [0, 1, 3.16, 10, 31.6, 100]
b_s = np.linspace(0.0, 0.7, 8)
RZ = np.zeros((len(Reynolds), len(b_s)))

alpha = np.load("analyzed_data/RZ_prop.npy")
b = np.linspace(1e-2, 1.9, 20)
alpha *= 2/b 
alpha[:, 0] = 1
for i in range(len(Reynolds)):
	plt.plot(b, alpha[i, :])
plt.show()
for i in range(len(Reynolds)):
	plt.plot(b, b +b*b*(2-np.log10(Reynolds[i]))/2)
plt.show()

for i in range(len(Reynolds)):
	plt.plot(b, 1e-3*(3*alpha[i, :]*b-4)/(alpha[i,:]*b*(4*b-2*alpha[i, :]*b*b)))
#plt.yscale("log")
plt.show()