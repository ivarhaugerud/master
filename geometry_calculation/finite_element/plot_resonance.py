import os
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
import scipy.interpolate as scp
import scipy.signal as scs

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'	
root = "../../../master_latex/results/"



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


D_parallels = np.load("data/D_para_vary_kappa_F0s_D_nu1.25_tau3.0.npy")


F0s      = np.array([12, 24, 36, 48, 60, 72])/1.25
Ds       = np.array([0.6, 0.8, 1.0, 1.2, 1.4]) 
kappas   = np.arange(0.1, 1.501, 0.1) 

kapppa_cont = np.linspace(min(kappas), max(kappas), int(1e4))
kappa_res = np.zeros((len(F0s), len(Ds)))
slopes = np.zeros((len(F0s), 2))

for i in range(len(F0s)):
	for j in range(len(Ds)):
		interpoo = scp.interp1d(kappas, D_parallels[:, i, j], "cubic")(kapppa_cont)
		#plt.plot(kappas, D_parallels[:, i, j])
		if len(scs.argrelmax(D_parallels[:,i,j])[0]) == 1:
			kappa_res[i, j] = kappas[scs.argrelmax(D_parallels[:,i,j])]
			#plt.plot(kappa_res[i, j], D_parallels[scs.argrelmax(D_parallels[:,i,j]),i,j], "ko")
			kappa_res[i, j] = kapppa_cont[np.argmax([interpoo])]
		if len(scs.argrelmax(D_parallels[:,i,j])[0]) > 1:
			indxs = scs.argrelmax(D_parallels[:,i,j])[0]
			max_index = np.argmax(D_parallels[indxs,i,j])
			kappa_res[i, j] = kappas[indxs[max_index]]
			kappa_res[i, j] = kapppa_cont[np.argmax(interpoo)]
			#plt.plot(kappa_res[i, j], D_parallels[indxs[max_index],i,j], "ko")
			print("two maxima for F0=", F0s[i], " and D=", Ds[j], "for kappa=", kappas[indxs])


for i in range(len(Ds)):
	plt.plot(F0s, (2*np.pi/kappa_res[:, i]), label=r"$D=%3.2f$" % Ds[i])
#plt.xscale("log")
plt.xlabel("F0")
plt.ylabel("lambda")
plt.legend(loc="best")
plt.show()
M = np.zeros((len(Ds), 2))

for i in range(len(Ds)):
	plt.plot(F0s[1:]*Ds[i]**(0.3333), kappa_res[1:, i], "o", label=r"$D=%3.2f$" % Ds[i])
	m, c, delta_m, delta_c = linear_regresion(np.log(F0s[1:]*Ds[i]**(0.3333)), np.log(kappa_res[1:, i]))
	plt.plot(F0s[1:]*Ds[i]**(0.3333), np.exp(c+m*np.log(F0s[1:]*Ds[i]**(0.3333))))
	print(m, c, delta_m, delta_c)
	M[i,0] = m
	M[i, 1] = delta_m
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$ F_0 D^{1/3}$")
plt.ylabel(r"$\kappa$")
plt.legend(loc="best")
plt.show()

for i in range(len(Ds)):
	plt.plot( (F0s[1:]*Ds[i]**(0.3333))**(-0.6666), kappa_res[1:, i], "o", label=r"$D=%3.2f$" % Ds[i])
plt.legend(loc="best")
plt.xlabel(r"$F_0^{-2/3} D^{-2/9}$")
plt.ylabel(r"\kappa")
plt.show()

for i in range(len(Ds)):
	x = (F0s[1:]*Ds[i]**(0.3333))**(0.6666)
	plt.plot( x, x, "-", label=r"$D=%3.2f$" % Ds[i])
	plt.plot( x, 2*np.pi/kappa_res[1:, i], "o")
plt.legend(loc="best")
plt.xlabel(r"$F_0^{-2/3} D^{-2/9}$")
plt.show()

plt.errorbar(Ds, M[:,0], yerr=M[:,1], fmt="o")
plt.show()


fig = plt.figure(1)
x_, y_ = np.meshgrid(F0s, Ds)
Map = matplotlib.cm.get_cmap('Spectral_r')
kappa_res[np.where(kappa_res == 0)] = -1
ax1 = plt.contourf(x_,y_, np.transpose((2*np.pi/(kappa_res))), levels=np.linspace(np.amin(2*np.pi/kappa_res), np.amax(2*np.pi/kappa_res), 15), cmap=Map)
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Resonance wave length $\lambda_{res}$', fontsize=8)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"External force  $F_0$", fontsize=8)
plt.ylabel(r"Diffusion coefficient $D$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/D_0_eff.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
















D_parallels = np.load("data/D_para_vary_kappa_F0s_D_nu1.25_tau3.0.npy")


F0s      = np.array([12, 24, 36, 48, 60, 72])
Ds       = np.array([0.6, 0.8, 1.0, 1.2, 1.4]) #np.logspace(-3, 0, 20) #np.arange(0.6, 2.7, 0.1)
kappas   = np.arange(0.1, 1.501, 0.1) #np.array([0.85, 0.90, 0.95, 1.0, 1.05, 1.10, 1.15])#np.arange(0.2, 1.75, 0.1)
"""
"""
nu = 2.25
F0s      = np.logspace(0.5,  3, 25)/2.25
Ds       = np.logspace(-2, 1, 25) #np.logspace(-3, 0, 20) #np.arange(0.6, 2.7, 0.1)
kappas   = np.array([0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6, 4.4, 5.0, 6, 7.0, 8, 9]) #np.array([0.85, 0.90, 0.95, 1.0, 1.05, 1.10, 1.15])#np.arange(0.2, 1.75, 0.1)
D_parallels = np.load("data/D_para_vary_kappa_F0s_D_nu2.25_tau3.0.npy")

kapppa_cont = np.linspace(min(kappas), max(kappas), int(1e4))
nu = 1.25


kappa_res = np.zeros((len(F0s), len(Ds)))
slopes = np.zeros((len(F0s), 2))

for i in range(len(F0s)):
	for j in range(len(Ds)):
		interpoo = scp.interp1d(kappas, D_parallels[:, i, j], "cubic")(kapppa_cont)
		plt.plot(kappas, D_parallels[:, i, j])
		if len(scs.argrelmax(D_parallels[:,i,j])[0]) == 1:
			kappa_res[i, j] = kappas[scs.argrelmax(D_parallels[:,i,j])]
			plt.plot(kappa_res[i, j], D_parallels[scs.argrelmax(D_parallels[:,i,j]),i,j], "ko")

		if len(scs.argrelmax(D_parallels[:,i,j])[0]) > 1:
			indxs = scs.argrelmax(D_parallels[:,i,j])[0]
			max_index = np.argmax(D_parallels[indxs,i,j])
			kappa_res[i, j] = kappas[indxs[max_index]]
			plt.plot(kappa_res[i, j], D_parallels[indxs[max_index],i,j], "ko")
			print("two maxima for F0=", F0s[i], " and D=", Ds[j], "for kappa=", kappas[indxs])
	#plt.plot(Ds, kappa_res[i, :], color="C"+str(i))
	#m, c, delta_m, delta_c = linear_regresion(np.log(Ds), np.log(kappa_res[i, :]))
	#plt.plot(Ds, np.exp(c+m*np.log(Ds)), color="C"+str(i))
	#slopes[i, 0] = m
	#slopes[i, 1] = delta_m
	#plt.title("F0=%3.2f" % F0s[i])
	#plt.xlabel("wave number")
	#plt.show()
#plt.plot(Ds, np.sqrt(1/Ds), "ro")
plt.xscale("log")
#plt.yscale("log")
plt.show()


for i in range(len(Ds)):
	plt.plot(F0s, kappa_res[:, i])
#plt.plot(Ds, np.sqrt(1/Ds), "ro")
plt.xscale("log")
plt.show()


fig = plt.figure(1)
x_, y_ = np.meshgrid(F0s, Ds)
Map = matplotlib.cm.get_cmap('Spectral_r')
kappa_res[np.where(kappa_res == 0)] = -1
ax1 = plt.contourf(x_,y_, np.transpose(((kappa_res))) , levels=np.linspace(0, np.amax(kappa_res), 15), cmap=Map)
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Resonance wave number $\kappa_{res}$', fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"External force  $F_0$", fontsize=8)
plt.ylabel(r"Diffusion coefficient $D$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/D_0_eff.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()






nu = 2.25
F0s      = np.logspace(0.5,  3, 25)/2.25
Ds       = np.logspace(-2, 1, 25) #np.logspace(-3, 0, 20) #np.arange(0.6, 2.7, 0.1)
kappas   = np.array([0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6, 4.4, 5.0, 6, 7.0, 8, 9]) #np.array([0.85, 0.90, 0.95, 1.0, 1.05, 1.10, 1.15])#np.arange(0.2, 1.75, 0.1)
D_parallels = np.load("data/D_para_vary_kappa_F0s_D_nu2.25_tau3.0.npy")

kapppa_cont = np.linspace(min(kappas), max(kappas), int(1e4))
nu = 1.25


kappa_res = np.zeros((len(F0s), len(Ds)))
slopes = np.zeros((len(F0s), 2))

for i in range(len(F0s)):
	for j in range(len(Ds)):
		interpoo = scp.interp1d(kappas, D_parallels[:, i, j], "cubic")(kapppa_cont)
		plt.plot(kappas, D_parallels[:, i, j])
		if len(scs.argrelmax(interpoo)[0]) == 1:
			kappa_res[i, j] = kapppa_cont[scs.argrelmax(interpoo)]
			#plt.plot(kappa_res[i, j], D_parallels[scs.argrelmax(D_parallels[:,i,j]),i,j], "ko")

		if len(scs.argrelmax(interpoo)[0]) > 1:
			indxs = scs.argrelmax(interpoo)[0]
			max_index = np.argmin(interpoo[indxs])
			kappa_res[i, j] = kapppa_cont[indxs[max_index]]
			#plt.plot(kappa_res[i, j], interpoo[indxs[max_index]], "ko")
			print("two maxima for F0=", F0s[i], " and D=", Ds[j], "for kappa=", kapppa_cont[indxs])
	#plt.plot(Ds, kappa_res[i, :], color="C"+str(i))
	#m, c, delta_m, delta_c = linear_regresion(np.log(Ds), np.log(kappa_res[i, :]))
	#plt.plot(Ds, np.exp(c+m*np.log(Ds)), color="C"+str(i))
	#slopes[i, 0] = m
	#slopes[i, 1] = delta_m
	#plt.title("F0=%3.2f" % F0s[i])
	#plt.xlabel("wave number")
	#plt.show()
#plt.plot(Ds, np.sqrt(1/Ds), "ro")
plt.xscale("log")
#plt.yscale("log")
plt.show()


for i in range(len(Ds)):
	plt.plot(F0s, kappa_res[:, i])
#plt.plot(Ds, np.sqrt(1/Ds), "ro")
plt.xscale("log")
plt.show()


fig = plt.figure(1)
x_, y_ = np.meshgrid(F0s, Ds)
Map = matplotlib.cm.get_cmap('Spectral_r')
kappa_res[np.where(kappa_res == 0)] = -1
ax1 = plt.contourf(x_,y_, np.transpose(((kappa_res))) , levels=np.linspace(0, np.amax(kappa_res), 15), cmap=Map)
cbar = fig.colorbar(ax1, format='%1.2f')
cbar.ax.set_ylabel(r'Resonance wave number $\kappa_{res}$', fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"External force  $F_0$", fontsize=8)
plt.ylabel(r"Diffusion coefficient $D$", fontsize=8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
#filename = root + "figures/D_0_eff.pdf"
#plt.savefig(filename, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()
