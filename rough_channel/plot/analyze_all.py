import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import os
import scipy.interpolate as sci 

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

root = "../master_latex/results/figures/"
plt.style.use(['science','no-latex', 'grid'])
#sns.color_palette("hls", 1)
import matplotlib
#matplotlib.rc('xtick', labelsize=10)
#matplotlib.rc('ytick', labelsize=10)
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'

#plt.style.use("bmh")
#sns.color_palette("hls", 1)

def func2(x, a, b):
    return a*np.exp(-b*x)

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


data_strings = ["data_square/Pe0_square_rough_mesh_res550_bmax0.4.dat",
                     "data_square/Pe0_square_rough_mesh_res380_bmax0.82.dat",
                     "data_square/Pe0_square_rough_mesh_res300_bmax1.24.dat",
                     "data_square/Pe0_square_rough_mesh_res200_bmax1.66.dat",
                     "data_square/Pe0_square_rough_mesh_res100_bmax1.98.dat"]

Brenner = np.zeros(int( (len(data_strings)-1)*21 + 16 + 1))
b = np.linspace(0, 2, int( (len(data_strings)-1)*21 + 16 + 1))

for i in range(len(data_strings)-1):
    data = np.loadtxt(data_strings[i])
    b[i*21:(i+1)*21]     = data[:, 0]
    Brenner[i*21:(i+1)*21] = data[:, 1]

data = np.loadtxt(data_strings[-1])
b[-17:-1]     = data[:, 0]
Brenner[-17:-1] = data[:, 1]

b[-1] = 2.0
Brenner[-1] = 0

interpool_Pe0 = sci.interp1d(b, Brenner)

cont_b = np.linspace(0, 2, int(1e4))
cont_b_short = np.linspace(0, 1, int(1e3))
m, c, delta_m, delta_c = linear_regresion(cont_b_short, interpool_Pe0(cont_b_short))


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


smoothd_Pe = np.logspace(0, 3, int(1e4))
Pe = np.logspace(0, 3, 40)
D_eff = np.zeros((len(data_strings), len(smoothd_Pe), 3))
b = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9])
colors_b = plt.cm.Spectral(np.linspace(0,1,int(len(b))))#jet is also good
max_index = np.zeros(len(data_strings)-1, dtype='int')
slopes_with_Pe = np.zeros((len(b), 4))
for i in range(len(data_strings)):
    data = np.loadtxt(data_strings[i])
    interpool = sci.interp1d(data[:, 0], data[:, 1], kind='cubic')

    Pe = data[:, 0]
    D_eff[i, :, 0] = (interpool(smoothd_Pe) - interpool_Pe0(b[i]))*105/(2*smoothd_Pe*smoothd_Pe)
    D_eff[i, :, 1] = interpool(smoothd_Pe)
    D_eff[i, :, 2] = (interpool(smoothd_Pe)-interpool_Pe0(b[i]))*105/(2*smoothd_Pe*smoothd_Pe)
    slopes_with_Pe[i, 0], slopes_with_Pe[i, 1], slopes_with_Pe[i, 2], slopes_with_Pe[i, 3] = linear_regresion(np.log(smoothd_Pe), np.log(D_eff[i,:,2]))

    if i % 3 == 0 and i > 0:
        fig = plt.figure(1)
        plt.plot(smoothd_Pe, D_eff[i, :, 2], label=r"$b=$"+str(b[i]), color=colors_b[i])
        #plt.plot(smoothd_Pe, np.exp(slopes_with_Pe[i,1]+np.log(smoothd_Pe)*slopes_with_Pe[i,0]), "--", color=sns.color_palette()[int(i/3)])

plt.xlabel(r"Peclet number Pe", fontsize=8)
plt.ylabel(r"Geometric factor $g$", fontsize=8)
plt.xscale("log")
#plt.axis([0.8, 1200, 0.8, np.amax(D_eff[:,:,2])*1.4])
plt.legend(loc="best", fontsize=8, ncol=3, labelspacing=0.2, borderpad=0.3, handlelength=0.8, framealpha=0.5)#edgecolor="k")#, frameon=False)
plt.yscale("log")

plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)

name = root+"g_tilde.pdf"
fig.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))

fig = plt.figure(2)
Peclet_analyze   = np.array([1, 3.16, 10, 31.6, 100, 316, 1000])
Peclet_exponents = [0, 0.5, 1, 1.5, 2, 2.5, 3]
colors = plt.cm.Spectral(np.linspace(0,1,len(Peclet_analyze)))#jet is also good

for i in range(len(Peclet_analyze)):
    current_data = D_eff[:,np.argmin(abs(Peclet_analyze[i]-smoothd_Pe)),2]
    plt.plot(b, current_data, label="Pe$=10^{%2.1f}$" % Peclet_exponents[i], color=colors[i])

plt.plot(b[:7], np.exp(2.85*b[:7]), "--", color="k", label=r"$e^{2.85b}$")
plt.xlabel(r"Rougness $b$", fontsize=8)
plt.ylabel(r"Geometric factor $g$", fontsize=8)

plt.legend(loc="upper left", fontsize=8, ncol=2, labelspacing=0.2, borderpad=0.3, handlelength=0.8, framealpha=0.5)#edgecolor="k")#, frameon=False)
plt.yscale("log")
plt.axis([-0.1, 1.99, 0.88, 93.8])
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=2, labelspacing=0.4, borderpad=0.6, handlelength=0.8, frameon=False)
name = root+"g_tilde_vs_b.pdf"
fig.savefig(name, bbox_inches="tight")
os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()

fig = plt.figure(113)
for i in range(len(Peclet_analyze)):
    current_data = D_eff[:,np.argmin(abs(Peclet_analyze[i]-smoothd_Pe)),1]
    D_aris = (1+2*Peclet_analyze[i]*Peclet_analyze[i]/105)
    plt.plot(b, (current_data-D_aris)/(D_aris), label="Pe$=10^{%2.1f}$" % Peclet_exponents[i])

plt.xlabel(r"Rougness $b$", fontsize=8)
plt.ylabel(r"Rel change in effective dispersion $\frac{D_\parallel-D_\parallel^{aris}}{D_\parallel^{aris}}$", fontsize=8)

plt.legend(loc="best", fontsize=8, ncol=1, labelspacing=0.2, borderpad=0.3, handlelength=0.8, framealpha=0.5)#edgecolor="k")#, frameon=False)

plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.legend(loc="best", fontsize=8, ncol=1, labelspacing=0.5, borderpad=0.8, handlelength=0.9)
name = "figures/rel_change_aris.pdf"
fig.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()


"""
print("For large b")
plt.figure(8)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
start_index = np.argmin(abs(b-0.3))
plt.errorbar(b, slopes_with_Pe[:,0], yerr=slopes_with_Pe[:,2], ls="none")
plt.scatter(b, slopes_with_Pe[:,0], s=5)

m, c, delta_m, delta_c = linear_regresion((b[start_index:]), slopes_with_Pe[start_index:, 0])
print("For slope: ", m, c, delta_m, delta_c)
plt.plot(b[start_index:], c+m*(b[start_index:]), label=r"$\beta = 0.06-0.211 b$")
#plt.xscale("log")
plt.xlabel(r"Rougness $b$", fontsize=8)
plt.legend(loc=3, fontsize=8)
plt.ylabel(r"Slope of geometric factor wrt. Pe $[\beta]$", fontsize=8)
plt.axis([-0.05, 2.0, -0.35, 0.02])
 # location for the zoomed portion 
sub_axes = plt.axes([0.7, 0.44, .2, .44]) 

# plot the zoomed portion
sub_axes.plot(b, slopes_with_Pe[:,1], c="k")
sub_axes.set_ylabel(r"Constant $\alpha$", fontsize=8)
sub_axes.set_xlabel(r"Roughness $b$", fontsize=8)
sub_axes.tick_params(axis='both', which='major', labelsize=8)
sub_axes.tick_params(axis='both', which='minor', labelsize=8)
sub_axes.axis("equal")

#sub_axes.xlabel(r"Rougness $b$", fontsize=14)
#name = "figures/slopePe_Re0.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))

plt.figure(8)
start_index = np.argmin(abs(b-1.1))
plt.errorbar(b, slopes_with_Pe[:,1], yerr=slopes_with_Pe[:,3], fmt="o")
m, c, delta_m, delta_c = linear_regresion((b[start_index:]), (slopes_with_Pe[start_index:, 1]))
print("For constant: ", m, c, delta_m, delta_c)
plt.plot(b[start_index:], c+m*(b[start_index:]), label=r"$\alpha = 2.890(3)+2.26(1)\log{(b)}$", color=sns.color_palette()[int(1)])
#plt.xscale("log")
plt.xlabel(r"Rougness $b$", fontsize=14)
plt.ylabel(r"Constant term of geometric factor wrt. Pe", fontsize=14)

#plt.legend(loc="best", fontsize=12)
#name = "figures/constantPe_Re0.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()

print("For small b")
plt.figure(8)
start_index = np.argmin(abs(b-0.4))
plt.errorbar(b, slopes_with_Pe[:,0], yerr=slopes_with_Pe[:,2], fmt="o")
plt.xlabel(r"Rougness $b$", fontsize=14)
plt.ylabel(r"Slope of geometric factor wrt. Pe", fontsize=14)

plt.figure(9)
plt.errorbar(b, slopes_with_Pe[:,1], yerr=slopes_with_Pe[:,3], fmt="o")

#late slope
start_index = np.argmin(abs(b-0.5))
m, c, delta_m, delta_c = linear_regresion(np.log(b[start_index:]), (slopes_with_Pe[start_index:, 1]))
print(m, c, delta_m, delta_c)
print("For constant: ", m, c, delta_m, delta_c)
plt.plot(b[start_index:], c+m*np.log(b[start_index:]), label=r"$\alpha = 2.91(1)+2.20(2)\log{(b)}$")

#early slope
end_index = np.argmin(abs(b-1.0))
start_index = np.argmin(abs(b-0.7))
m, c, delta_m, delta_c = linear_regresion((b[:end_index]), (slopes_with_Pe[:end_index, 1]))
print("For constant: ", m, c, delta_m, delta_c)
plt.plot(b[:end_index], c+m*(b[:end_index]), "--", label=r"$\alpha=-0.03(1)+2.99(3)b$", color=sns.color_palette()[3])

plt.xlabel(r"Rougness $b$", fontsize=14)
plt.ylabel(r"Constant term of geometric factor wrt. Pe", fontsize=14)
plt.legend(loc="best", fontsize=12)
#name = "figures/constant_of_Re0_smallb.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
#plt.show()

plt.figure(2)
for i in range(len(b)):
    if i < 6:
        plt.plot(smoothd_Pe, D_eff[i, :, 2], label=r"$b=$"+str(b[i]), color=sns.color_palette()[int(i)])
        #plt.plot(smoothd_Pe, np.exp(3*b[i]*np.ones(len(smoothd_Pe))), "--", color=sns.color_palette()[int(i)])

plt.xlabel(r"Peclet number Pe")
plt.ylabel(r"Geometric factor $g$")
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
#plt.yscale("log")
#name = "figures/g_tilde.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
#plt.show()


plt.figure(1)
for i in range(len(b)):
    if i > 10:
        plt.plot(smoothd_Pe, D_eff[i, :, 2], label=r"$b=$"+str(b[i]), color=sns.color_palette()[int(i-11)])
        plt.plot(smoothd_Pe, 18*b[i]**(9/4)*smoothd_Pe**(-(np.log(b[i])/4 +1/6)), "--", color=sns.color_palette()[int(i-11)])

plt.xlabel(r"Peclet number Pe")
plt.ylabel(r"Geometric factor $g$")
plt.xscale("log")
plt.legend(loc="best", fontsize=14)
#plt.yscale("log")
#name = "figures/g_tilde.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
#plt.show()
"""
"""
plt.figure(3)
m, c, delta_m, delta_c = linear_regresion( np.log(b[1:-5]), np.log(smoothd_Pe[max_index[:-5]]) )
plt.plot(b[1:-5], np.exp(c+m*np.log(b[1:-5])), label=r"Slope: "+str(m)[0:5] + "(" + str(delta_m*100)[:1] + ")")
print(m, delta_m)
m, c, delta_m, delta_c = linear_regresion( np.log(b[-5:]), np.log(smoothd_Pe[max_index[-5:]]) )
print(m, delta_m)#smoothd_Pe

plt.plot(b[-5:], np.exp(c+m*np.log(b[-5:])), label=r"Slope: "+str(m)[0:5] + "(" + str(delta_m*100)[:1] + ")")
plt.plot(b[1:], smoothd_Pe[max_index], "ko")
plt.yscale("log")
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
plt.ylabel(r"Optimal Peclet number Pe", fontsize=14)
plt.xlabel(r"Roughness $b$", fontsize=14)
#name = "figures/optimal_Pe.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()

"""
"""
fig = plt.figure(4)
x_, y_ = np.meshgrid(smoothd_Pe, b)
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
#plt.savefig("figures/g_tilde_contour.pdf")
"""
"""
plt.figure(5)
Peclet_analyze = np.array([1, 3.16, 10, 31.6, 100, 316, 1000])
Peclet_exponents = [0, 0.5, 1, 1.5, 2, 2.5, 3]
slopes = np.zeros((len(Peclet_analyze), 4))
max_index = np.argmin(abs(1.2-np.array(b)))
stop_index = max_index*np.ones(len(Peclet_analyze), dtype="int")# - np.linspace(0, len(Peclet_analyze)-1, len(Peclet_analyze), dtype="int")

for i in range(len(Peclet_analyze)):
    current_data = D_eff[:,np.argmin(abs(Peclet_analyze[i]-smoothd_Pe)),2]
    plt.plot(b, current_data, label="Pe = $10^{%2.1f}$" % Peclet_exponents[i] )
    slopes[i, 0], slopes[i,2], slopes[i,1], slopes[i,3] = linear_regresion(b[:stop_index[i]-i], np.log(current_data[:stop_index[i]-i]))
    print(b[:stop_index[i]-i])
    #slopes[i, 0], c, slopes[i,1], delta_c = linear_regresion(b, np.log(current_data))

plt.yscale("log")
plt.legend(loc="best", fontsize=12)
plt.xlabel(r"Roughness $b$", fontsize=14)
plt.ylabel(r"Geometric factor $g$", fontsize=14)
#name = "figures/gtilde_vs_b_Re0.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))

plt.figure(6)
plt.errorbar(Peclet_analyze, slopes[:, 0], yerr=slopes[:,1], fmt="o")
plt.xscale("log")
plt.xlabel(r"Peclet number Pe", fontsize=14)
plt.ylabel(r"Slope of geometric factor $\frac{\partial g}{\partial b}$", fontsize=14)
m, c, delta_m, delta_c = linear_regresion(np.log(Peclet_analyze), slopes[:,0])
print(m, c, delta_m, delta_c)
plt.plot(Peclet_analyze, c+m*np.log(Peclet_analyze))
#name = "figures/slopes_Re0.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))


plt.figure(8)
plt.errorbar(Peclet_analyze, slopes[:, 2], yerr=slopes[:,3], fmt="o")
plt.xscale("log")
plt.xlabel(r"Peclet number Pe", fontsize=14)
plt.ylabel(r"Slope of geometric factor $\frac{\partial g}{\partial b}$", fontsize=14)
#m, c, delta_m, delta_c = linear_regresion(np.log(Peclet_analyze), slopes[:,0])
#print(m, c, delta_m, delta_c)
#plt.plot(Peclet_analyze, c+m*np.log(Peclet_analyze))
#name = "figures/slopes_Re0.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))

plt.figure(7)
for i in range(len(data_strings)):
    if i < 9:
        plt.plot(smoothd_Pe, D_eff[i, :, 2], label=r"$b=$"+str(b[i]), color=sns.color_palette()[i])
        plt.plot(smoothd_Pe, smoothd_Pe**(m*b[i])*np.exp(c*b[i]), "--", color=sns.color_palette()[i])

plt.xlabel(r"Peclet number Pe")
plt.ylabel(r"Geometric factor $g$(Re $=0$)")
plt.xscale("log")
plt.legend(loc="best", fontsize=12)
#name = "figures/g_tilde.pdf"
#plt.savefig(name, bbox_inches="tight")
#os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show(   )
"""