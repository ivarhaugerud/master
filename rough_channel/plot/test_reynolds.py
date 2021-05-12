import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
import os
import scipy.interpolate as sci 

plt.style.use(['science','no-latex', 'grid'])
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
root = "../master_latex/results/figures/rough/"

datafiles = {}
datafiles["Re0"] = ["data_square/Re0.0_b0.0_res700.dat",
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
                    "data_square/Re0.0_b1.4_res200.dat",
                    "data_square/Re0.0_b1.5_res200.dat",
                    "data_square/Re0.0_b1.6_res210.dat",
                    "data_square/Re0.0_b1.7_res240.dat"]

datafiles["Re0.01"] = ["data_square/Re0.010_b0.0_res400.dat",
		               "data_square/Re0.009_b0.1_res300.dat",
		               "data_square/Re0.009_b0.2_res200.dat",
		               "data_square/Re0.009_b0.3_res300.dat",
		               "data_square/Re0.009_b0.4_res300.dat",
		               "data_square/Re0.009_b0.5_res280.dat",
		               "data_square/Re0.009_b0.6_res260.dat",
		               "data_square/Re0.009_b0.7_res240.dat",
                   "data_square/Re0.009_b0.8_res220.dat",
                   "data_square/Re0.009_b0.9_res210.dat",
                   "data_square/Re0.009_b1.0_res200.dat",
                   "data_square/Re0.009_b1.1_res195.dat",
                   "data_square/Re0.009_b1.2_res190.dat",
                   "data_square/Re0.009_b1.3_res210.dat",
                   "data_square/Re0.009_b1.4_res200.dat",
                   "data_square/Re0.009_b1.5_res200.dat",
                   "data_square/Re0.009_b1.6_res210.dat",
                   "data_square/Re0.009_b1.7_res210.dat"]

datafiles["Re1"]   = ["data_square/Re0.999_b0.0_res400.dat",
		                  "data_square/Re0.999_b0.1_res300.dat",
		                  "data_square/Re0.999_b0.2_res200.dat",
		                  "data_square/Re0.999_b0.3_res300.dat",
		                  "data_square/Re0.999_b0.4_res300.dat",
		                  "data_square/Re0.986_b0.5_res280.dat",
		                  "data_square/Re0.999_b0.6_res260.dat",
		                  "data_square/Re0.999_b0.7_res240.dat",
                      "data_square/Re0.999_b0.8_res220.dat",
                      "data_square/Re0.999_b0.9_res210.dat",
                      "data_square/Re0.999_b1.0_res200.dat",
                      "data_square/Re0.999_b1.1_res195.dat",
                      "data_square/Re0.999_b1.2_res190.dat",
                      "data_square/Re1.014_b1.3_res200.dat",
                      "data_square/Re0.999_b1.4_res210.dat",
                      "data_square/Re0.999_b1.5_res210.dat",
                      "data_square/Re0.998_b1.6_res210.dat",
                      "data_square/Re0.999_b1.7_res210.dat"]

datafiles["Re3.16"]  = ["data_square/Re3.162_b0.0_res700.dat",
                       "data_square/Re3.162_b0.1_res700.dat",
                       "data_square/Re3.215_b0.2_res300.dat",
                       "data_square/Re3.160_b0.3_res400.dat",
                       "data_square/Re3.162_b0.4_res400.dat",
                       "data_square/Re3.162_b0.5_res375.dat",
                       "data_square/Re3.162_b0.6_res325.dat",
                       "data_square/Re3.161_b0.7_res300.dat",
                       "data_square/Re3.159_b0.8_res300.dat",
                       "data_square/Re3.160_b0.9_res250.dat",
                       "data_square/Re3.160_b1.0_res250.dat",
                       "data_square/Re3.160_b1.1_res250.dat",
                       "data_square/Re3.158_b1.2_res250.dat",
                       "data_square/Re3.156_b1.3_res250.dat",
                       "data_square/Re3.143_b1.4_res225.dat",
                       "data_square/Re3.144_b1.5_res215.dat",
                       "data_square/Re3.161_b1.6_res215.dat",
                       "data_square/Re3.134_b1.7_res210.dat"]

datafiles["Re10"]   = ["data_square/Re10.05_b0.0_res700.dat",
                       "data_square/Re9.999_b0.1_res400.dat",
                       "data_square/Re9.999_b0.2_res400.dat",
                       "data_square/Re10.25_b0.3_res370.dat",
                       "data_square/Re9.999_b0.4_res320.dat",
                       "data_square/Re9.997_b0.5_res280.dat",
                       "data_square/Re10.00_b0.6_res260.dat",
                       "data_square/Re10.24_b0.7_res240.dat",
                       "data_square/Re10.00_b0.8_res220.dat",
                       "data_square/Re9.978_b0.9_res200.dat",
                       "data_square/Re9.964_b1.0_res190.dat",
                       "data_square/Re9.989_b1.1_res180.dat",
                       "data_square/Re9.916_b1.2_res180.dat",
                       "data_square/Re9.469_b1.3_res210.dat",
                       "data_square/Re9.408_b1.4_res200.dat",
                       "data_square/Re9.363_b1.5_res200.dat",
                       "data_square/Re10.62_b1.6_res210.dat",
                       "data_square/Re9.472_b1.7_res210.dat"]

datafiles["Re31.6"]   = ["data_square/Re31.62_b0.0_res700.dat",
                         "data_square/Re31.62_b0.1_res700.dat",
                         "data_square/Re31.62_b0.2_res700.dat",
                         "data_square/Re31.62_b0.3_res500.dat",
                         "data_square/Re31.62_b0.4_res450.dat",
                         "data_square/Re31.61_b0.5_res400.dat",
                         "data_square/Re31.61_b0.6_res350.dat",
                         "data_square/Re31.52_b0.7_res350.dat",
                         "data_square/Re31.55_b0.8_res325.dat",
                         "data_square/Re31.49_b0.9_res300.dat",
                         "data_square/Re31.49_b1.0_res275.dat",
                         "data_square/Re31.42_b1.1_res250.dat",
                         "data_square/Re31.10_b1.2_res250.dat",
                         "data_square/Re31.19_b1.3_res230.dat",
                         "data_square/Re31.01_b1.4_res210.dat",
                         "data_square/Re31.09_b1.5_res210.dat",
                         "data_square/Re30.18_b1.6_res240.dat",
                         "data_square/Re31.43_b1.7_res230.dat"]

datafiles["Re100"] = ["data_square/Re98.69_b0.0_res700.dat",
			                "data_square/Re99.98_b0.1_res300.dat",
			                "data_square/Re99.88_b0.2_res200.dat",
			                "data_square/Re99.16_b0.3_res150.dat",
			                "data_square/Re98.80_b0.4_res150.dat",
			                "data_square/Re97.12_b0.5_res280.dat",
			                "data_square/Re99.91_b0.6_res120.dat",
                      "data_square/Re96.01_b0.7_res240.dat",
                      "data_square/Re99.88_b0.8_res220.dat",
                      "data_square/Re99.80_b0.9_res210.dat",
                      "data_square/Re99.69_b1.0_res200.dat",
                      "data_square/Re99.70_b1.1_res195.dat",
                      "data_square/Re99.33_b1.2_res190.dat",
                      "data_square/Re99.34_b1.3_res200.dat",
                      "data_square/Re99.12_b1.4_res200.dat"]

datafiles["Re200"]   = ["data_square/Re200.0_b0.0_res700.dat",
                        "data_square/Re199.9_b0.1_res400.dat",
                        "data_square/Re199.6_b0.2_res400.dat",
                        "data_square/Re199.9_b0.3_res370.dat",
                        "data_square/Re199.9_b0.4_res320.dat",
                        "data_square/Re199.9_b0.5_res280.dat",
                        "data_square/Re199.9_b0.6_res260.dat",
                        "data_square/Re199.8_b0.7_res240.dat",
                        "data_square/Re199.7_b0.8_res220.dat",
                        "data_square/Re200.2_b0.9_res200.dat",
                        "data_square/Re200.1_b1.0_res230.dat",
                        "data_square/Re202.0_b1.1_res220.dat",
                        "data_square/Re201.8_b1.2_res210.dat",
                        "data_square/Re199.7_b1.3_res200.dat"]

Re = [0, 1, 3.16, 10, 31.6, 100]
Pe = np.logspace(0, 3, int(1e4))
b = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]

D_eff_Pe0 = interpool_Pe0(b)
D = np.zeros( (len(Re), len(b), len(Pe)))
g_tilde = np.zeros( (len(Re), len(b), len(Pe)))
max_Pe = np.zeros((len(Re), len(b)))

for j in range(len(Re)):
    for i in range(len(b)):
        try:
          data = np.loadtxt(datafiles["Re"+str(Re[j])][i])
          current_b = np.linspace(0, (len(data)-1)/10, len(data))

          interpool = sci.interp1d(data[:, 0], data[:, 1], kind='cubic')
          D[j, i, :] = interpool(Pe)
          g_tilde[j, i, :]    = (interpool(Pe)-interpool_Pe0(b[i]))*105/(2*Pe*Pe)#105*(D[j, i, :]-1)/(2*Pe*Pe)

        except:
          a = 1

Pe_analyze = [1, 3.16, 10, 31.6, 100, 316, 1000]
Pe_analyze_indecs = np.zeros(len(Pe_analyze), dtype="int")

for i in range(len(Pe_analyze)):
    Pe_analyze_indecs[i] = np.argmin(abs(Pe-Pe_analyze[i]))

analyze_bs = [0, 3, 6, 9, 12]

###
#plot of g tilde vs Pe for different Re numbers for specific b
###
"""

for i in range(len(b)):
    for j in range(len(Re)-1):
        plt.figure(i)
        plt.plot(Pe, g_tilde[j+1, i, :], color=sns.color_palette()[j], label=r"Re = "+str(Re[j+1]))
    
    plt.title(r"Roughness $b=$"+str(b[i]))
    plt.xscale("log")
    #plt.yscale("log")
    plt.legend(loc="best")
    plt.ylabel(r"Geometric factor $\tilde{g}$")
    plt.xlabel(r"Peclet number Pe")
    plt.savefig("figures/g_vs_Pe_different_Re_for_b"+str(b[i])+".pdf")

plt.show()

for j in range(len(b)):
    if j % 3 == 0:
      plt.figure(j)
      for i in range(len(Re)):
          plt.plot(Pe, (D[i, j, :]-D[0, j, :])/D[0, j, :], "-", label=r"$Re=$"+str(Re[i]))
      plt.xscale("log")
      plt.xlabel(r"Peclet number Pe", fontsize=14)
      plt.ylabel(r"Relative change in $D_{\parallel}$  $\left[\frac{ D_{\parallel}(Re)- D_{\parallel}(Re=0)}{D_{\parallel}(Re=0)}\right]$", fontsize=14)
      plt.legend(loc="best", fontsize=12)
      plt.title(r"Roughness $b=$" + str(b[j]), fontsize=14)
      name = "figures/rel_change_in_g_each_b_vs_Pe_for_b"+str(b[j]*10)[:1]+".pdf"
      plt.savefig(name, bbox_inches="tight")
      os.system('pdfcrop %s %s &> /dev/null &'%(name, name))
plt.show()
"""
for i in range(len(Pe_analyze)):
  for j in range(len(Re)-1):
    j += 1
    fig = plt.figure(i)
    plt.title("Peclet number Pe $=$ " + str(Pe_analyze[i]), fontsize=8)
    current_len_b = len(np.trim_zeros(g_tilde[j,:,Pe_analyze_indecs[i]]))
    interpool_b = np.linspace(0, b[current_len_b]-0.1, int(1e4))
    f = (D[j, :current_len_b, Pe_analyze_indecs[i]]-D[0, :current_len_b, Pe_analyze_indecs[i]])/D[0, :current_len_b, Pe_analyze_indecs[i]]
    interpool = sci.interp1d(b[:current_len_b], f, kind='cubic')
    plt.plot(interpool_b, interpool(interpool_b), label="Re = " + str(Re[j]))
  
  plt.ylabel(r"Rel. change in $D_{\parallel}$  $\left[ \frac{D_{\parallel}(Pe, Re, b) - D_{\parallel}(Pe, 0, b)}{D_{\parallel}(Pe, 0, b)}  \right]$", fontsize=8)
  plt.xlabel(r"Roughness $b$", fontsize=8)
  plt.legend(loc="best", fontsize=8)
  filename = "rel_change_vs_b_for_Pe=" + str(Pe_analyze[i])
  filename = filename.replace(".", "_")
  filename = root + filename
  filename += ".pdf"
  plt.tick_params(axis='both', which='major', labelsize=8)
  plt.tick_params(axis='both', which='minor', labelsize=8)
  plt.savefig(filename, bbox_inches="tight")
  os.system('pdfcrop %s %s &> /dev/null &'%(filename, filename))
plt.show()

###
#plot g-tilde vs Pe for different b's for a specific Reynolds-number 
###
for i in range(len(Re)):
    for j in range(len(b)):
        if j%4 == 0:
            plt.figure(i)
            plt.plot(Pe, (D[i, j, :]-D_eff_Pe0[j])/(D[i, 0, :]-D_eff_Pe0[j]), label=r"$b=$"+str(b[j]))

    plt.figure(i)
    plt.xscale("log")
    plt.xlabel(r"Peclet number Pe", fontsize=14)
    plt.ylabel(r"Change in $D_{\parallel}$  $\left[ \frac{D_{eff}(Pe, Re, b) - D_{eff}(0, 0, b)}{D_{eff}(Pe, Re, 0) - D_{eff}(0, 0, 0)}  \right]$", fontsize=14)
    plt.legend(loc="best", fontsize=12)
    plt.title(r"Reynoldsnumber Re $=$" + str(Re[i]), fontsize=14)
    plt.savefig("figures/relative_change_for_Re"+str(Re[i])+".pdf")

plt.show()