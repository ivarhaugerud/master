import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as sci
from scipy.interpolate import griddata

#general variabels for all runs
#U  = 1
#dx = 1e-4
#dt = dx/U
datafiles = {}
Np = 1000

datafiles["Re0"] = ["data_square/U_Re0.0_b0.0.h5",
                    "data_square/U_Re0.0_b0.1.h5",
                    "data_square/U_Re0.0_b0.2.h5",
                    "data_square/U_Re0.0_b0.3.h5",
                    "data_square/U_Re0.0_b0.4.h5",
                    "data_square/U_Re0.0_b0.5.h5",
                    "data_square/U_Re0.0_b0.6.h5",
                    "data_square/U_Re0.0_b0.7.h5",
                    "data_square/U_Re0.0_b0.8.h5",
                    "data_square/U_Re0.0_b0.9.h5",
                    "data_square/U_Re0.0_b1.0.h5",
                    "data_square/U_Re0.0_b1.1.h5",
                    "data_square/U_Re0.0_b1.2.h5",
                    "data_square/U_Re0.0_b1.3.h5",
                    "data_square/U_Re0.0_b1.4.h5",
                    "data_square/U_Re0.0_b1.5.h5",
                    "data_square/U_Re0.0_b1.6.h5",
                    "data_square/U_Re0.0_b1.7.h5",
                    "data_square/U_Re0.0_b1.8.h5",
                    "data_square/U_Re0.0_b1.9.h5"]

datafiles["Re3.16"]  = ["data_square/U_Re3.1622992954684266_b0.0.h5",
					   "data_square/U_Re3.1622824459933745_b0.1.h5",
					   "data_square/U_Re3.1743667288554556_b0.2.h5",
					   "data_square/U_Re3.1389545179169165_b0.3.h5",
					   "data_square/U_Re3.16289307602155_b0.4.h5",
					   "data_square/U_Re3.162250167750368_b0.5.h5",
					   "data_square/U_Re3.1620604605281137_b0.6.h5",
					   "data_square/U_Re3.1613413262353163_b0.7.h5",
					   "data_square/U_Re3.1599334969768544_b0.8.h5",
					   "data_square/U_Re3.1599334969768544_b0.8.h5",
					   "data_square/U_Re3.1603960417556274_b0.9.h5",
					   "data_square/U_Re3.160744387020638_b1.0.h5",
					   "data_square/U_Re3.160318000306653_b1.1.h5",
					   "data_square/U_Re3.1584782127278257_b1.2.h5",
					   "data_square/U_Re3.1565632191172868_b1.3.h5",
					   "data_square/U_Re3.143656837514886_b1.4.h5",
					   "data_square/U_Re3.144955122968455_b1.5.h5",
					   "data_square/U_Re3.13089885490002_b1.6.h5",
					   "data_square/U_Re3.1344305866432576_b1.7.h5"]

datafiles["Re200"]=["data_square/U_Re200.0000694519953_b0.0.h5",
					"data_square/U_Re199.94372793156856_b0.1.h5",
					"data_square/U_Re199.61049255447293_b0.2.h5",
					"data_square/U_Re198.83064666477193_b0.3.h5",
					"data_square/U_Re199.98499284126197_b0.4.h5",
                    "data_square/U_Re199.96463226764294_b0.5.h5",
                    "data_square/U_Re199.9288531368267_b0.6.h5",
                    "data_square/U_Re199.86895377643867_b0.7.h5",
                    "data_square/U_Re199.7758564428973_b0.8.h5",
                    "data_square/U_Re200.20258987512474_b0.9.h5",
					"data_square/U_Re200.14729497191382_b1.0.h5",
                    "data_square/U_Re202.06420758993707_b1.1.h5",
                    "data_square/U_Re201.87980508787433_b1.2.h5",
                    "data_square/U_Re199.75965102299836_b1.3.h5"]

datafiles["Re100"]=["data_square/U_Re99.99516550653246_b0.0.h5",
					"data_square/U_Re99.98722541016618_b0.1.h5",
					"data_square/U_Re99.8891039092023_b0.2.h5",
					"data_square/U_Re99.63282436417836_b0.3.h5",
					"data_square/U_Re99.14138058731952_b0.4.h5",
                    "data_square/U_Re97.12638496567374_b0.5.h5",
                    "data_square/U_Re100.21614117614246_b0.6.h5",
                    "data_square/U_Re96.01492224789249_b0.7.h5",
                    "data_square/U_Re99.8855949408189_b0.8.h5",
                    "data_square/U_Re99.80651625648649_b0.9.h5",
					"data_square/U_Re99.69198187524063_b1.0.h5",
                    "data_square/U_Re99.7061142188546_b1.1.h5",
                    "data_square/U_Re99.33730855250279_b1.2.h5",
                    "data_square/U_Re99.34428447407375_b1.3.h5",
                    "data_square/U_Re99.1296244009446_b1.4.h5"]

datafiles["Re10"]=[ "data_square/U_Re10.051804966552208_b0.0.h5",
					"data_square/U_Re9.999997770621917_b0.1.h5",
					"data_square/U_Re9.999930337643791_b0.2.h5",
					"data_square/U_Re10.255257538646818_b0.3.h5",
					"data_square/U_Re9.771864894326647_b0.4.h5",
                    "data_square/U_Re9.997659266298484_b0.5.h5",
                    "data_square/U_Re10.001642453706866_b0.6.h5",
                    "data_square/U_Re10.247256841667147_b0.7.h5",
                    "data_square/U_Re10.004834460729734_b0.8.h5",
                    "data_square/U_Re9.978422162880927_b0.9.h5",
					"data_square/U_Re9.96448746834171_b1.0.h5",
                    "data_square/U_Re9.989318864931873_b1.1.h5",
                    "data_square/U_Re9.916423413075721_b1.2.h5",
                    "data_square/U_Re9.98573547811867_b1.3.h5",
                    "data_square/U_Re9.408779838364293_b1.4.h5",
                    "data_square/U_Re9.363471196840745_b1.5.h5",
                    "data_square/U_Re10.62267836105853_b1.6.h5",
             		"data_square/U_Re9.472213212753651_b1.7.h5"]

datafiles["Re1"] = ["data_square/U_Re0.9999999816977955_b0.0.h5",
					"data_square/U_Re0.9999999931630003_b0.1.h5",
					"data_square/U_Re0.9999996844273361_b0.2.h5",
					"data_square/U_Re0.9999973075406687_b0.3.h5",
					"data_square/U_Re0.9999873849533011_b0.4.h5",
                    "data_square/U_Re0.9869973750262079_b0.5.h5",
                    "data_square/U_Re0.9999104949298099_b0.6.h5",
                    "data_square/U_Re0.9998381202804935_b0.7.h5",
                    "data_square/U_Re0.9997425912249307_b0.8.h5",
                    "data_square/U_Re0.9996180157388926_b0.9.h5",
					"data_square/U_Re0.9994614176461558_b1.0.h5",
                    "data_square/U_Re0.9992814521837956_b1.1.h5",
                    "data_square/U_Re0.9990974203048193_b1.2.h5",
                    "data_square/U_Re1.014350817085215_b1.3.h5",
                    "data_square/U_Re0.9994783032179493_b1.4.h5",
                    "data_square/U_Re0.9987890132974687_b1.5.h5",
                    "data_square/U_Re0.9988508601151445_b1.6.h5",
                    "data_square/U_Re0.9990163022468528_b1.7.h5"]

datafiles["Re001"]=["data_square/U_Re0.010000151146402576_b0.0.h5",
					"data_square/U_Re0.0099999999979114_b0.1.h5",
					"data_square/U_Re0.009999999998958571_b0.2.h5",
					"data_square/U_Re0.009999999996608364_b0.3.h5",
					"data_square/U_Re0.009999999990544274_b0.4.h5",
                    "data_square/U_Re0.009870353581223542_b0.5.h5",
                    "data_square/U_Re0.009999999979643247_b0.6.h5",
                    "data_square/U_Re0.009999999963368226_b0.7.h5",
                    "data_square/U_Re0.009999999921646041_b0.8.h5",
                    "data_square/U_Re0.009999999901382194_b0.9.h5",
					"data_square/U_Re0.009999999881619656_b1.0.h5",
                    "data_square/U_Re0.009999999290071534_b1.1.h5",
                    "data_square/U_Re0.00999999970433836_b1.2.h5",
                    "data_square/U_Re0.010591256566707232_b1.3.h5",
                    "data_square/U_Re0.009999998811874288_b1.4.h5",
                    "data_square/U_Re0.009999998780497767_b1.5.h5",
                    "data_square/U_Re0.010000000712796187_b1.6.h5",
                    "data_square/U_Re0.009999998756266462_b1.7.h5"]

names = ["Re0", "Re001", "Re1", "Re3.16", "Re10", "Re100", "Re200"]
Reynolds = [0, 0.01, 1, 3.16, 10, 100, 200]
bs = np.linspace(0, 1.9, 20)

prop_RZs = np.zeros((len(Reynolds), len(bs)))

for c in range(len(names)):
	for a in range(len(datafiles[names[c]])):
		a += 15
		f = h5py.File(datafiles[names[c]][a], 'r')
		vector = np.array(list(f["VisualisationVector"]["0"]))
		geometry = np.array(list(f["Mesh"]["mesh"]["geometry"]))
		Re = Reynolds[c]
		print(Reynolds[c], bs[a])
		non_RZ = []
		RZ = []
		
		b = 0.5*max(geometry[:,0])
		if  a == 0:
			b = 0

		for i in range(len(geometry[:,0])):
			if -1 < geometry[i, 1] < -1+b/2 and 1e-3 < geometry[i,0] < b/2 - 1e-3:
				if vector[i, 1] > 0:
					non_RZ.append(i)
				else:
					RZ.append(i)

			elif -1 < geometry[i, 1] < -1+b/2 and 2*b-1e-3 > geometry[i,0] > 3*b/2 + 1e-3:
				if vector[i, 1] < 0:
					non_RZ.append(i)
				else:
					RZ.append(i)
			elif geometry[i, 1] < -1+b/2 and geometry[i,0] > 1e-3 and (geometry[i,0] < 0.5*b-1e-3 or 3*b/2 + 1e-3 < geometry[i,0] < 2*b-1e-3):
				RZ.append(i)


		#plt.scatter(geometry[:,0], geometry[:,1])	
		if a == 0:
			prop_RZ = 0
			#print(b, prop_RZ, len(non_RZ), len(RZ), (1) )

		else:
			prop_RZ = b*(1-len(non_RZ)/len(RZ))/2
			A = len(RZ)/(len(RZ)+len(non_RZ))
			new_prop = b*(1+A)/4
			test = len(RZ)/len(geometry[:,0])
			prop_RZ = new_prop
			#print("NEW PROP = ", new_prop, a)
			print(b, Re, prop_RZ, len(non_RZ), len(RZ))

		if a == 0:
			b = 0.0

		#np.save("analyzed_data/non_RZ_Re"+str(Re)+"_b"+str(b)[:3], geometry[non_RZ,:])
		#np.save("analyzed_data/RZ_Re"+str(Re)+"_b"+str(b)[:3],     geometry[RZ,:])
		#for i in range(len(non_RZ)):
		#	plt.plot(geometry[non_RZ[i], 0], geometry[non_RZ[i], 1], "ro")
		
		
		x_axis = geometry[:,0]
		y_axis = geometry[:,1]

		Nx = int(1e3)
		Ny = int(1e3)
		x = np.linspace(min(x_axis), max(x_axis), Nx)
		y = np.linspace(min(y_axis), max(y_axis), Ny)
		X, Y = np.meshgrid(x,y)

		# Interpolate using three different methods and plot
		ux = griddata((x_axis, y_axis), vector[:,0], (X, Y), method='cubic')
		uy = griddata((x_axis, y_axis), vector[:,1], (X, Y), method='cubic')
		speed = np.sqrt(np.square(ux) + np.square(uy))
		CS = plt.contourf(X, Y, speed, 10)
		
		plt.scatter(geometry[RZ,0], geometry[RZ, 1], color="r", alpha=0.5)
		plt.streamplot(X, Y, ux, uy, density=0.6+0.3, color='k', arrowsize=0.6)
		plt.show()
		#print(c, a, prop_RZ, np.shape(prop_RZs))
		prop_RZs[c, a] = prop_RZ 
		
#np.save("analyzed_data/RZ_prop_try2", prop_RZs)
for i in range(len(prop_RZs[:,0])):
	plt.plot(bs, prop_RZs[i,:])
plt.show()
