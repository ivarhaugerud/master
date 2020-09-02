import numpy as np 
import matplotlib.pyplot as plt 


b = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
Re = [0]

RZ_area = np.zeros((len(Re), len(b)))
for j in range(len(Re)):
	for i in range(len(b)):
		pos = np.load("analyzed_data/RZ_b="+str(b[i])+"_Re="+str(Re[j])+".npy")

		nr_of_converged = 0
		Np = int(len(pos[0,:,0]))

		plt.scatter(pos[0,:,0], pos[0,:,1])
		for k in range(Np):
			if abs(pos[-1, k, 0] - pos[0,k,0]) < 2*b[i]:
				nr_of_converged += 1
		RZ_area[j,i] = nr_of_converged/Np

plt.plot(b, RZ_area, "o")
plt.show()
		
