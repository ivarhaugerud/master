import numpy as np 
import matplotlib.pyplot as plt 

N = int(1e4)
t = np.linspace(0, 1, N)
omega = 4*2*np.pi

arg = omega*t 

cosinus =  np.cos(arg)
sinus   =  np.sin(arg)

fracs_to_test = (0.75+np.linspace(0, 3, 4))/4
experimental = [0.22, 0.47, 0.72, 0.968]
for i in range(len(fracs_to_test)):
	plt.title(str(fracs_to_test[i]))
	plt.plot(t[:int(fracs_to_test[i]*N)], cosinus[:int(fracs_to_test[i]*N)])
	plt.plot(t[:int(fracs_to_test[i]*N)], -np.flipud(sinus[:int(fracs_to_test[i]*N)]), "--")
	plt.show()