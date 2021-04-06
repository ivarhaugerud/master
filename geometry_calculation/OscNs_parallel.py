from joblib import Parallel, delayed
import numpy as np 
import os 

epsilon = np.logspace(-2, np.log10(0.6), 8)
kappa = 0.5
Lx = 2*np.pi/kappa

def processInput(i, eps):
	res = int(150*(1+2*float(eps)))
	#print('mpiexec -n 4 python3 OscNS.py -res ' + str(res) + ' -dt 0.02 -T 25 -tau 5 -nu 1 -D 0.3 -f0 0 -f1 3 -Lx ' + i + " -epsilon " + str(eps))
	#res = int(25+np.sqrt(2*150*150/i))
	#print(int(i*res/2), res)
	#print('python3 OscNS.py -res ' + str(res) + ' -dt 0.02 -T 25 -tau 5 -nu 1 -D 0.3 -f0 0 -f1 3 -Lx ' + str(i)[:7] + " -epsilon " + eps)
	#os.system('python3 OscNS.py -res 200 -dt 0.01 -T 1 -tau 0.5 -Lx ' + str(i) + "-epsilon " + str(eps))
	os.system('python3 OscNS.py -res ' + str(res) + ' -dt 0.02 -T 25 -tau 5 -nu 1 -D 0.3 -f0 0 -f1 3 -Lx ' + i + " -eps " + str(eps))

num_cores = 4

Parallel(n_jobs=num_cores)(delayed(processInput)(str(Lx)[:6], str(e)[:5]) for e in epsilon)
