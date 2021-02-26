import numpy as np 
import matplotlib.pyplot as plt 

#simulation parameters
dt           = 0.01
tau          = 3.0 
periods      = 50
timesteps    = int(tau/dt)
datafiles    = 1000
N            = int(15*1e3)
RW_timesteps = 20
D            = 0#1e-6
alpha        = np.sqrt(2*D*dt/RW_timesteps)
var          = np.zeros(int(periods*timesteps))
t            = np.linspace(0, tau-dt, timesteps)+tau/2
pos_saves    = np.linspace(0, timesteps*periods, datafiles, dtype="int")
saves_counter= 0
kappa = 1.05
l = 2*np.pi/kappa
epsilon = 0.1 

eta = np.linspace(-1*np.pi/kappa, 3*np.pi/kappa, int(1e3))
B_up =   1+epsilon*np.sin(kappa*eta)
B_down =-1-epsilon*np.sin(kappa*eta)

dirr = "data/analytic_D0_eps"+str(epsilon)+"_kappa"+str(kappa)






from matplotlib.collections import LineCollection

t = np.linspace(0,1, N) # your "time" variable



#colors = cm.rainbow(t)

for i in range(datafiles):
	pos = np.load(dirr + "/pos_"+str(pos_saves[i])+".npy")
	plt.clf()
	"""

	# set up a list of (x,y) points
	points = np.array([pos[0,:],pos[1,:]]).transpose().reshape(-1,1,2)
	#print points.shape  # Out: (len(x),1,2)

	# set up a list of segments
	segs = np.concatenate([points[:-1],points[1:]],axis=1)
	#print segs.shape  # Out: ( len(x)-1, 2, 2 )
	                  # see what we've done here -- we've mapped our (x,y)
	                  # points to an array of segment start/end coordinates.
	                  # segs[i,0,:] == segs[i-1,1,:]
	print(np.shape(segs))

	# make the collection of segments
	norm = plt.Normalize(0, 1)
	lc = LineCollection(t, cmap=plt.get_cmap('viridis'), norm=norm)
	lc.set_array(t) # color the segments by our parameter

	# plot the collection
	#plt.gca().add_collection(lc) # add the collection to the plot
	plt.xlim(eta.min(), eta.max()) # line collections don't auto-scale the plot
	plt.ylim(-1-epsilon, 1+epsilon)
	"""
	plt.scatter(pos[0, :], pos[1, :], c=t, cmap="plasma", s=0.5)
	plt.plot(eta, B_up, "k")
	plt.plot(eta, B_down, "k")
	plt.xlabel(r"Horinzontal position [$a$]", fontsize=8)
	plt.ylabel(r"Vertical position [$a$]", fontsize=8)
	plt.pause(0.01)

