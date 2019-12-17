import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class LBM():
    def __init__(self, Nx, Ny, tau, F, tolerence):
        #size and time
        self.Nx = Nx
        self.Ny = Ny
        self.epsilon = tolerence

        #always true for D2Q9
        self.occupations = 9
        occupations      = self.occupations
        self.delta_t     = 1
        self.counter     = 0

        #to save data
        self.velocity = np.zeros((Ny, Nx, 2))
        self.density  = np.zeros((Ny, Nx))

        #occupations, velocity, density, source and force
        self.f        = np.zeros((Ny, Nx, occupations))
        self.f_prev   = np.zeros((Ny, Nx, occupations))
        self.f_equil  = np.zeros((Ny, Nx, occupations))
        self.u        = np.zeros((Ny ,Nx, 2))
        self.rho      = np.zeros((Ny, Nx))
        self.S        = np.zeros((Ny,Nx, occupations))
        self.F        = F

        #create probability array
        self.omega = np.zeros(occupations)
        self.set_omega()

        #create velocity vectors
        self.c = np.zeros((2, occupations), dtype='int')
        self.set_c()     #initialize

        #constants that can be calculated outside of for-loops
        self.alpha = 1 - self.delta_t/tau
        self.beta  = self.delta_t/tau
        self.gamma = 3*(1-self.delta_t/(2*tau))

        self.initialize(1)

    def initialize(self, rho):
    	self.f[:,:,0] = 4*rho/9
    	self.f[:,:,1] = rho/9
    	self.f[:,:,2] = rho/9
    	self.f[:,:,3] = rho/9
    	self.f[:,:,4] = rho/9
    	self.f[:,:,5] = rho/36
    	self.f[:,:,6] = rho/36
    	self.f[:,:,7] = rho/36
    	self.f[:,:,8] = rho/36

    def initialize_other(self, ny, nx, occupations, value):
    	self.f[ny, nx, occupations] = value

    def set_omega(self):
    	self.omega[0] = 4/9
    	self.omega[1:5] = 1/9
    	self.omega[5:]  = 1/36

    def set_c(self):
    	self.c[1, 1] = 1
    	self.c[0, 2] = 1
    	self.c[1, 3] = -1
    	self.c[0, 4] = -1
    	self.c[:, 5] = 1
    	self.c[1, 6] = -1
    	self.c[0, 6] = 1
    	self.c[:, 7] = -1
    	self.c[1, 8] = 1
    	self.c[0, 8] = -1

    def f_eq(self):
    	u   = self.u
    	rho = self.rho

    	u_squared = np.multiply(u[:,:,0], u[:,:,0]) + np.multiply(u[:,:,1], u[:,:,1])
    	self.f_equil[:,:,0] = np.multiply(rho, 2*(2 - 3*u_squared)/9)
    	self.f_equil[:,:,1] = np.multiply(rho, (2 + 6*u[:,:,1] + 9*np.multiply(u[:,:,1], u[:,:,1]) - 3*u_squared)/18)
    	self.f_equil[:,:,2] = np.multiply(rho, (2 + 6*u[:,:,0] + 9*np.multiply(u[:,:,0], u[:,:,0]) - 3*u_squared)/18)
    	self.f_equil[:,:,3] = np.multiply(rho, (2 - 6*u[:,:,1] + 9*np.multiply(u[:,:,1], u[:,:,1]) - 3*u_squared)/18)
    	self.f_equil[:,:,4] = np.multiply(rho, (2 - 6*u[:,:,0] + 9*np.multiply(u[:,:,0], u[:,:,0]) - 3*u_squared)/18)
    	self.f_equil[:,:,5] = np.multiply(rho, (1 + 3*(u[:,:,1]+u[:,:,0]) + 9*np.multiply(u[:,:,0], u[:,:,1]) + 3*u_squared)/36)
    	self.f_equil[:,:,6] = np.multiply(rho, (1 - 3*(u[:,:,1]-u[:,:,0]) - 9*np.multiply(u[:,:,0], u[:,:,1]) + 3*u_squared)/36)
    	self.f_equil[:,:,7] = np.multiply(rho, (1 - 3*(u[:,:,1]+u[:,:,0]) + 9*np.multiply(u[:,:,0], u[:,:,1]) + 3*u_squared)/36)
    	self.f_equil[:,:,8] = np.multiply(rho, (1 + 3*(u[:,:,1]-u[:,:,0]) - 9*np.multiply(u[:,:,0], u[:,:,1]) + 3*u_squared)/36)

    def rho_and_u(self):
    	self.rho = np.sum(self.f, axis=2)
    	self.u[:,:,0] = np.sum(np.multiply(self.c[0,:], self.f[:,:,:]), axis=2)/self.rho + np.divide(self.F[0], (2*self.rho))
    	self.u[:,:,1] = np.sum(np.multiply(self.c[1,:], self.f[:,:,:]), axis=2)/self.rho + np.divide(self.F[1], (2*self.rho))


    def source_term(self):
    	F = self.F
    	c = self.c
    	u = self.u

    	self.S[:,:, 0] = self.gamma*self.omega[0]*( np.dot(F, c[:, :])[0]*(1+3*np.tensordot(c, u, axes=[0,2])[0,:,:]) - np.dot(u, F))
    	self.S[:,:, 1] = self.gamma*self.omega[1]*( np.dot(F, c[:, :])[1]*(1+3*np.tensordot(c, u, axes=[0,2])[1,:,:]) - np.dot(u, F))
    	self.S[:,:, 2] = self.gamma*self.omega[2]*( np.dot(F, c[:, :])[2]*(1+3*np.tensordot(c, u, axes=[0,2])[2,:,:]) - np.dot(u, F))
    	self.S[:,:, 3] = self.gamma*self.omega[3]*( np.dot(F, c[:, :])[3]*(1+3*np.tensordot(c, u, axes=[0,2])[3,:,:]) - np.dot(u, F))
    	self.S[:,:, 4] = self.gamma*self.omega[4]*( np.dot(F, c[:, :])[4]*(1+3*np.tensordot(c, u, axes=[0,2])[4,:,:]) - np.dot(u, F))
    	self.S[:,:, 5] = self.gamma*self.omega[5]*( np.dot(F, c[:, :])[5]*(1+3*np.tensordot(c, u, axes=[0,2])[5,:,:]) - np.dot(u, F))
    	self.S[:,:, 6] = self.gamma*self.omega[6]*( np.dot(F, c[:, :])[6]*(1+3*np.tensordot(c, u, axes=[0,2])[6,:,:]) - np.dot(u, F))
    	self.S[:,:, 7] = self.gamma*self.omega[7]*( np.dot(F, c[:, :])[7]*(1+3*np.tensordot(c, u, axes=[0,2])[7,:,:]) - np.dot(u, F))
    	self.S[:,:, 8] = self.gamma*self.omega[8]*( np.dot(F, c[:, :])[8]*(1+3*np.tensordot(c, u, axes=[0,2])[8,:,:]) - np.dot(u, F))

    def algorithm(self):
            #calculate

            self.source_term()
            self.f_eq()

            current_max_velocity = np.max(np.max(np.abs(self.u)))

            #save data
            self.velocity[:,:, :]  = self.u[:,:,:]
            self.density[:, :]     = self.rho[:,:]

            #stream
            self.f[:,:,:] = self.f[:,:,:]*self.alpha + self.f_equil*self.beta + self.S

            #propegate
            self.f_prev[:, :, :] = np.copy(self.f[:, :, :])
            self.f[:,:,1] = np.roll(self.f[:,:,1],  1, axis=1)
            self.f[:,:,3] = np.roll(self.f[:,:,3], -1, axis=1)
            self.f[:,:,2] = np.roll(self.f[:,:,2],  1, axis=0)
            self.f[:,:,4] = np.roll(self.f[:,:,4], -1, axis=0)
            self.f[:,:,5] = np.roll(self.f[:,:,5],  1, axis=0)
            self.f[:,:,5] = np.roll(self.f[:,:,5],  1, axis=1)
            self.f[:,:,6] = np.roll(self.f[:,:,6],  1, axis=0)
            self.f[:,:,6] = np.roll(self.f[:,:,6], -1, axis=1)
            self.f[:,:,7] = np.roll(self.f[:,:,7], -1, axis=0)
            self.f[:,:,7] = np.roll(self.f[:,:,7], -1, axis=1)
            self.f[:,:,8] = np.roll(self.f[:,:,8], -1, axis=0)
            self.f[:,:,8] = np.roll(self.f[:,:,8],  1, axis=1)

            self.f[0, :, 2]    = self.f_prev[0, :, 4]
            self.f[self.Ny-1, :, 4] = self.f_prev[self.Ny-1, :, 2]
            self.f[0, :, 5]    = self.f_prev[0, :, 7]
            self.f[self.Ny-1, :, 7] = self.f_prev[self.Ny-1, :, 5]
            self.f[0, :, 6]    = self.f_prev[0, :, 8]
            self.f[self.Ny-1, :, 8] = self.f_prev[self.Ny-1, :, 6]

            """
            if self.counter % 1000 == 0:
                np.save('../data/pouseill_150_50_'+str(self.counter), self.u[:,:,:])
                print(self.counter)
            """

            if abs(current_max_velocity - previous_max_velocity) < self.epsilon*np.max(np.max(self.u[:,:,1])):
                equil = True
                print(self.counter)
                #np.save('../data/velocity_dist_test' + str(int(self.Ny)), self.u[:,:,:])
                np.save('../data/pouseill_29_Ny_'+str(self.Ny), self.u[:,:,1])

            previous_max_velocity = current_max_velocity
            self.counter += 1

    def plot_final_velocity(self):
        x_axis = np.linspace(0, self.Nx-1, self.Nx)
        y_axis = np.linspace(0, self.Ny-1, self.Ny)

        plt.quiver(x_axis, y_axis, self.velocity[:,:,1], self.velocity[:,:,0])
        plt.show()
        #np.save('../data/test_pouseill_2', self.velocity[:,:,:])

NY = 150
NX = 50
f = np.array([0, 5*1e-8])
tau = 2
epsilon = 1e-9

instance = LBM(NX, NY, tau, f, epsilon)
instance.algorithm()

"""
NY = 150
NX = 50
f = np.array([0, 5*1e-8])
tau = 2
epsilon = 1e-9
height_list = np.arange(5, 195, 10, dtype='int')
print(height_list)
for heigts in height_list:
    print(heigts)
    instance = LBM(NX, heigts, tau, f, epsilon)
    instance.algorithm()
"""
