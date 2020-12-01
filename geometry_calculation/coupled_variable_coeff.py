import numpy as np 
import matplotlib.pyplot as plt 

def coupled_finite_element_solver(N, n, x, alpha, couple_forward, couple_backward, f, Bc0, Bc1):
	#define general parameters
	phi = np.zeros((N, len(x)))
	N_pos = np.linspace(min(x), max(x), N)
	Delta_x = N_pos[1]-N_pos[0]

	Nn = int(N*n)
	A   = np.zeros((Nn, Nn),  dtype="complex")
	A_p = np.zeros((Nn, Nn),  dtype="complex")
	b = np.zeros(Nn,          dtype="complex")
	u = np.zeros((n, len(x)), dtype="complex")

	phi = np.zeros((N, len(x)))

	#find form of phi's given our interval

	#find form of phi's given our interval
	for i in range(N):
		sta_index = np.argmin(abs(x-(N_pos[i]-Delta_x)))
		top_index = np.argmin(abs(x-(N_pos[i]        )))
		end_index = np.argmin(abs(x-(N_pos[i]+Delta_x)))
		phi[i, sta_index:top_index]   = np.linspace(0, 1, top_index-sta_index)
		phi[i, top_index:end_index]   = np.linspace(1, 0, end_index-top_index)
	
	phi[-1, -1] = 1 #some times the last element is not included, and is therefore set to zero if this is not done

	#calculate matrix elements using analytical results
	# A   = phi_i  * phi_j  matrix 
	# A_p = phi_i' * phi_j' matrix 

	for j in range(n):
		for i in range(N-1):
			A_p[j*N+i,     j*N+i] += 1
			A_p[j*N+i+1,   j*N+i] -= 1
			A_p[j*N+i,   j*N+i+1] -= 1
			A_p[j*N+i+1, j*N+i+1] += 1

			A[j*N+ i, j*N+ i]     += 2*alpha[j]
			A[j*N+ i+1, j*N+i]    += 1*alpha[j]
			A[j*N+ i, j*N+ i+1]   += 1*alpha[j]
			A[j*N+ i+1, j*N+i+1]  += 2*alpha[j]

	A_p *= 1/Delta_x
	A   *= Delta_x/6

	for i in range(n-1):
		A[(i+1)*N : (i+2)*N, i*N    :(i+1)*N] += couple_forward
		A[i*N     : (i+1)*N, (i+1)*N:(i+2)*N] += couple_backward

	#calculate source vector
	for j in range(n):
		for i in range(N):
			b[j*N+i] = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[i]+Delta_x/2)))] + f[j, np.argmin(abs(x-(N_pos[i]-Delta_x/2)))])/2

		b[j*N+0]   = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[0] + Delta_x/2)))])/2
		b[j*N+N-1] = -Delta_x*(f[j, np.argmin(abs(x-(N_pos[-1]- Delta_x/2)))])/2

		#if derivative is non-zero at boundary
		b[N*j+0]    -= Bc0[j]
		b[N*j+N-1]  += Bc1[j]

	sol = np.linalg.solve(A+A_p, b)

	#transfer back solution to regular basis
	for j in range(n):
		for i in range(N):
			u[j,:] += sol[j*N+i]*phi[i, :]
	return u

xi = np.linspace(-1, 1, int(1e5))
x = xi
n = 2 #number of vectors
N = 500
N_pos = np.linspace(-1, 1, N)
Delta = N_pos[1]-N_pos[0]

h = (1+4*np.pi*np.pi)*np.cos(np.pi*x)
f = -np.sin(2*np.pi*x)/2 - 1
source = (1+np.pi*np.pi)*(1+4*np.pi*np.pi)*np.cos(np.pi*x)/(1+np.sin(2*np.pi*x)/2) #np.gradient(np.gradient(h, xi), xi) - h )/f 
g = np.sin(np.pi*x) + h

#coupleing between vectors
couple_backward =  np.zeros((N, N))
couple_forward  =  np.zeros((N, N))

couple1 = np.sin(np.pi*xi)
for i in range(N-1):
	couple_forward[i, i]   = Delta*(couple1[np.argmin(abs(xi-(N_pos[i]-Delta/2)))] + 2*couple1[np.argmin(abs(xi-(N_pos[i])))] + couple1[np.argmin(abs(xi-(N_pos[i]+Delta/2)))])/6
	couple_forward[i, i+1] = Delta*couple1[np.argmin(abs(xi-(N_pos[i]+N_pos[i-1])/2))]/6
	couple_forward[i+1, i] = Delta*couple1[np.argmin(abs(xi-(N_pos[i+1]+N_pos[i])/2))]/6

couple_forward[0, 0]    = Delta*(couple1[np.argmin(abs(xi-(N_pos[0])))]  + couple1[np.argmin(abs(xi-(N_pos[0]  + Delta/2)))])/6
couple_forward[-1, -1]  = Delta*(couple1[np.argmin(abs(xi-(N_pos[-1])))] + couple1[np.argmin(abs(xi-(N_pos[-1] - Delta/2)))])/6

couple1 = source
for i in range(N-1):
	couple_backward[i, i]   = Delta*(couple1[np.argmin(abs(xi-(N_pos[i]-Delta/2)))] + 2*couple1[np.argmin(abs(xi-(N_pos[i])))] + couple1[np.argmin(abs(xi-(N_pos[i]+Delta/2)))])/6
	couple_backward[i, i+1] = Delta*couple1[np.argmin(abs(xi-(N_pos[i]+Delta/2)))]/6
	couple_backward[i+1, i] = Delta*couple1[np.argmin(abs(xi-(N_pos[i]+Delta/2)))]/6

couple_backward[0, 0]    = Delta*(couple1[np.argmin(abs(xi-(N_pos[0])))]  + couple1[np.argmin(abs(xi-(N_pos[0]  + Delta/2)))])/6
couple_backward[-1, -1]  = Delta*(couple1[np.argmin(abs(xi-(N_pos[-1])))] + couple1[np.argmin(abs(xi-(N_pos[-1] - Delta/2)))])/6


#boundary conditions
Bc0 = np.zeros(n) #BC at xi = -1
Bc1 = np.zeros(n) #BC at xi =  1

Bc0[0] = -np.pi
Bc1[0] = -np.pi

Bc0[1] = -np.pi
Bc1[1] = -np.pi

q = np.zeros((n, len(xi)))
q[0, :] = np.cos(np.pi*xi)*np.cos(np.pi*xi)
q[1, :] = -(1+np.pi*np.pi)*np.sin(np.pi*xi)

coeff = np.ones(n)


sol = coupled_finite_element_solver(N, n, xi, coeff, couple_backward, couple_forward, q, Bc0, Bc1)

plt.plot(xi, sol[0, :], label="f")
plt.plot(xi, sol[1, :], label="g")
plt.plot(xi, f, "--")
plt.plot(xi, g, "--")
plt.legend(loc="best")
plt.show()
