import numpy as np 
import matplotlib.pyplot as plt 

def coupled_finite_element_solver(N, n, x, alpha, couple_forward, couple_backward, f, Bc0, Bc1):
	#define general parameters
	#phi = np.zeros((N, len(x)))
	N_pos = np.linspace(min(x), max(x), N)
	Delta_x = N_pos[1]-N_pos[0]

	Nn = int(N*n)
	A   = np.zeros((Nn, Nn), dtype="complex")
	A_p = np.zeros((Nn, Nn), dtype="complex")
	b = np.zeros(Nn, dtype="complex")
	u = np.zeros((n, len(x)), dtype="complex")

	phi = np.zeros((N, len(x)), dtype="complex")
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
		A[(i+1)*N : (i+2)*N, i*N    :(i+1)*N] += couple_forward  *np.identity(N)*Delta_x
		A[i*N     : (i+1)*N, (i+1)*N:(i+2)*N] += couple_backward *np.identity(N)*Delta_x
	print(A)
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


#works for differential equation with constant terms, now just need coupeling to work as well
n = 2  #number of vectors 
N = 300 #length of each vector 
x = np.linspace(0, 1, int(1e5))

alpha = np.zeros(n) #self-interaction

alpha[0] =  1
alpha[1] =  3

#coupleing between vectors
couple_forward  =  2*np.ones(N)
couple_backward =  4*np.ones(N)

sol1 =  1j*np.exp(1j*x) - 1j*np.exp(-1j*x) + np.sqrt(5)*np.exp(np.sqrt(5)*x)/10 - np.sqrt(5)*np.exp(-np.sqrt(5)*x)/10
sol2 = -1j*np.exp(1j*x) + 1j*np.exp(-1j*x) + np.sqrt(5)*np.exp(np.sqrt(5)*x)/5  - np.sqrt(5)*np.exp(-np.sqrt(5)*x)/5

#since integral over last and first basis function is half the value of the others
couple_backward[0]  *= 0.5
couple_backward[-1] *= 0.5 

couple_forward[0]  *= 0.5
couple_forward[-1] *= 0.5 

#boundary conditions
Bc0 = np.zeros(n, dtype="complex")
Bc1 = np.zeros(n, dtype="complex")

Bc0[0] = np.gradient(sol1, x)[0]
Bc1[0] = np.gradient(sol1, x)[-1]

Bc0[1] = np.gradient(sol2, x)[0]
Bc1[1] = np.gradient(sol2, x)[-1]

print(Bc0, Bc1)
f = np.zeros((n, len(x)))
#f[0, :] =  (1+np.pi*np.pi)*np.sin(np.pi*x)
#f[1, :] =  (1+np.pi*np.pi)*np.sin(np.pi*x)


sol = coupled_finite_element_solver(N, n, x, alpha, couple_backward, couple_forward, f, Bc0, Bc1)

for i in range(len(sol[:,0])):
	plt.plot(x, np.real(sol[i,:]), label="num"+str(i))

plt.plot(x, np.real(sol1), "--")
plt.plot(x, np.real(sol2), "--")
plt.legend(loc="best")
plt.show()


for i in range(len(sol[:,0])):
	plt.plot(x, np.imag(sol[i,:]), label="num"+str(i))

plt.plot(x, np.imag(sol1), "--")
plt.plot(x, np.imag(sol2), "--")
plt.legend(loc="best")
plt.show()