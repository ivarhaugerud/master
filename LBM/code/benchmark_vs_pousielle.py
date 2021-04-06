import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci 
"""
f0 = 5*1e-8
u = np.loadtxt("../data/benchmark_testing_benchmark_pousielle_u.txt")

print(np.shape(u))

y = np.linspace(-1, 1, len(u[0, :]))
a = len(y)/2

U_num = sci.trapz(u[0, :], y)/2
U_ana = a*a*f0/(3*0.5)

print((U_num-U_ana)/U_ana)
print(U_num, U_ana)
plt.plot(y, u[0, :])
plt.plot(y, 3*U_num*(1-y*y)/2, "--")
plt.show()
"""
fig,ax = plt.subplots(figsize=(6,6))

ux = np.loadtxt("../data/benchmark_testing_stokes_ux.txt")
uy = np.loadtxt("../data/benchmark_testing_stokes_uy.txt")


y = np.linspace(0, len(ux[0, :]), len(ux[0, :]))
x = np.linspace(0, len(ux[:, 0]), len(ux[:, 0]))
X, Y = np.meshgrid(x, y)

ax.contourf(X, Y, np.transpose(np.sqrt(ux*ux+uy*uy)))
ax.axis("equal")
ax.streamplot(X, Y, np.transpose(ux), np.transpose(uy))

plt.show()