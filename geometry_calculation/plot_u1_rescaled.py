import numpy as np 
import matplotlib.pyplot as plt 
#import seaborn as sns
import pandas as pd

#plt.style.use("bmh")
#sns.color_palette("hls", 1)
import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

N = int(500)
epsilon = 0.3

Sc = 1
F0 = 1
pi = np.pi 
omega = 10
kappa = 1
Nt = 10
T = np.linspace(0, np.pi/omega, Nt)
eta = np.linspace(-2*np.pi, 2*np.pi/kappa, 300)
xi = np.linspace(-1, 1, 100)

u_x = np.zeros(len(xi), dtype="complex")
u_y = np.zeros(len(xi), dtype="complex")
flow_field = np.zeros((len(xi), len(eta), len(T), 2), dtype="complex")

a = 1
x = np.zeros((len(xi), len(eta)))
y = np.zeros((len(xi), len(eta)))

for i in range(len(eta)):
	y[:, i] = xi*a*(1+epsilon*np.sin(kappa*eta[i]))

for i in range(len(x)):
	x[i, :] = a*eta

X, Y = np.meshgrid(eta, xi)

for i in range(len(T)):
	t = T[i]
	psi1 = 0
	phi1 = 0

	for n in range(N):
		n = 2*n + 1
		alpha_n = n*n*pi*pi*Sc/4

		an = 4*Sc*F0*((-1)**(0.5*(n-1)))*(alpha_n-1j*omega)/(n*pi*(alpha_n*alpha_n + omega*omega))
		bn = -2*kappa*kappa*F0*((-1)**(0.5*(n-1)))/(kappa*kappa + n*n*pi*pi/4)

		dn = bn - an*n*pi*kappa*kappa/2 
		cn = an*n*n*pi*pi/2 + bn*n*pi/(kappa*kappa + n*n*pi*pi/4)

		theta_n = -2*kappa*F0*((-1)**(0.5*(n-1)))/(kappa*kappa + n*n*pi*pi/4)
		denominator = 1j*omega/Sc + kappa*kappa + n*n*pi*pi/4
		lambda_n = (1-n*n*pi*pi*0.5*(1/denominator + 1/(kappa*kappa+n*n*pi*pi/4)))

		psi1 += ((-1)**(0.5*(n-1)))*(theta_n*lambda_n*np.sinh(kappa) - dn*np.cosh(kappa))/denominator
		phi1 += ((-1)**(0.5*(n-1)))*(theta_n*lambda_n*np.cosh(kappa) + dn*np.sinh(kappa))/denominator

		u_x += cn*np.cos(n*pi*xi/2)/(1j*omega/Sc + kappa*kappa + n*n*pi*pi/4)
		u_x += dn*(xi*np.sin(n*pi*xi/2)/denominator + n*pi*np.cos(n*pi*xi/2)/(denominator*denominator))


		u_y += theta_n*np.sin(n*pi*xi/2)*(1-n*n*pi*pi*0.5/(kappa*kappa+n*n*pi*pi/4))/denominator
		u_y += theta_n*n*pi*0.5*(xi*np.cos(n*pi*xi/2)/denominator - n*pi*np.sin(n*pi*xi/2)/(denominator*denominator))
	
	psi1 *= 1j*omega/(Sc*kappa*np.cosh(2*kappa))*np.exp(1j*omega*t)
	phi1 *= 1j*omega/(Sc*kappa*np.cosh(2*kappa))*np.exp(1j*omega*t)

	u_x *= np.exp(1j*omega*t)
	u_y *= np.exp(1j*omega*t)

	for j in range(len(xi)):
		#flow_field[j, :, i, 0] = np.real(np.sin(kappa*eta)*u_x[j] + Sc*kappa*np.cosh(kappa*xi[j])*(psi1*np.sin(kappa*eta) - phi1*np.cos(kappa*eta))/(1j*omega))
		flow_field[j, :, i, 1] = np.real(np.cos(kappa*eta)*u_y[j] - Sc*kappa*np.sinh(kappa*xi[j])*(psi1*np.cos(kappa*eta) + phi1*np.sin(kappa*eta))/(1j*omega))

	u_x = np.zeros(len(xi), dtype="complex")
	u_y = np.zeros(len(xi), dtype="complex")
"""
for t in range(len(T)):
	X, Y = np.meshgrid(eta, xi)
	plt.quiver(X, Y, flow_field[:,:,t,0], flow_field[:,:,t,1])
	plt.draw()
	plt.pause(0.01)
	plt.clf()
plt.show()
"""
u_x = np.zeros(len(xi), dtype=complex)
u_x0 = np.zeros((len(xi), len(eta), len(T)), dtype=complex)

for i in range(len(T)):
	t = T[i]
	for n in range(N):
		n = 2*n +1
		alpha_n = n*n*pi*pi*Sc/4
		u_x += (-1)**(0.5*(n-1))*np.cos(n*np.pi*xi/2)*(alpha_n*np.cos(omega*t) + omega*np.sin(omega*t))/(n*(alpha_n**2 + omega**2))

	u_x *= 4/np.pi

	for j in range(len(xi)):
		u_x0[j, :, i] = u_x[j]

	u_x = np.zeros(len(xi))
"""
for t in range(len(T)):
	X, Y = np.meshgrid(eta, xi)
	Z = np.sqrt((epsilon*flow_field[:,:,t,0]+u_x0[:, :, t])*(epsilon*flow_field[:,:,t,0]+u_x0[:, :, t]) + epsilon*flow_field[:,:,t,1]*epsilon*flow_field[:,:,t,1])

	cp = plt.contourf(X, Y, Z, levels=np.linspace(0, 0.55, 30))
	cb = plt.colorbar(cp)

	plt.quiver(X, Y, epsilon*flow_field[:,:,t,0]+u_x0[:, :, t], epsilon*flow_field[:,:,t,1], linewidths=Z.flatten()/60)
	plt.draw()
	#plt.axis([min(eta), max(eta), -1.01, 1.1])
	plt.pause(0.01)
	plt.clf()
plt.show()
"""
max_vel = np.amax(np.sqrt(np.real( (epsilon*flow_field[:,:,:,0]+u_x0[:, :, :])*(epsilon*flow_field[:,:,:,0]+u_x0[:, :, :]) + flow_field[:,:,:,1]*flow_field[:,:,:,1])))

from scipy.interpolate import griddata

def streams(ax,xx,yy,u,v,speed,base_map=False):
	 x = np.linspace(np.real(xx).min(), xx.max(), 50)
	 y = np.linspace(np.real(yy).min(), yy.max(), 50)

	 xi, yi = np.meshgrid(x,y)

	 #then, interpolate your data onto this grid:

	 px = xx.flatten()
	 py = yy.flatten()
	 pu = u.flatten()
	 pv = v.flatten()
	 pspeed = speed.flatten()

	 gu = griddata(zip(px,py), pu, (xi,yi))
	 gv = griddata(zip(px,py), pv, (xi,yi))
	 gspeed = griddata(zip(px,py), pspeed, (xi,yi))

	 lw = np.real(2*gspeed/np.nanmax(gspeed))
	 #now, you can use x, y, gu, gv and gspeed in streamplot:

	 if base_map:
	    xx,yy = ax(xx,yy)
	    xi,yi = ax(xi,yi)

	 #ax.contour(xx,yy,speed, alpha=0.4)
	 #ax.plot(xx,yy,'-k',alpha=0.3)
	 #ax.plot(xx.T,yy.T,'-k',alpha=0.3)
	 #ax.plot(xi,yi,'-b',alpha=0.1)
	 #ax.plot(xi.T,yi.T,'-b',alpha=0.1)
	 ax.streamplot(x,y,np.real(gu),np.real(gv), density=2, color="k")

slicing_x = 10
slicing_y = 10
total_flow = np.zeros((len(xi), len(eta), len(T), 2))
for t in range(len(T)):
	for i in range(len(eta)):
		total_flow[:, i, t, 0] = np.real(epsilon*flow_field[:,i,t,0]+u_x0[:, i, t])
		total_flow[:, i, t, 1] = np.real(epsilon*flow_field[:,i,t,1])
"""
lin_x_axis = np.linspace(min(eta), max(eta), 1e3)
lin_y_axis = np.linspace(-1+epsilon, 1-epsilon, 1e3)
XX, YY = np.meshgrid(lin_x_axis, lin_y_axis)

import scipy.interpolate as sci 
fig, ax = plt.subplots()

points = np.zeros((len(xi)*len(eta), 2))
for i in range(len(points[:, 0])):
	points[i, 0] = x.flatten()[i]
	points[i, 1] = y.flatten()[i]

values = np.zeros((len(xi)*len(eta), len(T), 2))

for t in range(len(T)):
	for i in range(len(values[:, 0])-1):
		#print(np.concatenate(total_flow[:,: t+5, 0]).ravel())
		values[i, t, 0] = np.concatenate(total_flow[:,: t, 0]).ravel()[i]
		values[i, t, 1] = np.concatenate(total_flow[:,: t, 1]).ravel()[i]
"""
for i in range(len(T)):
	#u_x_inter = sci.griddata(points, total_flow[:,: t, 0], (XX, YY))
	#u_y_inter = sci.griddata(points, total_flow[:,:,t, 1], (XX, YY))
	#plt.streamplot(XX, YY, u_x_inter, u_y_inter, density=0.5, color='k', linewidth=lw)
	plt.clf()
	plt.title(r"$T=$"+str(T[i]*omega)[:4] + r"$[\omega^{-1}]$")
	Z = np.real(np.sqrt((epsilon*flow_field[:,:,i,0]+u_x0[:, :, i])*(epsilon*flow_field[:,:,i,0]+u_x0[:, :, i]) + epsilon*flow_field[:,:,i,1]*epsilon*flow_field[:,:,i,1]))
	#Z = np.real(np.sqrt((epsilon*flow_field[:,:,i,0])*(epsilon*flow_field[:,:,i,0]) + epsilon*flow_field[:,:,i,1]*epsilon*flow_field[:,:,i,1]))

	cp = plt.contourf(x, y, Z, levels=np.linspace(0, max_vel*1.001, 40))
	cb = plt.colorbar(cp)
	plt.quiver(x[1::slicing_x, 1::slicing_y], y[1::slicing_x, 1::slicing_y], np.real(epsilon*flow_field[1::slicing_x,1::slicing_y,i,0]+u_x0[1::slicing_x, 1::slicing_y, i]), np.real(epsilon*flow_field[1::slicing_x,1::slicing_y,i,1]))
	#streams(ax, x, y, np.real(epsilon*flow_field[:,:,i,0]+u_x0[:, :, i]), np.real(epsilon*flow_field[:,:,i,1]), Z)#, color='k')#, linewidth=lw, density=0.5,)

	plt.axis([min(eta), max(eta), -1 - epsilon*1.2, 1+epsilon*1.2])

	plt.draw()
	plt.pause(0.01)
	plt.xlabel(r"Horizontal position $\eta$ $[\kappa^{-1}]$", fontsize=12)
	plt.ylabel(r"Vertical position   $\xi$ ", fontsize=12)
	plt.savefig("figures/streamplot_scaled_epsilon1_nr"+str(i)+".pdf")
plt.show()
