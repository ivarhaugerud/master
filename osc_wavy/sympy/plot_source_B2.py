from numpy import *
import numpy as np  
import matplotlib.pyplot as plt 
import os
import scipy.integrate as sci

plt.style.use(['science','no-latex'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'

nu = 3.6
omega = (2*np.pi)/3
F0 = 12/nu
Sc = nu
Pe = 6
kappa = 0.4# np.arange(0.1, 2.5, 0.4)
xi = np.linspace(-1, 1, int(1e4))

rho = np.sqrt(1j*omega)
gamma = np.sqrt(1j*omega/Sc)
kappa_p = np.sqrt(kappa*kappa+gamma*gamma)
P_1  = (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p)))

B0_deriv       = Pe*F0*(sinh(rho*xi)/sinh(rho)-xi)/(2*rho*rho)
B0_deriv_deriv = Pe*F0*(rho*cosh(rho*xi)/sinh(rho)-1)/(2*rho*rho)
B1 = F0*(1-xi*xi - 2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(4*rho*rho) + Pe*F0*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)*(1+rho/tanh(rho))/(2*rho*sinh(rho)))/(rho*rho)

linear = (0.5*kappa*kappa*xi*B0_deriv + (3+kappa*kappa*xi*xi)*B0_deriv_deriv/2 - kappa*kappa*xi*np.gradient(B1, xi)/2 - np.gradient(np.gradient(B1, xi), xi))/2 + Pe*(F0*gamma*gamma*(2*xi*xi-xi*xi*xi*xi-1)/8 + F0*(1-xi*xi)/4 + P_1*kappa*(1-xi*xi)*(2-gamma*gamma)/4 + F0*gamma*gamma/15 - P_1*kappa*(1-gamma*gamma/2)/3)/2
#plt.plot(xi, linear)

sol1 = xi*xi*cosh(rho*xi)*F0/(sinh(rho)*16*rho*rho*rho)*(kappa*kappa+Pe*(2*rho*rho-kappa*kappa)-Pe*rho*kappa*kappa/tanh(rho)) + xi*sinh(rho*xi)/sinh(rho)*F0/(4*rho*rho)*(1-Pe/2-Pe*rho/tanh(rho) - (kappa*kappa+Pe*(2*rho*rho-kappa*kappa)-Pe*rho*kappa*kappa/(tanh(rho)))/(4*rho*rho)) - F0*kappa*kappa*(2*Pe-1)*(xi*xi+2/(rho*rho))/(4*rho*rho*rho*rho) - F0*(3*Pe/2-1)/(2*rho*rho*rho*rho)
sol2   = Pe*F0*gamma*gamma*(24/(rho*rho*rho*rho) + 12*xi*xi/(rho*rho) - 4/(rho*rho) + xi*xi*xi*xi -2*xi*xi + 1)/(16*rho*rho) + Pe*(2/(rho*rho) + xi*xi - 1)*(F0 + P_1*kappa*(2-gamma*gamma))/(8*rho*rho) - Pe*(F0*gamma*gamma/15 - P_1*kappa*(1-gamma*gamma/2)/3)/(2*rho*rho)
sol = sol1+sol2
sol += -cosh(rho*xi)*(-F0*Pe*kappa**2*rho**2 + 2*F0*Pe*rho**4*tanh(rho)**2 - 4*F0*Pe*rho**4 + F0*kappa**2*rho*tanh(rho) - F0*rho**3*(Pe*kappa**2 + 4*Pe - 4)*tanh(rho) + F0*(24*Pe*gamma**2 - 15*Pe*kappa**2 + 7*kappa**2)*tanh(rho)**2 + rho**2*(-F0*Pe*kappa**2 + F0*kappa**2 + 4*F0 - 4*P_1*Pe*gamma**2*kappa + 8*P_1*Pe*kappa)*tanh(rho)**2)/(sinh(rho)*16*rho**5*tanh(rho)**2)

#kappa = np.linspace(0.01, 5, int(1e3))
rho = np.sqrt(1j*np.linspace(0.01, 100, int(1e3)))
D_para = F0*(15*F0*Pe**2*kappa**2*rho**4 - 30*F0*Pe**2*rho**6*tanh(rho)**4 + 120*F0*Pe**2*rho**6*tanh(rho)**2 - 90*F0*Pe**2*rho**6 + 15*F0*Pe*kappa**2*rho**3*(4*Pe - 3)*tanh(rho) + 10*F0*Pe*rho**5*(Pe*kappa**2 + 15*Pe - 12)*tanh(rho)**3 - 10*F0*Pe*rho**5*(Pe*kappa**2 + 15*Pe - 12)*tanh(rho) + 15*F0*rho**2*(24*Pe**2*gamma**2 + 16*Pe**2*kappa**2 - 7*Pe*kappa**2 + 2*kappa**2)*tanh(rho)**2 + 15*F0*rho*(120*Pe**2*gamma**2 - 13*Pe**2*kappa**2 + 2*Pe*kappa**2 + 10*kappa**2)*tanh(rho)**3 - 60*F0*(36*Pe**2*gamma**2 + 2*Pe**2*kappa**2 - 2*Pe*kappa**2 + 3*kappa**2)*tanh(rho)**4 - rho**4*(15*F0*Pe**2*kappa**2 + 150*F0*Pe**2 - 10*F0*Pe*kappa**2 - 240*F0*Pe + 30*F0 + 60*P_1*Pe**2*gamma**2*kappa - 120*P_1*Pe**2*kappa)*tanh(rho)**2 + rho**4*(16*F0*Pe**2*gamma**2 + 130*F0*Pe**2 - 2*F0*Pe*kappa**2 - 280*F0*Pe + 8*F0*kappa**2 + 50*F0 + 100*P_1*Pe**2*gamma**2*kappa - 200*P_1*Pe**2*kappa)*tanh(rho)**4 - rho**3*(105*F0*Pe**2*kappa**2 + 90*F0*Pe**2 - 45*F0*Pe*kappa**2 - 360*F0*Pe + 90*F0 + 180*P_1*Pe**2*gamma**2*kappa - 360*P_1*Pe**2*kappa)*tanh(rho)**3 - rho**2*(840*F0*Pe**2*gamma**2 + 45*F0*Pe**2*kappa**2 - 480*F0*Pe**2 - 50*F0*Pe*kappa**2 + 720*F0*Pe + 70*F0*kappa**2 - 120*F0 - 240*P_1*Pe**2*gamma**2*kappa + 480*P_1*Pe**2*kappa)*tanh(rho)**4)/(480*rho**8*tanh(rho)**4)


plt.plot(rho, D_para)
#plt.xscale("log")
plt.show()

plt.plot(xi, sol)
plt.plot(xi, sol*F0*(1 - xi**2)/4)
plt.show()
"""

ux0 = F0*(1 - xi**2)/2
uy1 = F0*kappa*xi*(1 - xi**2)/2
ux1 = F0*(xi**2 - 1)/2

A = F0*(1-Pe-rho*Pe/tanh(rho))/(2*rho*rho*rho*sinh(rho))
B2 = (Pe*F0*F0*kappa/(32*rho*rho*rho*rho))*(xi*xi-1+3/(rho*rho)) + (Pe*F0*F0*kappa/(16*rho*rho*rho*rho))*(1/(2*rho*rho) - xi*xi*xi*xi/2+xi*xi/2-4/(rho*rho*rho*rho)-4*xi*xi/(rho*rho))
B2 += A*Pe*F0*kappa/(4*rho*rho)*(10*cosh(rho*xi)/(rho*rho) + 4*xi*sinh(rho*xi)/rho + cosh(rho*xi)*(xi*xi-1))
B2 +=  Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(1/(2*rho*rho*rho*rho) + 36*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) + xi*xi/(2*rho*rho) + 15*xi*sinh(rho*xi)/(rho*rho*sinh(rho)) + 3*xi*xi*cosh(rho*xi)/(rho*sinh(rho)) + xi*xi*xi*sinh(rho*xi)/(2*sinh(rho)) )
B2 -= Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(1/(2*rho*rho) + cosh(rho*xi)/(rho*sinh(rho)) + xi*sinh(rho*xi)/(2*sinh(rho)) )
B2 += -(F0**2*Pe*kappa*(-3*Pe*rho**2/sinh(rho)**2 + 11*Pe*rho/tanh(rho) + 45*Pe + 3*rho/tanh(rho) + 5)/(4*rho**6))*cosh(sqrt(2)*rho*xi)/(sqrt(2)*rho*sinh(sqrt(2)*rho))



B0  = Pe*F0*F0*kappa/(16*rho*rho)*((xi*xi/2)*(1-2/(rho*rho))+xi*xi*xi*xi*(1/(rho*rho)-1)/6 + xi*xi*xi*xi*xi*xi/30 + 2*cosh(rho*xi)/(rho*rho*rho*sinh(rho))*(1-xi*xi-6/(rho*rho)) + 8*xi*sinh(rho*xi)/(rho*rho*rho*rho*sinh(rho)))
B0 += Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(cosh(rho*xi)/(rho*sinh(rho))*(15/(rho*rho)+7*xi*xi/2-3/2) + xi*sinh(rho*xi)/(2*sinh(rho))*(1-xi*xi-22/(rho*rho)-4/(rho*tanh(rho)))+xi*xi*(1-xi*xi/6)/2  + cosh(rho*xi)/(2*tanh(rho)*sinh(rho))*(xi*xi-1+6/(rho*rho)) )
A = -kappa*(5*F0**2*Pe**2*rho**2/12 - F0**2*Pe**2*rho**2/tanh(rho)**2/4 - 3/4*F0**2*Pe**2*rho/tanh(rho) + F0**2*Pe**2 + F0**2*Pe*rho**4/30 - F0**2*Pe*rho**2/12 + F0**2*Pe*rho/tanh(rho)/4 - F0**2*Pe/4 + rho**6)/rho**6
B0 += A*cosh(kappa*xi)/(kappa*sinh(kappa))

B_minus_1 = (B0 + B2)/4
plt.plot(xi, B_minus_1)
plt.show()

non_linear =  0.5*kappa*xi*ux0*np.gradient(B_minus_1) - 0.5*uy1*np.gradient(B_minus_1) + 0.5*kappa*ux1*B_minus_1 
plt.plot(xi, non_linear)
plt.show()
"""