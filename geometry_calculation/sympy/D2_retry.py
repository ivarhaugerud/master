from sympy import * 
#import sympy as sy
"""
from numpy import *
import scipy.integrate as scp
import matplotlib.pyplot as plt 
F0 = 6
tau = 3.0
D = 0.1
Pe = 1/D
omega = 2*pi/tau
nu = 1e3
xi = linspace(-1, 1, int(1e4))

kappa = 0.01
omega = logspace(-1.5, 1.5, 60)
rho = sqrt(1j*omega/D)

one = F0**2*Pe**2*(0.125*rho**7*tanh(conjugate(rho))**3 - 0.416666666666667*rho**6*tanh(rho)**3*tanh(conjugate(rho))**3 + 0.5*rho**6*tanh(rho)**3*tanh(conjugate(rho)) + 0.5*rho**6*tanh(rho)*tanh(conjugate(rho))**3 - rho**5*(0.125*rho**2 + 0.375)*tanh(rho)**2*tanh(conjugate(rho))**3 - (0.03125*rho**4 - 0.0625*rho**2*conjugate(rho)**2 + 0.03125*conjugate(rho)**4)*tanh(rho)**3*conjugate(rho)**3 + 0.03125*(rho**4*conjugate(rho)**2 + 3.5*rho**4 - 2*rho**2*conjugate(rho)**4 - 8*rho**2*conjugate(rho)**2 + conjugate(rho)**6 + 0.5*conjugate(rho)**4)*tanh(rho)**3*tanh(conjugate(rho))**2*conjugate(rho))/(rho**10*tanh(rho)**3*tanh(conjugate(rho))**3)
two = F0**2*Pe**2*(0.125*rho**7*tanh(conjugate(rho))**3 - 0.416666666666667*rho**6*tanh(rho)**3*tanh(conjugate(rho))**3 + 0.5*rho**6*tanh(rho)**3*tanh(conjugate(rho)) + 0.5*rho**6*tanh(rho)*tanh(conjugate(rho))**3 - rho**5*(0.125*rho**2 + 0.375)*tanh(rho)**2*tanh(conjugate(rho))**3 - (0.03125*rho**4 - 0.0625*rho**2*conjugate(rho)**2 + 0.03125*conjugate(rho)**4)*tanh(rho)**3*conjugate(rho)**3 + 0.03125*(rho**4*conjugate(rho)**2 + 3.5*rho**4 - 2*rho**2*conjugate(rho)**4 - 8*rho**2*conjugate(rho)**2 + conjugate(rho)**6 + 0.5*conjugate(rho)**4)*tanh(rho)**3*tanh(conjugate(rho))**2*conjugate(rho))/(rho**10*tanh(rho)**3*tanh(conjugate(rho))**3)
two = F0**2*Pe**2*(0.125*rho**7*tanh(conjugate(rho))**3 - 0.416666666666667*rho**6*tanh(rho)**3*tanh(conjugate(rho))**3 + 0.5*rho**6*tanh(rho)**3*tanh(conjugate(rho)) + 0.5*rho**6*tanh(rho)*tanh(conjugate(rho))**3 - rho**5*(0.125*rho**2 + 0.375)*tanh(rho)**2*tanh(conjugate(rho))**3 - (rho**4/8)*tanh(rho)**3*conjugate(rho)**3 + 0.03125*(-rho**6 + 3.5*rho**4 - 2*rho**6 + 8*rho**4 - (rho)**6 + 0.5*(rho)**4)*tanh(rho)**3*tanh(conjugate(rho))**2*conjugate(rho))/(rho**10*tanh(rho)**3*tanh(conjugate(rho))**3)
two = rho/tanh(rho) - 4/(tanh(rho)*tanh(rho)) - rho/(tanh(rho)**3) + 3/(rho*tanh(rho))
two += conjugate(two) + 10/3 
two *= - F0*F0*Pe*Pe/(8*rho**4)

plt.plot(omega, one)
plt.plot(omega, two, "--")
plt.xscale("log")
plt.show()
"""
#define variables
xi = symbols("xi", real=True)
kappa, F0, Pe = symbols("kappa F0 Pe", real=True)
rho, gamma = symbols("rho gamma")

#print(simplify(F0**2*Pe**2*(0.125*rho**7*tanh(conjugate(rho))**3 - 0.416666666666667*rho**6*tanh(rho)**3*tanh(conjugate(rho))**3 + 0.5*rho**6*tanh(rho)**3*tanh(conjugate(rho)) + 0.5*rho**6*tanh(rho)*tanh(conjugate(rho))**3 - rho**5*(0.125*rho**2 + 0.375)*tanh(rho)**2*tanh(conjugate(rho))**3 - (rho**4/8)*tanh(rho)**3*conjugate(rho)**3 + 0.03125*(-rho**6 + 3.5*rho**4 - 2*rho**6 + 8*rho**4 - (rho)**6 + 0.5*(rho)**4)*tanh(rho)**3*tanh(conjugate(rho))**2*conjugate(rho))/(rho**10*tanh(rho)**3*tanh(conjugate(rho))**3)))


D0 = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))
#D0 = sy.cos(gamma)
a = series(D0, gamma, n=4)
print(simplify(a))
# = series(expr=D0, x=gamma, x0=0, n=2, dir="+")

"""
B0 =  F0*Pe*cosh(rho*xi)/(2*rho*rho*rho*sinh(rho)) - F0*Pe*xi*xi/(4*rho*rho) #
B0_deriv = diff(B0, xi)
B0_deriv_deriv = diff(B0_deriv, xi)

beta1 = 3/(2*rho**4) + (xi*xi - 1)/(4*rho*rho) + (xi*sinh(rho*xi) - cosh(rho*xi)/tanh(rho) - 2*cosh(rho*xi)/rho )/(2*rho*rho*sinh(rho))
beta1 *= Pe*F0 


beta2 = (1/4 - 2/(rho*rho) + 3/(rho**4) + xi*xi*(3/(rho*rho) - 1/2) + xi*xi*xi*xi/4) + (cosh(rho*xi)*(72/(rho*rho) + 6*xi*xi - 2)/rho + sinh(rho*xi)*(30*xi/(rho*rho) + xi**3 - xi) )/sinh(rho) + (2/rho + 1/tanh(rho))*(cosh(rho*xi)*(1-xi*xi - 10/(rho*rho)) - 4*xi*sinh(rho*xi)/rho)/sinh(rho)
beta2 *= Pe*Pe*F0*F0*kappa/(8*rho**4)
beta2 -= (F0*F0*Pe*Pe*kappa/(4*sqrt(2)*sinh(sqrt(2)*rho)*rho**7))*(40 + 8*rho/tanh(rho) - 3*rho*rho/(sinh(rho)**2))*cosh(sqrt(2)*rho*xi)

B2 = 5*xi*xi - 1 -3*xi*sinh(rho*xi)/sinh(rho) - 2*cosh(rho)*rho*xi*sinh(rho*xi)/(sinh(rho)*sinh(rho)) + rho*xi*xi*cosh(rho*xi)/sinh(rho) - xi*sinh(rho*xi)/sinh(rho)
B2 *= F0*Pe/(8*rho*rho)
B2 -= F0*Pe*cosh(rho*xi)*(rho*rho + 6 - 2*rho*rho/(tanh(rho)**2) -4*rho/tanh(rho) )/(8*rho*rho*rho*sinh(rho))
"""
"""
#full_sol = F0**2*Pe**2*(1.0*(xi*sinh(rho) - sinh(rho*xi))*(xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*sinh(rho)*sinh(conjugate(rho))*tanh(rho)**2*tanh(conjugate(rho))**2 - (xi*sinh(rho) - sinh(rho*xi))*(-(-10*xi*sinh(conjugate(rho))**2 + 2*(xi*cosh(xi*conjugate(rho))*conjugate(rho) + sinh(xi*conjugate(rho)))*cosh(conjugate(rho))*conjugate(rho) + (-xi**2*sinh(xi*conjugate(rho))*conjugate(rho)**2 + 2*xi*cosh(xi*conjugate(rho))*conjugate(rho) + 4*sinh(xi*conjugate(rho)))*sinh(conjugate(rho)))*tanh(conjugate(rho))**2 + (-(conjugate(rho)**2 + 6)*tanh(conjugate(rho))**2 + 4*tanh(conjugate(rho))*conjugate(rho) + 2*conjugate(rho)**2)*sinh(xi*conjugate(rho))*sinh(conjugate(rho)) - 2*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(conjugate(rho))*tanh(conjugate(rho)))*sinh(rho)*tanh(rho)**2 - (xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*((2*rho**2 + 4*rho*tanh(rho) - (rho**2 + 6)*tanh(rho)**2)*sinh(rho)*sinh(rho*xi) - 2*(-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*sinh(rho)*tanh(rho) - (2*rho*(rho*xi*cosh(rho*xi) + sinh(rho*xi))*cosh(rho) - 10*xi*sinh(rho)**2 + (-rho**2*xi**2*sinh(rho*xi) + 2*rho*xi*cosh(rho*xi) + 4*sinh(rho*xi))*sinh(rho))*tanh(rho)**2)*sinh(conjugate(rho))*tanh(conjugate(rho))**2 + (-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(rho)*sinh(conjugate(rho))*tanh(rho)*tanh(conjugate(rho)))/(8*rho**2*sinh(rho)**2*sinh(conjugate(rho))**2*tanh(rho)**2*tanh(conjugate(rho))**2*conjugate(rho)**2)
#integrated_full = integrate(full_sol, (xi, -1, 1), conds="none")
#print("result: ", integrated_full)
#print("\n\n\n simplified result: ", simplify(integrated_full))
"""
#rho = a + b*I
full = 1/2*F0**2*Pe**2*(7/12*rho**8*conjugate(rho)**2 + rho**8/2 + rho**8*conjugate(rho)**3/(2*tanh(conjugate(rho))) + 7/4*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 7/4*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1/2*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 1/2*rho**7*conjugate(rho)**4/tanh(rho) + 1/2*rho**7*conjugate(rho)**4/tanh(rho)**3 - 11/4*rho**6*conjugate(rho)**4 - rho**6*conjugate(rho)**2 - rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 9/4*rho**6*conjugate(rho)**4/tanh(rho)**2 + rho**5*conjugate(rho)**6/tanh(rho) - 1/4*rho**5*conjugate(rho)**4/tanh(rho) - rho**5*conjugate(rho)**6/tanh(rho)**3 + 11/4*rho**4*conjugate(rho)**6 + 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 1/4*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 9/4*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1/2*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4*rho**4*conjugate(rho)**6/tanh(rho)**2 - 1/2*rho**3*conjugate(rho)**8/tanh(rho) + 4*rho**3*conjugate(rho)**6/tanh(rho) + 1/2*rho**3*conjugate(rho)**8/tanh(rho)**3 - 7/12*rho**2*conjugate(rho)**8 + rho**2*conjugate(rho)**6 + 7/4*rho**2*conjugate(rho)**8/tanh(rho)**2 - 7/4*rho*conjugate(rho)**8/tanh(rho) - 1/2*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4)
#full -= -1/(kappa*kappa) + 1/(kappa*tanh(kappa)) + 1 #approx 8/3 + O(kappa^2)
full = simplify(full.subs(conjugate(rho)**2 , -rho*rho))
#print((full))
#print(simplify(full))
#maybe look at beta0 including kappa ^2




#print(full)
#print(simplify(full)) #
#print(series(full, r, n=4)) # -1.38777878078145e-17*F0**2*Pe**2/r**2 - 0.00687830687830688*F0**2*Pe**2 - 5.42101086242752e-20*I*F0**2*Pe**2*r + 1.60333493666827e-5*F0**2*Pe**2*r**2 + 4.2351647362715e-22*I*F0**2*Pe**2*r**3 + O(r**4
#print(limit(full, r, 0))
#print(limit(full, r, oo)) 

#print(simplify(full.subs(rho, 1000)))
#print("taylor small rho: ", simplify(series(full, rho, n=4)))
"""
from numpy import *
import scipy.integrate as scp
import matplotlib.pyplot as plt 
F0 = 6
tau = 3.0
D = 0.1
Pe = 1/D
omega = 2*pi/tau
nu = 1e3
xi = linspace(-1, 1, int(1e4))


kappa = 0.01
omega = logspace(-1.5, 1.5, 60)
D0 = zeros(len(omega), dtype="complex")
D_eff = zeros(len(omega), dtype="complex")
D_eff2 = zeros(len(omega), dtype="complex")

for i in range(len(omega)):
	rho = sqrt(1j*omega[i]/D)
	gamma = sqrt(1j*omega[i]/nu)
	D0[i] = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))
	D_eff[i] = 0.5*scp.trapz(F0**2*Pe**2*(1.0*(xi*sinh(rho) - sinh(rho*xi))*(xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*sinh(rho)*sinh(conjugate(rho))*tanh(rho)**2*tanh(conjugate(rho))**2 - (xi*sinh(rho) - sinh(rho*xi))*(-(-10*xi*sinh(conjugate(rho))**2 + 2*(xi*cosh(xi*conjugate(rho))*conjugate(rho) + sinh(xi*conjugate(rho)))*cosh(conjugate(rho))*conjugate(rho) + (-xi**2*sinh(xi*conjugate(rho))*conjugate(rho)**2 + 2*xi*cosh(xi*conjugate(rho))*conjugate(rho) + 4*sinh(xi*conjugate(rho)))*sinh(conjugate(rho)))*tanh(conjugate(rho))**2 + (-(conjugate(rho)**2 + 6)*tanh(conjugate(rho))**2 + 4*tanh(conjugate(rho))*conjugate(rho) + 2*conjugate(rho)**2)*sinh(xi*conjugate(rho))*sinh(conjugate(rho)) - 2*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(conjugate(rho))*tanh(conjugate(rho)))*sinh(rho)*tanh(rho)**2 - (xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*((2*rho**2 + 4*rho*tanh(rho) - (rho**2 + 6)*tanh(rho)**2)*sinh(rho)*sinh(rho*xi) - 2*(-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*sinh(rho)*tanh(rho) - (2*rho*(rho*xi*cosh(rho*xi) + sinh(rho*xi))*cosh(rho) - 10*xi*sinh(rho)**2 + (-rho**2*xi**2*sinh(rho*xi) + 2*rho*xi*cosh(rho*xi) + 4*sinh(rho*xi))*sinh(rho))*tanh(rho)**2)*sinh(conjugate(rho))*tanh(conjugate(rho))**2 + (-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(rho)*sinh(conjugate(rho))*tanh(rho)*tanh(conjugate(rho)))/(8*rho**2*sinh(rho)**2*sinh(conjugate(rho))**2*tanh(rho)**2*tanh(conjugate(rho))**2*conjugate(rho)**2), xi)#-F0**2*Pe**2*(4*rho**12*conjugate(rho)**2 + 12*rho**12 - 12*rho**12*conjugate(rho)/tanh(conjugate(rho)) + 3*rho**11*conjugate(rho)**4/tanh(rho) - 3*rho**11*conjugate(rho)**4/tanh(rho)**3 - 20*rho**10*conjugate(rho)**4 - 69*rho**10*conjugate(rho)**2 + 69*rho**10*conjugate(rho)**3/tanh(conjugate(rho)) - 15*rho**10*conjugate(rho)**4/sinh(rho)**2 - 12*rho**9*conjugate(rho)**6/tanh(rho) - 24*rho**9*conjugate(rho)**4/tanh(rho) + 12*rho**9*conjugate(rho)**6/tanh(rho)**3 + 40*rho**8*conjugate(rho)**6 + 165*rho**8*conjugate(rho)**4 - 123*rho**8*conjugate(rho)**5/tanh(conjugate(rho)) + 54*rho**8*conjugate(rho)**6/sinh(rho)**2 + 18*rho**7*conjugate(rho)**8/tanh(rho) + 57*rho**7*conjugate(rho)**6/tanh(rho) - 18*rho**7*conjugate(rho)**8/tanh(rho)**3 - 40*rho**6*conjugate(rho)**8 - 210*rho**6*conjugate(rho)**6 + 87*rho**6*conjugate(rho)**7/tanh(conjugate(rho)) - 72*rho**6*conjugate(rho)**8/sinh(rho)**2 - 12*rho**5*conjugate(rho)**10/tanh(rho) - 39*rho**5*conjugate(rho)**8/tanh(rho) + 12*rho**5*conjugate(rho)**10/tanh(rho)**3 + 20*rho**4*conjugate(rho)**10 + 150*rho**4*conjugate(rho)**8 - 21*rho**4*conjugate(rho)**9/tanh(conjugate(rho)) + 42*rho**4*conjugate(rho)**10/sinh(rho)**2 + 3*rho**3*conjugate(rho)**12/tanh(rho) + 3*rho**3*conjugate(rho)**10/tanh(rho) - 3*rho**3*conjugate(rho)**12/tanh(rho)**3 - 4*rho**2*conjugate(rho)**12 - 57*rho**2*conjugate(rho)**10 - 9*rho**2*conjugate(rho)**12/sinh(rho)**2 + 3*rho*conjugate(rho)**12/tanh(rho) + 9*conjugate(rho)**12)/(6*rho**4*(rho**4 - 2*rho**2*conjugate(rho)**2 + conjugate(rho)**4)*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4) - F0**2*Pe**2*(-4*rho**12*conjugate(rho)**2 + 9*rho**12 + 3*rho**12*conjugate(rho)**3/tanh(conjugate(rho)) + 3*rho**12*conjugate(rho)/tanh(conjugate(rho)) - 3*rho**12*conjugate(rho)**3/tanh(conjugate(rho))**3 - 9*rho**12*conjugate(rho)**2/sinh(conjugate(rho))**2 + 20*rho**10*conjugate(rho)**4 - 57*rho**10*conjugate(rho)**2 - 12*rho**10*conjugate(rho)**5/tanh(conjugate(rho)) + 3*rho**10*conjugate(rho)**3/tanh(conjugate(rho)) + 12*rho**10*conjugate(rho)**5/tanh(conjugate(rho))**3 + 42*rho**10*conjugate(rho)**4/sinh(conjugate(rho))**2 - 21*rho**9*conjugate(rho)**4/tanh(rho) - 40*rho**8*conjugate(rho)**6 + 150*rho**8*conjugate(rho)**4 + 18*rho**8*conjugate(rho)**7/tanh(conjugate(rho)) - 39*rho**8*conjugate(rho)**5/tanh(conjugate(rho)) - 18*rho**8*conjugate(rho)**7/tanh(conjugate(rho))**3 - 72*rho**8*conjugate(rho)**6/sinh(conjugate(rho))**2 + 87*rho**7*conjugate(rho)**6/tanh(rho) + 40*rho**6*conjugate(rho)**8 - 210*rho**6*conjugate(rho)**6 - 12*rho**6*conjugate(rho)**9/tanh(conjugate(rho)) + 57*rho**6*conjugate(rho)**7/tanh(conjugate(rho)) + 12*rho**6*conjugate(rho)**9/tanh(conjugate(rho))**3 + 54*rho**6*conjugate(rho)**8/sinh(conjugate(rho))**2 - 123*rho**5*conjugate(rho)**8/tanh(rho) - 20*rho**4*conjugate(rho)**10 + 165*rho**4*conjugate(rho)**8 + 3*rho**4*conjugate(rho)**11/tanh(conjugate(rho)) - 24*rho**4*conjugate(rho)**9/tanh(conjugate(rho)) - 3*rho**4*conjugate(rho)**11/tanh(conjugate(rho))**3 - 15*rho**4*conjugate(rho)**10/sinh(conjugate(rho))**2 + 69*rho**3*conjugate(rho)**10/tanh(rho) + 4*rho**2*conjugate(rho)**12 - 69*rho**2*conjugate(rho)**10 - 12*rho*conjugate(rho)**12/tanh(rho) + 12*conjugate(rho)**12)/(6*rho**4*(rho**4 - 2*rho**2*conjugate(rho)**2 + conjugate(rho)**4)*(-rho**6 + 3*rho**4*conjugate(rho)**2 - 3*rho**2*conjugate(rho)**4 + conjugate(rho)**6)*conjugate(rho)**4) + F0**2*Pe**2*(rho**2*(-0.25*rho**4*tanh(rho)*conjugate(rho) + 0.25*rho*tanh(conjugate(rho))*conjugate(rho)**4 + (0.0833333333333333*rho**4*conjugate(rho)**2 + 0.25*rho**4 - 0.0833333333333333*rho**2*conjugate(rho)**4 - 0.25*conjugate(rho)**4)*tanh(rho)*tanh(conjugate(rho)))*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*tanh(rho)*tanh(conjugate(rho))*conjugate(rho)**2 + (rho**2 - conjugate(rho)**2)*(0.25*rho**6*(rho**4 - conjugate(rho)**4)*tanh(rho)**2*conjugate(rho)**2 - rho**6*(-rho**4 + 1.5*rho**2*conjugate(rho)**2 + 1.5*conjugate(rho)**4)*tanh(rho)**2*tanh(conjugate(rho))*conjugate(rho) + 3.33066907387548e-16*rho**5*(rho**2 - conjugate(rho)**2)*tanh(rho)*tanh(conjugate(rho))*conjugate(rho)**5 + 0.25*rho**2*(rho**4 - conjugate(rho)**4)*tanh(conjugate(rho))**2*conjugate(rho)**6 + rho*(1.5*rho**4 + 1.5*rho**2*conjugate(rho)**2 - 1.0*conjugate(rho)**4)*tanh(rho)*tanh(conjugate(rho))**2*conjugate(rho)**6 + (0.0333333333333333*rho**10*conjugate(rho)**4 - 0.5*rho**10*conjugate(rho)**2 - 1.25*rho**10 - 0.1*rho**8*conjugate(rho)**6 + 0.5*rho**8*conjugate(rho)**4 + 1.5*rho**8*conjugate(rho)**2 + 0.1*rho**6*conjugate(rho)**8 + 1.75*rho**6*conjugate(rho)**4 - 0.0333333333333333*rho**4*conjugate(rho)**10 - 0.5*rho**4*conjugate(rho)**8 - 1.75*rho**4*conjugate(rho)**6 + 0.5*rho**2*conjugate(rho)**10 - 1.5*rho**2*conjugate(rho)**8 + 1.25*conjugate(rho)**10)*tanh(rho)**2*tanh(conjugate(rho))**2))/(rho**6*(rho**2 - conjugate(rho)**2)*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*tanh(rho)**2*tanh(conjugate(rho))**2*conjugate(rho)**6)
	D_eff2[i] = 0.5*F0**2*Pe**2*(0.583333333333333*rho**8*conjugate(rho)**2 + 0.5*rho**8 + 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho)) + 1.75*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 1.75*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 0.5*rho**7*conjugate(rho)**4/tanh(rho) + 0.5*rho**7*conjugate(rho)**4/tanh(rho)**3 - 2.75*rho**6*conjugate(rho)**4 - 1.0*rho**6*conjugate(rho)**2 - 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4.0*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4.0*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 2.25*rho**6*conjugate(rho)**4/tanh(rho)**2 + 1.0*rho**5*conjugate(rho)**6/tanh(rho) - 0.25*rho**5*conjugate(rho)**4/tanh(rho) - 1.0*rho**5*conjugate(rho)**6/tanh(rho)**3 + 2.75*rho**4*conjugate(rho)**6 + 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 0.25*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 2.25*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4.0*rho**4*conjugate(rho)**6/tanh(rho)**2 - 0.5*rho**3*conjugate(rho)**8/tanh(rho) + 4.0*rho**3*conjugate(rho)**6/tanh(rho) + 0.5*rho**3*conjugate(rho)**8/tanh(rho)**3 - 0.583333333333333*rho**2*conjugate(rho)**8 + 1.0*rho**2*conjugate(rho)**6 + 1.75*rho**2*conjugate(rho)**8/tanh(rho)**2 - 1.75*rho*conjugate(rho)**8/tanh(rho) - 0.5*conjugate(rho)**8)/(rho**4*(1.0*rho**6 - 3.0*rho**4*conjugate(rho)**2 + 3.0*rho**2*conjugate(rho)**4 - 1.0*conjugate(rho)**6)*conjugate(rho)**4)
	#D_eff2[i] -= -1/(kappa*kappa) + 1/(kappa*tanh(kappa)) + 1 #the contributions from the BCs

plt.xscale("log")
#plt.yscale("log")
plt.plot(omega, D_eff, label="without BC")
plt.plot(omega, D_eff2, "--", label="with BC")
plt.xlabel(r"frequency $\omega$")
plt.ylabel(r"Second order effective diffusion")
plt.legend(loc="best")
plt.savefig("../../analytic_D2.pdf")
plt.show()

eps = 0.2
plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r"frequency $\omega$")
plt.ylabel(r"Effective diffusion")
plt.plot(omega, D0 + eps*eps*D_eff2, label="second order")
plt.plot(omega, D0, "--",label="zeroth order")
plt.legend(loc="best")
plt.savefig("../../analytic_D.pdf")
plt.show()





import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sci
from scipy.signal import savgol_filter
import os 

plt.style.use(['science','no-latex', 'grid'])
import matplotlib
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


tau = np.logspace(-2, 2, 10)
kappa = np.array([0.2, 0.6, 1.0, 1.4, 1.8, 2.2])
epsilon = 0.3
dt  = 0.004
nu  = 8000
Dm   = 0.7
F0  = 20000/nu 
Sc = nu 
Pe = 1/Dm


D  = np.zeros((len(tau), len(kappa)))

T = int(750/2)
B = np.load("../data_test/vary_omega_taylor.npy")
print(np.shape(B))

D_0 = np.zeros(len(tau))
D_eff2 = np.zeros((len(tau), len(kappa)))
from numpy import *

for i in range(len(tau)):
	omega = 2*np.pi/tau[i]
	gamma = np.sqrt(1j*omega/Sc)
	gamma_c = np.conj(gamma)
	rho = np.sqrt(1j*omega/Dm)
	rho_c = np.conj(rho)
	print(omega, gamma, rho)
	D_0[i] = 1 + Pe*Pe*F0*F0*np.tanh(gamma)*np.tanh(gamma_c)/(4*gamma*gamma_c*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/np.tanh(gamma) - gamma_c/np.tanh(gamma_c)) - 1/(rho*rho)*(rho/np.tanh(rho) - rho_c/np.tanh(rho_c)))
	for j in range(len(kappa)):
		k = kappa[j]
		#print(sci.trapz(B[i, j, -T:, 8], B[i, j, -T:, 0])/tau)
		print(np.shape(B[i, j, :, 8]))
		if (abs(D[i, j] - sci.trapz(np.trim_zeros(B[i, j, -2*T:-T, 8]), np.trim_zeros(B[i, j, -2*T:-T, 0]))/(tau[i]/2))/D[i,j]) < 0.02:
			D[i, j] = sci.trapz(np.trim_zeros(B[i, j, -T:, 8]), np.trim_zeros(B[i, j, -T:, 0]))/(tau[i]/2)
		#plt.plot(np.trim_zeros(B[i, j, :, 0]), np.trim_zeros(B[i, j, :, 8]), label="%3.2f" % kappa[j])
		#plt.plot(np.trim_zeros(B[i, j, -T:, 0]), np.trim_zeros(B[i, j, -T:, 8]))
		#plt.show()
		D_eff2[i, j] = 0.5*F0**2*Pe**2*(0.583333333333333*rho**8*conjugate(rho)**2 + 0.5*rho**8 + 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho)) + 1.75*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 1.75*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 0.5*rho**8*conjugate(rho)**3/tanh(conjugate(rho))**3 - 0.5*rho**7*conjugate(rho)**4/tanh(rho) + 0.5*rho**7*conjugate(rho)**4/tanh(rho)**3 - 2.75*rho**6*conjugate(rho)**4 - 1.0*rho**6*conjugate(rho)**2 - 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho)) - 4.0*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) + 4.0*rho**6*conjugate(rho)**4/tanh(conjugate(rho))**2 + 1.0*rho**6*conjugate(rho)**5/tanh(conjugate(rho))**3 + 2.25*rho**6*conjugate(rho)**4/tanh(rho)**2 + 1.0*rho**5*conjugate(rho)**6/tanh(rho) - 0.25*rho**5*conjugate(rho)**4/tanh(rho) - 1.0*rho**5*conjugate(rho)**6/tanh(rho)**3 + 2.75*rho**4*conjugate(rho)**6 + 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho)) + 0.25*rho**4*conjugate(rho)**5/tanh(conjugate(rho)) - 2.25*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 0.5*rho**4*conjugate(rho)**7/tanh(conjugate(rho))**3 - 4.0*rho**4*conjugate(rho)**6/tanh(rho)**2 - 0.5*rho**3*conjugate(rho)**8/tanh(rho) + 4.0*rho**3*conjugate(rho)**6/tanh(rho) + 0.5*rho**3*conjugate(rho)**8/tanh(rho)**3 - 0.583333333333333*rho**2*conjugate(rho)**8 + 1.0*rho**2*conjugate(rho)**6 + 1.75*rho**2*conjugate(rho)**8/tanh(rho)**2 - 1.75*rho*conjugate(rho)**8/tanh(rho) - 0.5*conjugate(rho)**8)/(rho**4*(1.0*rho**6 - 3.0*rho**4*conjugate(rho)**2 + 3.0*rho**2*conjugate(rho)**4 - 1.0*conjugate(rho)**6)*conjugate(rho)**4)
		#D_eff2[i, j] -= -1/(k*k) + 1/(k*tanh(k)) + 1 #the contributions from the BCs

plt.plot(2*np.pi/tau, D_eff2[:, 0])
plt.xscale("log")
plt.show()
for i in range(len(kappa[:2])):
	plt.plot(2*np.pi/tau, D[:, i])
	plt.plot(2*np.pi/tau, epsilon*epsilon*D_eff2[:, i]+D_0, "--")
plt.plot(2*np.pi/tau, D_0, "k")
plt.xscale("log")
plt.show()

"""

#-F0**2*Pe**2*(0.125*rho**3/tanh(rho) - 0.125*rho**3/tanh(rho)**3 + 0.416666666666667*rho**2 + 0.125*rho**2*conjugate(rho)/tanh(conjugate(rho)) - 0.5*rho**2/tanh(conjugate(rho))**2 - 0.5*rho**2/tanh(rho)**2 + 0.375*rho/tanh(rho) - 0.375*conjugate(rho)/tanh(conjugate(rho)) + 0.125*conjugate(rho)**3/tanh(conjugate(rho))**3)/rho**6
