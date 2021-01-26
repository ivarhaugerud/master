from sympy import * 

"""
from numpy import *
import matplotlib.pyplot as plt 

Pe = 6
F0 = 3
kappa = 0.5
rho = linspace(0.01, 1e2, int(1e3))
answ = F0**2*(-30*Pe**2*kappa**2*rho**4 + 30*Pe**2*kappa**2*rho**4/tanh(rho)**4 - 120*Pe**2*kappa**2*rho**3/tanh(rho) + 120*Pe**2*kappa**2*rho**3/tanh(rho)**3 - 755*Pe**2*kappa**2*rho**2 + 615*Pe**2*kappa**2*rho**2/tanh(rho)**2 + 915*Pe**2*kappa**2*rho/tanh(rho) - 1680*Pe**2*kappa**2 - 10*Pe**2*rho**6 + 40*Pe**2*rho**6/tanh(rho)**2 - 30*Pe**2*rho**6/tanh(rho)**4 + 120*Pe**2*rho**5/tanh(rho) - 120*Pe**2*rho**5/tanh(rho)**3 + 365*Pe**2*rho**4 - 285*Pe**2*rho**4/tanh(rho)**2 - 525*Pe**2*rho**3/tanh(rho) + 960*Pe**2*rho**2 + 16*Pe*kappa**2*rho**4 + 60*Pe*kappa**2*rho**3/tanh(rho) - 60*Pe*kappa**2*rho**3/tanh(rho)**3 + 150*Pe*kappa**2*rho**2 - 270*Pe*kappa**2*rho**2/tanh(rho)**2 + 90*Pe*kappa**2*rho/tanh(rho) + 240*Pe*kappa**2 - 60*Pe*rho**5/tanh(rho) + 60*Pe*rho**5/tanh(rho)**3 - 290*Pe*rho**4 + 210*Pe*rho**4/tanh(rho)**2 + 450*Pe*rho**3/tanh(rho) - 720*Pe*rho**2 + 8*kappa**2*rho**4 - 70*kappa**2*rho**2 + 30*kappa**2*rho**2/tanh(rho)**2 + 150*kappa**2*rho/tanh(rho) - 180*kappa**2 + 50*rho**4 - 30*rho**4/tanh(rho)**2 - 90*rho**3/tanh(rho) + 120*rho**2)/(240*rho**8)

plt.yscale("log")
plt.xscale("log")
plt.plot(rho, answ)
plt.show()
"""
#define variables
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i, P_1 = symbols("omega gamma kappa Sc F0 Pe rho i P_1")

B0_deriv = Pe*F0*(sinh(rho*xi)/sinh(rho)-xi)/(2*rho*rho)
B_plus = F0*(1-xi*xi - 2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(4*rho*rho) + Pe*F0*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)*(1+rho/tanh(rho))/(2*rho*sinh(rho)))/(rho*rho)
difference = diff(diff(B_plus, xi), xi) - rho*rho*B_plus
difference = (simplify(expand(difference)))
RHS = 2*(diff(B0_deriv, xi)) - F0*(1-xi*xi)/4 
print("B_plus: ", simplify(expand(difference - RHS)))
print("B_plus boundary: ", simplify((diff(B_plus, xi).subs(xi, 1))))

sol1 = xi*xi*cosh(rho*xi)*F0/(sinh(rho)*16*rho*rho*rho)*(kappa*kappa+Pe*(2*rho*rho-kappa*kappa)-Pe*rho*kappa*kappa/tanh(rho)) + xi*sinh(rho*xi)/sinh(rho)*F0/(4*rho*rho)*(1-Pe/2-Pe*rho/tanh(rho) - (kappa*kappa+Pe*(2*rho*rho-kappa*kappa)-Pe*rho*kappa*kappa/(tanh(rho)))/(4*rho*rho)) - F0*kappa*kappa*(2*Pe-1)*(xi*xi+2/(rho*rho))/(4*rho*rho*rho*rho) - F0*(3*Pe/2-1)/(2*rho*rho*rho*rho)
sol2   = Pe*F0*gamma*gamma*(24/(rho*rho*rho*rho) + 12*xi*xi/(rho*rho) - 4/(rho*rho) + xi*xi*xi*xi -2*xi*xi + 1)/(16*rho*rho) + Pe*(2/(rho*rho) + xi*xi - 1)*(F0 + P_1*kappa*(2-gamma*gamma))/(8*rho*rho) - Pe*(F0*gamma*gamma/15 - P_1*kappa*(1-gamma*gamma/2)/3)/(2*rho*rho)
B2 = sol1+sol2
B2 += -cosh(rho*xi)*(-F0*Pe*kappa**2*rho**2 + 2*F0*Pe*rho**4*tanh(rho)**2 - 4*F0*Pe*rho**4 + F0*kappa**2*rho*tanh(rho) - F0*rho**3*(Pe*kappa**2 + 4*Pe - 4)*tanh(rho) + F0*(24*Pe*gamma**2 - 15*Pe*kappa**2 + 7*kappa**2)*tanh(rho)**2 + rho**2*(-F0*Pe*kappa**2 + F0*kappa**2 + 4*F0 - 4*P_1*Pe*gamma**2*kappa + 8*P_1*Pe*kappa)*tanh(rho)**2)/(sinh(rho)*16*rho**5*tanh(rho)**2)


D_para = B0_deriv*B0_deriv*(1+kappa*kappa*xi)/2 + (B_plus*B_plus*kappa*kappa/2 + diff(B_plus, xi)*diff(B_plus, xi)/2) - B0_deriv*(diff(B_plus, xi)+kappa*kappa*xi*B_plus -2*diff(B2, xi))
D_para = simplify(expand(simplify(integrate(D_para, (xi, -1, 1)))))/2 #4 since integrating over sin^2(kn) and xi
print(D_para)
"""

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

D_para += kappa*B0 + kappa*B2 + kappa*xi*diff(B0+B2, xi) + (diff(B0, xi)*diff(B0, xi) + diff(B2, xi)*diff(B2, xi))/2

D_para = simplify(expand(simplify(D_para)))
print(D_para, "\n\n")
integal = simplify(expandsimplify((integrate(D_para, (xi, -1, 1)))))
print(integal, "\n\n")
print(simplify(series(integal, rho)))
"""
#D_para = 0.5*diff(B)

# SOLVING FOR B2
"""
A = F0*(1-Pe-rho*Pe/tanh(rho))/(2*rho*rho*rho*sinh(rho))
RHS = (Pe*F0*kappa*(1-xi*xi)/(4))*(F0*(1-xi*xi-2/(rho*rho))/(4*rho*rho))# 
RHS += Pe*F0*kappa*A*(1-xi*xi)*cosh(rho*xi)/4 
RHS += (Pe*F0*kappa*(1-xi*xi)/(4))*(Pe*F0*(1/(rho*rho) + xi*sinh(rho*xi)/(2*sinh(rho)))/(rho*rho))


sol = (Pe*F0*F0*kappa/(32*rho*rho*rho*rho))*(xi*xi-1+3/(rho*rho)) + (Pe*F0*F0*kappa/(16*rho*rho*rho*rho))*(1/(2*rho*rho) - xi*xi*xi*xi/2+xi*xi/2-4/(rho*rho*rho*rho)-4*xi*xi/(rho*rho))
sol += A*Pe*F0*kappa/(4*rho*rho)*(10*cosh(rho*xi)/(rho*rho) + 4*xi*sinh(rho*xi)/rho + cosh(rho*xi)*(xi*xi-1))
sol +=  Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(1/(2*rho*rho*rho*rho) + 36*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) + xi*xi/(2*rho*rho) + 15*xi*sinh(rho*xi)/(rho*rho*sinh(rho)) + 3*xi*xi*cosh(rho*xi)/(rho*sinh(rho)) + xi*xi*xi*sinh(rho*xi)/(2*sinh(rho)) )
sol -= Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(1/(2*rho*rho) + cosh(rho*xi)/(rho*sinh(rho)) + xi*sinh(rho*xi)/(2*sinh(rho)) )
sol += -(F0**2*Pe*kappa*(-3*Pe*rho**2/sinh(rho)**2 + 11*Pe*rho/tanh(rho) + 45*Pe + 3*rho/tanh(rho) + 5)/(4*rho**6))*cosh(sqrt(2)*rho*xi)/(sqrt(2)*rho*sinh(sqrt(2)*rho))
difference = simplify(expand(diff(diff(sol, xi), xi) - 2*rho*rho*sol - RHS))
print("B2", difference)

# SOLVING FOR B0
RHS1 = Pe*kappa*F0*F0*(1-xi*xi)*(1-xi*xi-2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(16*rho*rho)
RHS2 = Pe*Pe*F0*F0*kappa*(1-xi*xi)*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)/(2*rho*sinh(rho))*(1+rho/tanh(rho)))/(4*rho*rho)
sol1 = Pe*F0*F0*kappa/(16*rho*rho)*((xi*xi/2)*(1-2/(rho*rho))+xi*xi*xi*xi*(1/(rho*rho)-1)/6 + xi*xi*xi*xi*xi*xi/30 + 2*cosh(rho*xi)/(rho*rho*rho*sinh(rho))*(1-xi*xi-6/(rho*rho)) + 8*xi*sinh(rho*xi)/(rho*rho*rho*rho*sinh(rho)))
sol2 = Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(cosh(rho*xi)/(rho*sinh(rho))*(15/(rho*rho)+7*xi*xi/2-3/2) + xi*sinh(rho*xi)/(2*sinh(rho))*(1-xi*xi-22/(rho*rho)-4/(rho*tanh(rho)))+xi*xi*(1-xi*xi/6)/2  + cosh(rho*xi)/(2*tanh(rho)*sinh(rho))*(xi*xi-1+6/(rho*rho)) )

print("B0 partic: ", simplify(expand(simplify(diff(diff(sol2+sol1, xi), xi)-(RHS2+RHS1)))))
A = -kappa*(5*F0**2*Pe**2*rho**2/12 - F0**2*Pe**2*rho**2/tanh(rho)**2/4 - 3/4*F0**2*Pe**2*rho/tanh(rho) + F0**2*Pe**2 + F0**2*Pe*rho**4/30 - F0**2*Pe*rho**2/12 + F0**2*Pe*rho/tanh(rho)/4 - F0**2*Pe/4 + rho**6)/rho**6
sol3 = A*cosh(kappa*xi)/(kappa*sinh(kappa))
print("B0 homo: ", simplify(expand(simplify(diff(diff(sol3, xi), xi)-kappa*kappa*sol3))))

sol2 += sol3
print("BOUNDARY CONDITIONS B0 at xi= 1: ", simplify(expand(simplify((diff(sol1+sol2, xi).subs(xi,  1))))))
print("BOUNDARY CONDITIONS B0 at xi=-1: ", simplify(expand(simplify((diff(sol1+sol2, xi).subs(xi, -1))))))
"""