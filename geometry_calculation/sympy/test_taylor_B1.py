from sympy import * 

#define variables
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i = symbols("omega gamma kappa Sc F0 Pe rho i")


sol = F0*(1-xi*xi - 2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(4*rho*rho) + Pe*F0*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)*(1+rho/tanh(rho))/(2*rho*sinh(rho)))/(rho*rho)
difference = diff(diff(sol, xi), xi) - rho*rho*sol
difference = (simplify(expand(difference)))
RHS = -F0*(1-xi*xi)/4 + Pe*F0*(rho*cosh(rho*xi)/sinh(rho) - 1)/(rho*rho)
print(simplify(expand(difference - RHS)))
print(simplify((diff(sol, xi).subs(xi, 1))))


# SOLVING FOR B2
"""
A = (F0/2 - Pe*F0*(1+rho/tanh(rho))/2)/(rho*rho*rho*sinh(rho))
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
"""

# SOLVING FOR B0

RHS1 = Pe*kappa*F0*F0*(1-xi*xi)*(1-xi*xi-2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(16*rho*rho)
RHS2 = Pe*Pe*F0*F0*kappa*(1-xi*xi)*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)/(2*rho*sinh(rho))*(1+rho/tanh(rho)))/(4*rho*rho)
#RHS2 = Pe*Pe*F0*F0*kappa*(1-xi*xi)*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)))/(4*rho*rho)

#rint(simplify(expand(simplify(Pe*kappa*F0*sol*(1-xi*xi)/4)-RHS1-RHS2)))
#print(simplify(expand(diff(diff(sol, xi),xi)-RHS)))

#sol1 = Pe*kappa*F0*F0*(xi*xi/2 - xi*xi*xi*xi/6 + xi*xi*xi*xi*xi*xi/30 - xi*xi/(rho*rho) + xi*xi*xi*xi/(6*rho*rho) + 2*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) - 2*xi*xi*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) + 8*xi*sinh(rho*xi)/(rho*rho*rho*rho*sinh(rho))-12*cosh(rho*xi)/(rho*rho*rho*rho*rho*sinh(rho)))/(16*rho*rho)
#sol2 = Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(cosh(rho*xi)/(rho*sinh(rho))*(12/(rho*rho)+3*xi*xi-1-9*xi*tanh(rho*xi)/(rho*rho*rho)-rho*xi*xi*xi*tanh(rho*xi)/2+rho*xi*tanh(rho*xi)/2) +xi*xi*(1-xi*xi/6)/2 + (1+rho/tanh(rho))/(2*rho*rho*rho*sinh(rho))*(cosh(rho*xi)*(xi*xi-1+6/(rho*rho))-4*xi*sinh(rho*xi)/rho))
#sol2 = Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(12*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) - 9*xi*sinh(rho*xi)/(rho*rho*sinh(rho))+3*xi*xi*cosh(rho*xi)/(rho*sinh(rho))-cosh(rho*xi)/(rho*sinh(rho))-xi*xi*xi*xi/12-xi*xi*xi*sinh(rho*xi)/(2*sinh(rho))+xi*xi/2+xi*sinh(rho*xi)/(2*sinh(rho)) + (1+rho/tanh(rho))/(2*rho*sinh(rho))*(6*cosh(rho*xi)/(rho*rho) - 4*xi*sinh(rho*xi)/rho + xi*xi*cosh(rho*xi)-cosh(rho*xi)))
sol1 = Pe*F0*F0*kappa/(16*rho*rho)*((xi*xi/2)*(1-2/(rho*rho))+xi*xi*xi*xi*(1/(rho*rho)-1)/6 + xi*xi*xi*xi*xi*xi/30 + 2*cosh(rho*xi)/(rho*rho*rho*sinh(rho))*(1-xi*xi-6/(rho*rho)) + 8*xi*sinh(rho*xi)/(rho*rho*rho*rho*sinh(rho)))
print(simplify(expand(simplify(diff(diff(sol1, xi), xi)-(RHS1)))))


sol2 = Pe*Pe*F0*F0*kappa/(4*rho*rho*rho*rho)*(cosh(rho*xi)/(rho*sinh(rho))*(12/(rho*rho)+3*xi*xi-1) + xi*sinh(rho*xi)/(2*sinh(rho))*(1-xi*xi-18/(rho*rho))+xi*xi*(1-xi*xi/6)/2  - xi*sinh(rho*xi)/sinh(rho)*2*(1+rho/tanh(rho))/(rho*rho) + (1+rho/tanh(rho))*cosh(rho*xi)/(2*rho*sinh(rho))*(xi*xi-1+6/(rho*rho)) )
print(simplify(expand(simplify(diff(diff(sol2, xi), xi)-(RHS2)))))
