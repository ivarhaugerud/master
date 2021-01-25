from sympy import * 

#define variables
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i = symbols("omega gamma kappa Sc F0 Pe rho i")

B0_deriv = Pe*F0*(sinh(rho*xi)/sinh(rho)-1)/(rho*rho)
B_plus = F0*(1-xi*xi - 2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(4*rho*rho) + Pe*F0*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)*(1+rho/tanh(rho))/(2*rho*sinh(rho)))/(rho*rho)
difference = diff(diff(B_plus, xi), xi) - rho*rho*B_plus
difference = (simplify(expand(difference)))
RHS = -F0*(1-xi*xi)/4 + Pe*F0*(rho*cosh(rho*xi)/sinh(rho) - 1)/(rho*rho)
#RHS = (diff(B0_deriv, xi)) - F0*(1-xi*xi)/4 
print("B_plus: ", simplify(expand(difference - RHS)))
print("B_plus boundary: ", simplify((diff(B_plus, xi).subs(xi, 1))))

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