from sympy import * 
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i, P_1, kappa_p, U2 = symbols("omega gamma kappa Sc F0 Pe rho i P_1 kappa_p U2")


order = 2
ux0  = F0*(1-cosh(gamma*xi)/cosh(gamma))/(gamma*gamma)
ux1 = ((P_1*kappa*cosh(kappa)/(gamma*gamma))*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + (F0*tanh(gamma)/gamma)*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = (kappa*P_1*sinh(kappa)/(gamma*gamma))*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
ux20   = F0*cosh(gamma*xi)/(4*cosh(gamma))*(1-xi*xi) + (P_1*sinh(kappa)/(2*gamma*gamma))*(kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa)+gamma*gamma*cosh(gamma*xi)/cosh(gamma)-kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p))


ux0 = simplify(series(ux0, gamma, n=order).removeO())
ux0 = simplify(series(ux0, kappa, n=order).removeO())
print("\n ux0 (1): ", ux0)

uy1 = uy1.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
uy1 = series(uy1, kappa_p, n=order).removeO()
uy1 = uy1.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
uy1 = series(uy1, kappa, n=order).removeO()
uy1 = simplify(series(uy1, gamma, n=order).removeO())
print("\n uy (1): ", uy1)

ux1 = ux1.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
ux1 = series(ux1, kappa_p, n=order).removeO()
ux1 = ux1.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
ux1 = series(ux1, kappa, n=order).removeO()
ux1 = simplify(series(ux1, gamma, n=order).removeO())
ux1 = simplify(series(ux1, kappa, n=order).removeO())
print("\n ux (1): ", ux1)


#ux20 = ux20.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
ux20 = series(ux20, kappa_p, n=order+1).removeO()
ux20 = ux20.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
ux20 = series(ux20, kappa, n=order+1).removeO()
ux20 = simplify(series(ux20, gamma, n=order+1).removeO())
ux20 = simplify(series(ux20, kappa, n=order+1).removeO())
print("\n ux (2): ", ux20)
print("boundary ux_0 (2): ", simplify(expand(simplify(ux20.subs(xi, 1)))))

"""
P1  = ( F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p)))
P1 = series(P1, kappa_p, n=order+1).removeO()
P1 = P1.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
P1 = series(P1, kappa, n=order+1).removeO()
P1 = simplify(series(P1, gamma, n=order+1).removeO())
P1 = simplify(series(P1, kappa, n=order+1).removeO())
print("\n P (1): ", P1)
print("boundary P (1): ", simplify(expand(simplify(P1.subs(xi, 1)))))
"""
#print(simplify(expand(simplify(ux20.subs(P_1, P1)))))

B_deriv = Pe*F0*(sinh(rho*xi)/sinh(rho)-xi)/(rho*rho)
B1 = F0*(1-xi*xi - 2/(rho*rho) + 2*cosh(rho*xi)/(rho*sinh(rho)))/(4*rho*rho) + Pe*F0*(1/(rho*rho)+xi*sinh(rho*xi)/(2*sinh(rho)) - cosh(rho*xi)*(1+rho/tanh(rho))/(2*rho*sinh(rho)))/(rho*rho)


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

#linear terms, not including veocity
RHS1 = kappa*kappa*xi*B_deriv/2 + (3+kappa*kappa*xi*xi)*diff(B_deriv, xi)/2 - kappa*kappa*xi*diff(B1, xi)/2 - diff(diff(B1, xi), xi)# + Pe*ux20 - Pe*U2

print("RHS1: ", simplify((simplify(RHS1))))
##non-linear terms
#RHS2 = (Pe/2)*(kappa*xi*ux0*)


my_RHS1 = xi*sinh(rho*xi)/sinh(rho)*F0*Pe/2*(  kappa*kappa/(2*rho*tanh(rho)) + kappa*kappa/(2*rho*rho) - 1/(2*rho*rho*Pe) - 1) + cosh(rho*xi)/sinh(rho)*F0*Pe/(2*rho)*(rho/tanh(rho) + 2 - 1/Pe) + xi*xi*F0*kappa*kappa*Pe/(rho*rho)*(rho/4-1+1/(4*Pe)) + F0/(2*rho*rho)  -3*Pe*F0/(2*rho*rho)
print("RHS1: ", simplify(expand(simplify(RHS1-my_RHS1))))
