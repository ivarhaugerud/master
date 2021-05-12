from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
A, B, C, D, P_1 = symbols("A B C D P_1")
i = symbols("i")

#define known quantities
u0  = A*(1-cosh(gamma*xi)/cosh(gamma))
ux1 = sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)*cos(kappa*eta)
P2_partic = (P_1*sin(2*kappa*eta)/3)*(kappa*xi*sinh(kappa*xi)/2 + cosh(kappa*xi)/3)

#define operators acting on known quantities
laplas1_u1y = -2*kappa*xi*cos(kappa*eta)*diff(uy1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(uy1, xi) - 2*sin(kappa*eta)*diff(uy1, xi, xi)
grad1y_P1   = -sin(kappa*eta)*diff(P1, xi)


#g = sinh(kappa*xi)*2*gamma*gamma*(3*kappa*kappa-gamma*gamma)/(9*(gamma*gamma+3*kappa*kappa)**2) + kappa*xi*cosh(kappa*xi)*(gamma*gamma-9*kappa*kappa)/(6*(gamma*gamma + 3*kappa*kappa))
gg = kappa_p*sinh(kappa)*xi*cosh(kappa_p*xi)/(2*sinh(kappa_p))

RHS = (sinh(kappa)/(sinh(kappa_p)))*(3*kappa*kappa*kappa_p*xi*cosh(kappa_p*xi)/2 - kappa_p*kappa_p*sinh(kappa_p*xi))


print(simplify((diff(gg, xi, xi) - (gamma*gamma+4*kappa*kappa)*gg + RHS)).subs(kappa_p, sqrt(kappa*kappa+gamma*gamma)))

#-(-2*gamma**4*sinh(kappa*xi)/9 + 2*gamma**2*kappa**2*sinh(kappa*xi)/3 + gamma**2*kappa*xi*cosh(kappa*xi) + 3*kappa**3*xi*cosh(kappa*xi) - 2*kappa**2*sinh(kappa*xi))/(gamma**2 + 3*kappa**2)

#-(gamma**2 + 3*kappa**2)**2*(gamma**2*kappa*xi*cosh(kappa*xi)/6 - 2*gamma**2*sinh(kappa*xi)/9 - 3*kappa**3*xi*cosh(kappa*xi)/2 + kappa**2*sinh(kappa*xi))/gamma**2
#-gamma**2*kappa*xi*cosh(kappa*xi)/6 + 2*gamma**2*sinh(kappa*xi)/9 + kappa**3*xi*cosh(kappa*xi) - kappa**2*sinh(kappa*xi)/3 + 9*kappa**5*xi*cosh(kappa*xi)/(2*gamma**2) - 3*kappa**4*sinh(kappa*xi)/gamma**2

#-P_1*kappa*(9*gamma**2*kappa**2*sinh(kappa*xi) - 81*kappa**5*xi*cosh(kappa*xi) + 54*kappa**4*sinh(kappa*xi))/(18*gamma**2*(gamma**2 + 3*kappa**2))

"""
#define all the terms
term1 = laplas1_u1y
term2 = grad1y_P1
term3 = diff(P2_partic, xi)

#add together, the substitute back
RHS = term1 + term2 + term3
RHS = RHS.subs(A, 1/(gamma*gamma))
#RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
#RHS = 0

#print(simplify(RHS))
#P_1*kappa*(3*gamma**2*kappa*xi*cosh(kappa*xi) - 4*gamma**2*sinh(kappa*xi) - 27*kappa**3*xi*cosh(kappa*xi) + 27*kappa**2*kappa_p*xi*sinh(kappa)*cosh(kappa_p*xi)/sinh(kappa_p) + 18*kappa**2*sinh(kappa*xi) - 18*kappa_p**2*sinh(kappa)*sinh(kappa_p*xi)/sinh(kappa_p))*sin(2*eta*kappa)/(18*gamma**2)

#Define known LHS terms 
LHS_y_term1 =  P_1*kappa* ( kappa*xi*cosh(kappa*xi)*(gamma*gamma*gamma*gamma/(6*kappa*kappa*kappa*kappa) - 3*gamma*gamma/(2*kappa*kappa)) + sinh(kappa*xi)*(2*gamma*gamma/(3*kappa*kappa) - 2*gamma*gamma*gamma*gamma/(9*kappa*kappa*kappa*kappa)) )/(9*gamma*gamma)
LHS_y_term2 = (P_1*kappa*kappa_p*sinh(kappa)/(2*gamma*gamma*sinh(kappa_p)))*xi*cosh(kappa_p*xi)
#LHS_y_term3 = By*sinh(kappa_pp*xi) + sinh(2*kappa*xi)*

#Add together terms, and take derivatives
LHS         = LHS_y_term1+LHS_y_term2
full_LHS    = -diff(LHS, xi, xi) + sqrt(gamma)*LHS + 4*kappa*kappa*LHS


sol = full_LHS
#sol = sol.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
sol = sol.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
sol = sol.subs(i, sqrt(-1))
sol = sol.subs(gamma, sqrt(i*omega/Sc))
"""
"""
sol = sol.subs(omega, 3)
sol = sol.subs(kappa, 2)
sol = sol.subs(Sc, 1)
sol = sol.subs(xi, 17)
"""
"""
expanded = expand(sol)
print("\n \n Expanded full: ", expanded)
print("\n \n Simplifty full: ", simplify(expanded))
#SIMP = simplify(expanded)
#print("\n \n Simplifty expanded RHS: ", simplify(expand(SIMP)))
"""