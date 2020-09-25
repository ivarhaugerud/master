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
grad1y_P1 = -sin(kappa*eta)*diff(P1, xi)

#define all the terms
term1 = laplas1_u1y
term2 = grad1y_P1
term3 = diff(P2_partic, xi)

#add together, the substitute back and simplify
RHS = term1 + term2 + term3
RHS = RHS.subs(A, 1/(gamma*gamma))
RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
#RHS = RHS.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
#RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
#RHS = RHS.subs(gamma, sqrt(i*omega/Sc))

print("added them together \n")
print("RHS: ", RHS)
expanded_RHS = expand(RHS)
print("\n \n Expanded RHS: ", expanded_RHS)
print("\n \n Simplifty RHS: ", simplify(RHS))
SIMP = simplify(expanded_RHS)
print("\n \n Simplifty expanded RHS: ", expand(SIMP))
print("\n \n Simplifty expanded RHS: ", latex(expand(SIMP)))




