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
"""
print("added them together \n")
print("RHS: ", RHS)
expanded_RHS = expand(RHS)
print("\n \n Expanded RHS: ", expanded_RHS)
print("\n \n Simplifty RHS: ", simplify(RHS))
SIMP = simplify(expanded_RHS)
print("\n \n Simplifty expanded RHS: ", expand(SIMP))
print("\n \n Simplifty expanded RHS: ", latex(expand(SIMP)))
"""

my_RHS = kappa*xi*cosh(kappa*xi)*(gamma*gamma-9*kappa*kappa)/6 +sinh(kappa*xi)*(kappa*kappa-2*gamma*gamma/9) + sinh(kappa)*(3*kappa*kappa*kappa_p*xi*cosh(kappa_p*xi)/2-kappa_p*kappa_p*sinh(kappa_p*xi))/sinh(kappa_p)
my_RHS *= sin(2*kappa*eta)*P_1*kappa/(gamma*gamma)
difference = RHS - my_RHS

difference = difference.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
#print(difference)
#print(expand(difference))
print("\n\n\n simplified", simplify(expand(difference)))


my_RHS = kappa*xi*cosh(kappa*xi)*(gamma*gamma-9*kappa*kappa)/6 +sinh(kappa*xi)*(kappa*kappa-2*gamma*gamma/9) + sinh(kappa)*(3*kappa*kappa*kappa_p*xi*cosh(kappa_p*xi)/2-kappa_p*kappa_p*sinh(kappa_p*xi))/sinh(kappa_p)
fy = 2*sinh(kappa*xi)*kappa*kappa*gamma*gamma*(1-gamma*gamma/(3*kappa*kappa))/(3*(gamma*gamma+3*kappa*kappa)**2) + kappa*xi*cosh(kappa*xi)*(gamma*gamma*gamma*gamma/6 - 9*kappa*kappa*kappa*kappa/2 - gamma*gamma*kappa*kappa)/((gamma*gamma+3*kappa*kappa)**2) + sinh(kappa)*kappa_p*xi*cosh(kappa_p*xi)/(2*sinh(kappa_p))

print(simplify( (diff(fy, xi, xi) - (gamma*gamma+4*kappa*kappa)*fy + my_RHS).subs(kappa_p, sqrt(gamma*gamma+kappa*kappa))))