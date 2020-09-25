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
laplas1_u1x = -2*kappa*xi*cos(kappa*eta)*diff(ux1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(ux1, xi) - 2*sin(kappa*eta)*diff(ux1, xi, xi)
grad1x_P1 = -kappa*xi*cos(kappa*eta)*diff(P1, xi)

#define all the terms
term1 = diff(u0, xi, xi) + diff(u0, eta, eta)
term2 = laplas1_u1x
term3 = grad1x_P1
term4 = diff(P2_partic, eta)

#add together, the substitute back and simplify
RHS = term1 + term2 + term3 + term4
RHS = RHS.subs(A, 1/(gamma*gamma))
RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
#RHS = RHS.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
RHS = RHS.subs(gamma, sqrt(i*omega/Sc))
RHS = RHS.subs(eta, pi/(2*kappa))

print("added them together \n")
print("RHS: ", RHS)
expanded_RHS = expand(RHS)
print("\n \n Expanded RHS: ", expanded_RHS)
print("\n \n Simplifty RHS: ", simplify(RHS))
SIMP = simplify(expanded_RHS)
print("\n \n Simplifty expanded RHS: ", expand(SIMP))
print("\n \n Simplifty expanded RHS: ", latex(expand(SIMP)))

 # + 2*kappa**2*xi*sinh(gamma*xi)/(gamma*cosh(gamma)) - 2*kappa_p**2*sin(eta*kappa)**2*cosh(kappa_p*xi)*tanh(gamma)/(gamma*cosh(kappa_p))

 #Simplifty expanded RHS:  P_1*kappa**2*xi*sin(eta*kappa)**2*sinh(kappa*xi)/3 - 2*P_1*kappa**2*xi*sinh(kappa*xi)/3 - 4*P_1*kappa*sin(eta*kappa)**2*cosh(kappa*xi)/9 + 2*P_1*kappa*cosh(kappa*xi)/9 + 3*P_1*kappa**4*xi*sin(eta*kappa)**2*sinh(kappa*xi)/gamma**2 - 2*P_1*kappa**4*xi*sinh(kappa*xi)/gamma**2 - 3*P_1*kappa**3*kappa_p*xi*sin(eta*kappa)**2*sinh(kappa_p*xi)*cosh(kappa)/(gamma**2*cosh(kappa_p)) + 2*P_1*kappa**3*kappa_p*xi*sinh(kappa_p*xi)*cosh(kappa)/(gamma**2*cosh(kappa_p)) - 2*P_1*kappa**3*sin(eta*kappa)**2*cosh(kappa*xi)/gamma**2 + 2*P_1*kappa*kappa_p**2*sin(eta*kappa)**2*cosh(kappa)*cosh(kappa_p*xi)/(gamma**2*cosh(kappa_p)) + 2*gamma*xi*sin(eta*kappa)**2*sinh(gamma*xi)/cosh(gamma) - 3*kappa**2*xi**2*sin(eta*kappa)**2*cosh(gamma*xi)/cosh(gamma) + 2*kappa**2*xi**2*cosh(gamma*xi)/cosh(gamma) + 4*sin(eta*kappa)**2*cosh(gamma*xi)/cosh(gamma) - cosh(gamma*xi)/cosh(gamma) + 3*kappa**2*kappa_p*xi*sin(eta*kappa)**2*sinh(kappa_p*xi)*tanh(gamma)/(gamma*cosh(kappa_p)) - 2*kappa**2*kappa_p*xi*sinh(kappa_p*xi)*tanh(gamma)/(gamma*cosh(kappa_p)) - 3*kappa**2*xi*sin(eta*kappa)**2*sinh(gamma*xi)/(gamma*cosh(gamma)) + 2*kappa**2*xi*sinh(gamma*xi)/(gamma*cosh(gamma)) - 2*kappa_p**2*sin(eta*kappa)**2*cosh(kappa_p*xi)*tanh(gamma)/(gamma*cosh(kappa_p))




