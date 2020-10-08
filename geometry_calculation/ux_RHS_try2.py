from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
A, B, C, D, P_1, psi = symbols("A B C D P_1 psi")
i = symbols("i")

#define known quantities
u0  = A*(1-cosh(gamma*xi)/cosh(gamma))
ux1 = sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)*cos(kappa*eta)
P2  = (P_1*sin(2*kappa*eta)/3)*(kappa*xi*sinh(kappa*xi)/2 + cosh(kappa*xi)/3) + sin(2*kappa*eta)*psi*cosh(2*kappa*xi)

#define operators acting on known quantities
laplas2_u0 = xi*kappa*kappa*(3*cos(2*kappa*eta) + 1)*diff(u0, xi)/2 + (3*sin(kappa*eta)*sin(kappa*eta)+kappa*kappa*xi*xi*cos(kappa*eta)*cos(kappa*eta))*diff(u0, xi, xi)
laplas1_u1x = -2*kappa*xi*cos(kappa*eta)*diff(ux1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(ux1, xi) - 2*sin(kappa*eta)*diff(ux1, xi, xi)
grad1x_P1 = -kappa*xi*cos(kappa*eta)*diff(P1, xi)
laplas2_u0x = xi*kappa*kappa*(3*cos(kappa*eta)*cos(kappa*eta)-1)*diff(u0, xi) + (3*sin(kappa*eta)*sin(kappa*eta) + xi*xi*kappa*kappa*cos(kappa*eta)*cos(kappa*eta))*diff(u0, xi, xi)

#define all the terms
term1 = laplas2_u0x# - diff(u0, xi, xi) - diff(u0, eta, eta)
term2 = laplas1_u1x
term3 = -grad1x_P1
term4 = -diff(P2, eta)

#add together, the substitute back and simplify
RHS = term1 + term2 + term3 + term4
RHS = RHS.subs(A, 1/(gamma*gamma))
RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
#RHS = RHS.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
print(simplify(RHS))
RHS_no_eta = RHS.subs(eta, pi/(4*kappa))

RHS_no_eta   = simplify(RHS_no_eta)
RHS_with_eta = simplify((RHS - RHS_no_eta)/(cos(2*kappa*eta)))
print("\n\n", RHS_no_eta)
print("\n\n", RHS_with_eta.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma)))


my_RHS1 = cos(2*kappa*eta)*((3-kappa*kappa*xi*xi)*cosh(gamma*xi)/(2*cosh(gamma)) - 3*kappa*kappa*xi*sinh(gamma*xi)/(2*gamma*cosh(gamma)))
my_RHS2 = -cosh(gamma*xi)*(kappa*kappa*xi*xi+1)/(2*cosh(gamma)) - kappa*kappa*xi*sinh(gamma*xi)/(2*gamma*cosh(gamma))
print("Difference my RHS and theirs: ", simplify(RHS - my_RHS1-my_RHS2))



my_sol1 = cosh(gamma*xi)*(kappa*kappa-gamma*gamma/2-xi*xi*kappa**4)/(8*cosh(gamma)*kappa**4) - xi*sinh(gamma*xi)*(3*kappa*kappa+gamma*gamma)/(8*kappa*kappa*gamma*cosh(gamma))
print("Test my solutions:")
print("cos(2*kappa*eta)-term: ", simplify(diff(my_sol1, xi, xi) - (4*kappa*kappa+gamma*gamma)*my_sol1 + my_RHS1/(cos(2*kappa*eta))))

my_sol2 = cosh(gamma*xi)*((gamma**4)/2 - 2*gamma*gamma*kappa*kappa + 16*kappa**4 + kappa*kappa*xi*xi*(8*kappa**4 -4*gamma**2 *kappa**2  + (gamma**4)/4))/(cosh(gamma)*(gamma*gamma-4*kappa*kappa)**3) + kappa*gamma*kappa*xi*sinh(gamma*xi)*((8*kappa**4)/(gamma*gamma) + 4*kappa**2 - 3*gamma*gamma/2)/(cosh(gamma)*(gamma*gamma-4*kappa*kappa)**3)
print("\n\n non cos(2*kappa*eta)-term: ", simplify(diff(my_sol2, xi, xi) - (4*kappa*kappa)*my_sol2 + my_RHS2))