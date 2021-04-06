from sympy import * 
xi, eta = symbols("xi eta")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
A, B, C, D, P_1 = symbols("A B C D P_1")
i = symbols("i")

u0  = A*(1-cosh(gamma*xi)/cosh(gamma))
ux1 = sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)*cos(kappa*eta)


grad1_u1 = -xi*kappa*cos(kappa*eta)*diff(ux1, xi) - sin(kappa*eta)*diff(uy1, xi)
grad2_u0 = sin(kappa*eta)*xi*kappa*cos(kappa*eta)*diff(u0, xi)
laplas2_u0x = xi*kappa*kappa*(3*cos(kappa*eta)*cos(kappa*eta)-1)*diff(u0, xi) + (3*sin(kappa*eta)*sin(kappa*eta) + xi*xi*kappa*kappa*cos(kappa*eta)*cos(kappa*eta))*diff(u0, xi, xi)
laplas1_u1x = -2*kappa*xi*cos(kappa*eta)*diff(ux1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(ux1, xi) - 2*sin(kappa*eta)*diff(ux1, xi, xi)
laplas1_u1y = -2*kappa*xi*cos(kappa*eta)*diff(uy1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(uy1, xi) - 2*sin(kappa*eta)*diff(uy1, xi, xi)
grad1x_P1 = -kappa*xi*cos(kappa*eta)*diff(P1, xi)
grad1y_P1 = -sin(kappa*eta)*diff(P1, xi)

term1 = diff(grad1_u1+grad2_u0, xi, xi) + diff(grad1_u1+grad2_u0, eta, eta)
term2 = -(grad1_u1 + grad2_u0)*(i*omega/Sc)
term3 = -diff(laplas2_u0x, eta)
term4 = -diff(laplas1_u1x, eta) - diff(laplas1_u1y, xi)
term5 = -diff(grad1x_P1, eta) - diff(grad1y_P1, xi)

RHS = term1 + term2 + term3 + term4 + term5
RHS = RHS.subs(A, 1/(gamma*gamma))
RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
RHS = RHS.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
RHS = RHS.subs(gamma, sqrt(i*omega/Sc))

RHS = simplify(RHS)
print("\n \n Simplified RHS: ", RHS)


my_RHS = -P_1*kappa*kappa*kappa*xi*sinh(kappa*xi)*sin(2*kappa*eta)/2
difference = RHS-my_RHS
difference = difference.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference = difference.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
difference = difference.subs(gamma, sqrt(i*omega/Sc))
print("\n\n Difference my RHS and proper RHS: ", simplify(difference))


g = sin(2*kappa*eta)*(P_1*kappa*xi*sinh(kappa*xi)/6 + P_1*cosh(kappa*xi)/9)
my_sol_diff = diff(g, xi, xi) - 4*kappa*kappa*g - RHS

my_sol_diff = my_sol_diff.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
my_sol_diff = my_sol_diff.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
my_sol_diff = my_sol_diff.subs(gamma, sqrt(i*omega/Sc))
print("\n\nTest my solution g(xi): ", simplify(my_sol_diff))