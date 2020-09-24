from sympy import * 
"""
xi, eta, t = symbols("xi eta t")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
A, B, C, D, P_1 = symbols("A B C D P_1")
i = symbols("i")

u0  = exp(i*omega*t)*(1-cosh(gamma*xi)/cosh(gamma))
ux1 = exp(i*omega*t)*sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = exp(i*omega*t)*D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*exp(i*omega*t)*cosh(kappa*xi)*cos(kappa*eta)


grad1_u1 = - xi*kappa*cos(kappa*eta)*diff(ux1, xi) - sin(kappa*eta)*diff(ux1, xi)
grad2_u0 = sin(kappa*eta)*xi*kappa*cos(kappa*eta)*diff(u0, xi)
laplas2_u0x = xi*kappa*kappa*(3*cos(kappa*eta)*cos(kappa*eta)-1)*diff(u0, xi) + (3*sin(kappa*eta)*sin(kappa*eta) + xi*xi*kappa*kappa*cos(kappa*eta)*cos(kappa*eta))*diff(u0, xi, xi)
laplas1_u1x = -2*kappa*xi*cos(kappa*eta)*diff(ux1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(ux1, xi) - 2*sin(kappa*eta)*diff(ux1, xi, xi)
laplas1_u1y = -2*kappa*xi*cos(kappa*eta)*diff(uy1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(uy1, xi) - 2*sin(kappa*eta)*diff(uy1, xi, xi)
grad1x_P1 = -kappa*xi*cos(kappa*eta)*diff(P1, xi)
grad1y_P1 = -sin(kappa*eta)*diff(P1, xi)

term1 = diff(grad1_u1+grad2_u0, xi, xi) + diff(grad1_u1+grad2_u0, eta, eta)
term2 = -diff(grad1_u1 + grad2_u0, t)/Sc
term3 = -diff(laplas2_u0x, eta)
term4 = -diff(laplas1_u1x, eta) - diff(laplas1_u1y, xi)
term5 = -diff(grad1x_P1, eta) - diff(grad1y_P1, xi)

RHS = term1 + term2 + term3 + term4 + term5
print("added them together")
print(simplify(RHS))

"""

### Do same calculation without time variable

xi, eta = symbols("xi eta")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
A, B, C, D, P_1 = symbols("A B C D P_1")
i = symbols("i")

u0  = (1-cosh(gamma*xi)/cosh(gamma))
ux1 = sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)*cos(kappa*eta)


grad1_u1 = - xi*kappa*cos(kappa*eta)*diff(ux1, xi) - sin(kappa*eta)*diff(ux1, xi)
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
print("added them together \n")
print("RHS: ", RHS)
expanded_RHS = expand(RHS)
print("\n \n Expanded RHS: ", expanded_RHS)
print("\n \n Simplifty RHS: ", simplify(RHS))
print("\n \n Simplifty expanded RHS: ", simplify(expanded_RHS))
