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
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
print("added them together \n")
print("RHS: ", latex(RHS))
expanded_RHS = expand(RHS)
print("\n \n Expanded RHS: ", expanded_RHS)
print("\n \n Simplifty RHS: ", simplify(RHS))
SIMP = simplify(expanded_RHS)
print("\n \n Simplifty expanded RHS: ", simplify(SIMP))
print("\n \n Simplifty expanded RHS: ", latex(simplify(SIMP)))


"""
Simplified expanded RHS
(-A*Sc*gamma**3*kappa*xi*sinh(gamma*xi)/cosh(gamma) - 2*A*Sc*gamma**2*kappa**3*xi**2*cosh(gamma*xi)/cosh(gamma) + 4*A*Sc*gamma**2*kappa*cosh(gamma*xi)/cosh(gamma) - 2*A*Sc*gamma*kappa**3*xi*sinh(gamma*xi)/cosh(gamma) + A*gamma*i*kappa*omega*xi*sinh(gamma*xi)/cosh(gamma) + B*Sc*gamma**2*kappa*xi*sqrt(gamma**2 + kappa**2)*sinh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) - 2*B*Sc*gamma**2*kappa*cosh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) - 3*B*Sc*kappa**4*xi*sinh(kappa*xi)/cosh(kappa) + 3*B*Sc*kappa**3*xi*sqrt(gamma**2 + kappa**2)*sinh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) - 2*B*Sc*kappa**3*cosh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) + 2*B*Sc*kappa**3*cosh(kappa*xi)/cosh(kappa) + B*i*kappa**2*omega*xi*sinh(kappa*xi)/cosh(kappa) - B*i*kappa*omega*xi*sqrt(gamma**2 + kappa**2)*sinh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) + C*Sc*gamma**3*kappa*xi**2*cosh(gamma*xi)/sinh(gamma) - C*Sc*gamma**2*kappa*xi*sqrt(gamma**2 + kappa**2)*sinh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) + C*Sc*gamma**2*kappa*xi*sinh(gamma*xi)/sinh(gamma) + 2*C*Sc*gamma**2*kappa*cosh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) + 2*C*Sc*gamma*kappa**3*xi**2*cosh(gamma*xi)/sinh(gamma) - 4*C*Sc*gamma*kappa*cosh(gamma*xi)/sinh(gamma) - 3*C*Sc*kappa**3*xi*sqrt(gamma**2 + kappa**2)*sinh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) + 2*C*Sc*kappa**3*xi*sinh(gamma*xi)/sinh(gamma) + 2*C*Sc*kappa**3*cosh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) - C*gamma*i*kappa*omega*xi**2*cosh(gamma*xi)/sinh(gamma) + C*i*kappa*omega*xi*sqrt(gamma**2 + kappa**2)*sinh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) - C*i*kappa*omega*xi*sinh(gamma*xi)/sinh(gamma) - 3*D*Sc*gamma**2*kappa**2*xi*sinh(xi*sqrt(gamma**2 + kappa**2))/sinh(sqrt(gamma**2 + kappa**2)) + D*Sc*gamma**2*sqrt(gamma**2 + kappa**2)*cosh(xi*sqrt(gamma**2 + kappa**2))/sinh(sqrt(gamma**2 + kappa**2)) - 3*D*Sc*kappa**4*xi*sinh(xi*sqrt(gamma**2 + kappa**2))/sinh(sqrt(gamma**2 + kappa**2)) + 3*D*Sc*kappa**4*xi*sinh(kappa*xi)/sinh(kappa) - 2*D*Sc*kappa**3*cosh(kappa*xi)/sinh(kappa) + 2*D*Sc*kappa**2*sqrt(gamma**2 + kappa**2)*cosh(xi*sqrt(gamma**2 + kappa**2))/sinh(sqrt(gamma**2 + kappa**2)) - D*i*kappa*omega*cosh(kappa*xi)/sinh(kappa) + D*i*omega*sqrt(gamma**2 + kappa**2)*cosh(xi*sqrt(gamma**2 + kappa**2))/sinh(sqrt(gamma**2 + kappa**2)) - 2*P_1*Sc*kappa**3*xi*sinh(kappa*xi) + P_1*Sc*kappa**2*cosh(kappa*xi))*sin(2*eta*kappa)/(2*Sc)

"""