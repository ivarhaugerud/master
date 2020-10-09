from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
B, C, D, P_1, F0, i, psi_2 = symbols("B C D P_1 F0 i psi_2")

#define known quantities
u0  = F0*(1-cosh(gamma*xi)/cosh(gamma))/(gamma*gamma)
ux1 = sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)*cos(kappa*eta)

grad1_u1 = simplify(-xi*kappa*cos(kappa*eta)*diff(ux1, xi) - sin(kappa*eta)*diff(uy1, xi))
grad2_u0 = simplify(sin(kappa*eta)*xi*kappa*cos(kappa*eta)*diff(u0, xi))
laplas2_u0x = simplify(xi*kappa*kappa*(3*cos(2*kappa*eta)/2+1/2)*diff(u0, xi) + (3*sin(kappa*eta)*sin(kappa*eta) + xi*xi*kappa*kappa*cos(kappa*eta)*cos(kappa*eta))*diff(u0, xi, xi))
laplas1_u1x = simplify(-2*kappa*xi*cos(kappa*eta)*diff(ux1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(ux1, xi) - 2*sin(kappa*eta)*diff(ux1, xi, xi))
laplas1_u1y = simplify(-2*kappa*xi*cos(kappa*eta)*diff(uy1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(uy1, xi) - 2*sin(kappa*eta)*diff(uy1, xi, xi))
grad1x_P1 = simplify(-kappa*xi*cos(kappa*eta)*diff(P1, xi))
grad1y_P1 = simplify(-sin(kappa*eta)*diff(P1, xi))
grad0_u2 = - simplify(grad1_u1+grad2_u0)


print("CALCULATE P2")
print("-----------------------------------------------------------------")
term1 = -gamma*gamma*(grad0_u2) + diff(grad0_u2, xi, xi) + diff(grad0_u2, eta, eta)
term2 =  diff(laplas2_u0x, eta)
term3 =  diff(laplas1_u1x, eta) + diff(laplas1_u1y, xi)
term4 = -diff(grad1x_P1, eta) - diff(grad1y_P1, xi)

RHS = term1 + term2 + term3 + term4

RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, F0*tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
RHS = RHS.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))

RHS = simplify(expand(RHS/(sin(2*kappa*eta))))
my_RHS = expand(P_1*kappa*kappa*(cosh(kappa*xi)-3*kappa*xi*sinh(kappa*xi)/2))
difference = RHS-my_RHS
difference = difference.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference = difference.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
difference = simplify(expand(difference))
print("Difference my RHS and proper RHS: ", difference)

my_sol_P2 = sin(2*kappa*eta)*(P_1*kappa*xi*sinh(kappa*xi)/2 + psi_2*cosh(2*kappa*xi))
difference_sol_P2 = diff(my_sol_P2, xi, xi) + diff(my_sol_P2, eta, eta) - my_RHS*sin(2*kappa*eta)
difference_sol_P2 = difference_sol_P2.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference_sol_P2 = difference_sol_P2.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
difference_sol_P2 = simplify(expand(difference_sol_P2))
print("Plugging back in my solution: ", difference_sol_P2)

P2 = my_sol_P2
print("P2 = ", P2)
print("-----------------------------------------------------------------")
print("\n\n\n")
print("CALCULATE U2_x")
print("-----------------------------------------------------------------")
#define all the terms
term1 = laplas2_u0x
term2 = laplas1_u1x
term3 = -grad1x_P1
term4 = -diff(P2, eta)
RHS = term1 + term2 + term3 + term4

RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, F0*tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
RHS = RHS.subs(F0, P_1*(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p)))*kappa*cosh(kappa)/(gamma*tanh(gamma)))
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))

RHS = simplify(expand(RHS))
RHS_no_eta = simplify(RHS.subs(eta, pi/(4*kappa)))
RHS_just_eta = simplify(RHS - RHS_no_eta)
print("\n\n no eta RHS: ", RHS_no_eta)
print("\n\n just eta RHS: ", RHS_just_eta)

my_RHS = F0*(gamma*xi*sinh(gamma*xi) + cosh(gamma*xi)/2)/cosh(gamma) + P_1*kappa*kappa*((gamma*gamma-kappa*kappa)*xi*sinh(kappa*xi) - 2*kappa*cosh(kappa*xi))/(2*gamma*gamma) + P_1*kappa*kappa*kappa_p*kappa_p*sinh(kappa)*(xi*sinh(kappa_p*xi) + 2*kappa_p*cosh(kappa_p*xi)/(kappa*kappa))/(2*gamma*gamma*sinh(kappa_p))
difference = my_RHS- RHS_no_eta

difference = difference.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference = difference.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
difference = simplify(expand(difference))

print("Difference my RHS and theirs (no eta)", difference)

my_sol_non_eta = -F0*xi*xi*cosh(gamma*xi)/(4*cosh(gamma)) + P_1*kappa*kappa*xi*sinh(kappa*xi)/(2*gamma*gamma) - P_1*kappa_p*kappa_p*sinh(kappa)*xi*sinh(kappa_p*xi)/(2*gamma*gamma*sinh(kappa_p))
difference_sol_x = diff(my_sol_non_eta, xi, xi) - gamma*gamma*my_sol_non_eta+ my_RHS
difference_sol_x = difference_sol_x.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference_sol_x = difference_sol_x.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
difference_sol_x = simplify(expand(difference_sol_x))
print(difference_sol_x)

my_RHS = -gamma*F0*xi*sinh(gamma*xi)/cosh(gamma) + F0*cosh(gamma*xi)*(kappa*kappa*xi*xi-1/2)/cosh(gamma) - P_1*xi*sinh(kappa*xi)*(3*kappa**4 + gamma*gamma*kappa*kappa)/(2*gamma*gamma) + P_1*kappa*kappa*kappa*cosh(kappa*xi)/(gamma*gamma) - 2*kappa*psi_2*cosh(2*kappa*xi) + 3*P_1*kappa*kappa*kappa_p*kappa_p*sinh(kappa)*xi*sinh(kappa_p*xi)/(2*gamma*gamma*sinh(kappa_p)) - P_1*kappa_p*kappa_p*kappa_p*cosh(kappa_p*xi)*sinh(kappa)/(gamma*gamma*sinh(kappa_p))
my_RHS *= cos(2*kappa*eta)
difference = my_RHS- RHS_just_eta
difference = difference.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference = difference.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
difference = simplify(expand(difference))
print("Difference my RHS and theirs (just eta)", difference)

my_sol_eta = xi*xi*cosh(gamma*xi)*F0/(4*cosh(gamma)) + P_1*kappa*kappa*kappa*cosh(kappa*xi)/(gamma*gamma*(gamma*gamma+3*kappa*kappa)) - 2*kappa*psi_2*cosh(2*kappa*xi)/(gamma*gamma)  + (P_1*kappa_p*kappa_p*sinh(kappa)/(gamma*gamma*sinh(kappa_p)))*(xi*sinh(kappa_p*xi)/2) + P_1*kappa*kappa*(xi*sinh(kappa*xi)/2 + kappa*cosh(kappa*xi)/(gamma*gamma+3*kappa*kappa))/(gamma*gamma)
difference_sol_x = diff(my_sol_eta, xi, xi) - (gamma*gamma+4*kappa*kappa)*my_sol_eta + my_RHS/cos(2*kappa*eta)

#difference_sol_x = difference_sol_x.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference_sol_x = difference_sol_x.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
difference_sol_x = simplify(expand(difference_sol_x))
print(difference_sol_x)

#-F0*kappa**2*xi*(3*gamma**2*kappa*sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))*tanh(sqrt(gamma**2 + kappa**2)) + 2*gamma**2*kappa*sinh(kappa*xi)*sinh(sqrt(gamma**2 + kappa**2))*tanh(sqrt(gamma**2 + kappa**2)) - 3*gamma**2*sqrt(gamma**2 + kappa**2)*sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))*tanh(kappa) - 2*gamma**2*sqrt(gamma**2 + kappa**2)*sinh(kappa*xi)*sinh(sqrt(gamma**2 + kappa**2))*tanh(kappa) + 3*kappa**3*sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))*tanh(sqrt(gamma**2 + kappa**2)) + 6*kappa**3*sinh(kappa*xi)*sinh(sqrt(gamma**2 + kappa**2))*tanh(sqrt(gamma**2 + kappa**2)) - 3*kappa**2*sqrt(gamma**2 + kappa**2)*sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))*tanh(kappa) - 6*kappa**2*sqrt(gamma**2 + kappa**2)*sinh(kappa*xi)*sinh(sqrt(gamma**2 + kappa**2))*tanh(kappa))*tanh(gamma)*tanh(sqrt(gamma**2 + kappa**2))/(2*gamma*(gamma**2*tanh(kappa)**2 + kappa**2*tanh(kappa)**2 + kappa**2*tanh(sqrt(gamma**2 + kappa**2))**2 - 2*kappa*sqrt(gamma**2 + kappa**2)*tanh(kappa)*tanh(sqrt(gamma**2 + kappa**2)))*sinh(sqrt(gamma**2 + kappa**2))*cosh(kappa))
