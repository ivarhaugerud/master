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
"""
grad1_u1 = simplify(-xi*kappa*cos(kappa*eta)*diff(ux1, xi) - sin(kappa*eta)*diff(uy1, xi))
grad2_u0 = simplify(sin(kappa*eta)*xi*kappa*cos(kappa*eta)*diff(u0, xi))
laplas2_u0x = simplify(xi*kappa*kappa*(3*cos(2*kappa*eta)/2+1/2)*diff(u0, xi) + (3*sin(kappa*eta)*sin(kappa*eta) + xi*xi*kappa*kappa*cos(kappa*eta)*cos(kappa*eta))*diff(u0, xi, xi))
laplas1_u1x = simplify(-2*kappa*xi*cos(kappa*eta)*diff(ux1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(ux1, xi) - 2*sin(kappa*eta)*diff(ux1, xi, xi))
laplas1_u1y = simplify(-2*kappa*xi*cos(kappa*eta)*diff(uy1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(uy1, xi) - 2*sin(kappa*eta)*diff(uy1, xi, xi))
grad1x_P1 = simplify(-kappa*xi*cos(kappa*eta)*diff(P1, xi))
grad1y_P1 = simplify(-sin(kappa*eta)*diff(P1, xi))
grad0_u2 = - simplify(grad1_u1+grad2_u0)
"""

print("\n\n-----------------------------------------------------------------")
print("CALCULATE P2")
P2 = sin(2*kappa*eta)*(P_1*kappa*xi*sinh(kappa*xi)/2 + psi_2*cosh(2*kappa*xi))
print("P2 = ", P2)
print("-----------------------------------------------------------------")
print("\n\n\n")



print("-----------------------------------------------------------------")
print("CALCULATE U2_x")
my_sol_non_eta   = F0*cosh(gamma*xi)/(4*cosh(gamma))*(1-xi*xi) + (P_1*sinh(kappa)/(2*gamma*gamma))*(kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa)+gamma*gamma*cosh(gamma*xi)/cosh(gamma)-kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p))
my_sol_eta = xi*xi*cosh(gamma*xi)*F0/(4*cosh(gamma)) - 2*kappa*psi_2*cosh(2*kappa*xi)/(gamma*gamma) + (P_1*kappa_p*kappa_p*sinh(kappa)/(gamma*gamma*sinh(kappa_p)))*(xi*sinh(kappa_p*xi)/2) - P_1*kappa*kappa*xi*sinh(kappa*xi)/(2*gamma*gamma)
print("ux(xi, eta) = ", my_sol_eta, "\n")
print("ux(xi) = ", my_sol_non_eta)
print("-----------------------------------------------------------------")
print("\n\n\n")

print("-----------------------------------------------------------------")
print("CALCULATE U2_y")
my_sol_y = P_1*kappa*kappa_p*sinh(kappa)*xi*cosh(kappa_p*xi)/(2*gamma*gamma*sinh(kappa_p)) - P_1*kappa*kappa*xi*cosh(kappa*xi)/(2*gamma*gamma) - 2*kappa*psi_2*sinh(2*kappa*xi)/(gamma*gamma)
print("uy = ", my_sol_y)
print("\n\n\n")
print("-----------------------------------------------------------------")
print("CALCULATE undetermined coeffs")

kappa_pp = sqrt(gamma*gamma+4*kappa*kappa)
Ax, Ay = symbols("Ax Ay")
BCy = -kappa*F0*tanh(gamma)/(2*gamma) -2*kappa*psi_2*sinh(2*kappa)/(gamma*gamma)+Ay*sinh(kappa_pp)
BCy = -P_1*kappa*kappa*cosh(kappa)/(2*gamma*gamma) + P_1*kappa*kappa_p*sinh(kappa)/(2*gamma*gamma*tanh(kappa_p)) - 2*kappa*psi_2*sinh(2*kappa)/(gamma*gamma) + Ay*sinh(kappa_pp)
BCx = F0/4 + P_1*sinh(kappa)/2 + Ax*cosh(kappa_pp) - 2*kappa*psi_2*cosh(2*kappa)/(gamma*gamma)
divergence_condition = Ay*kappa_pp -2*kappa*Ax 
system = [BCy, BCx, divergence_condition]

print("[Ax, Ay, psi_2] = ", linsolve(system, [Ax, Ay, psi_2]))

print("\n\n\n")
print("-----------------------------------------------------------------")
print("Check solution with BCs")
sol_x_no_eta_BC = P_1*sinh(kappa)*(kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa) - kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p) + gamma*gamma*cosh(gamma*xi)/cosh(gamma))/(2*gamma*gamma) + cosh(gamma*xi)*F0*(1-xi*xi)/(4*cosh(gamma))

P_2 = symbols("P_2")
P2 = (P_1*sqrt(gamma**2 + 4*kappa**2)*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(sqrt(gamma**2 + 4*kappa**2)) + gamma**2*(F0 + 2*P_1*sinh(kappa))*sinh(sqrt(gamma**2 + 4*kappa**2))*tanh(kappa_p))/(4*(2*kappa*sinh(sqrt(gamma**2 + 4*kappa**2))*cosh(2*kappa) - sqrt(gamma**2 + 4*kappa**2)*sinh(2*kappa)*cosh(sqrt(gamma**2 + 4*kappa**2)))*tanh(kappa_p))
sol_x_eta = P_1*sinh(kappa)*(kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p) -kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa) - gamma*gamma*cosh(kappa_pp*xi)/cosh(kappa_pp))/(2*gamma*gamma) + F0*(xi*xi*cosh(gamma*xi)/cosh(gamma) - cosh(kappa_pp*xi)/cosh(kappa_pp))/4 - P_2*cosh(2*kappa)*(cosh(2*kappa*xi)/cosh(2*kappa) - cosh(kappa_pp*xi)/cosh(kappa_pp))
sol_y = (P_1*kappa*sinh(kappa)/(2*gamma*gamma*tanh(kappa_p)))*(kappa_p*xi*cosh(kappa_p*xi)/cosh(kappa_p) - kappa_p*sinh(kappa_pp*xi)/sinh(kappa_pp) -kappa*(tanh(kappa_p)/tanh(kappa))*(xi*cosh(kappa*xi)/cosh(kappa) - sinh(kappa_pp*xi)/sinh(kappa_pp))) - P_2*sinh(2*kappa)*(sinh(2*kappa*xi)/sinh(2*kappa) - sinh(kappa_pp*xi)/sinh(kappa_pp))


print("\n\n\n")
print("-----------------------------------------------------------------")
print("series expansion in kappa and gamma")
order = 2

my_sol_non_eta = my_sol_non_eta.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
my_sol_non_eta = series(my_sol_non_eta, kappa_p, n=order).removeO()
my_sol_non_eta = my_sol_non_eta.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
my_sol_non_eta = simplify(series(my_sol_non_eta, kappa, n=order).removeO())
my_sol_non_eta = simplify(series(my_sol_non_eta, gamma, n=order).removeO())
print("taylor ux no eta: ", my_sol_non_eta)
#print(simplify(my_sol_non_eta.subs(xi, 1)))


sol_x_eta = sol_x_eta.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
sol_x_eta = series(sol_x_eta, kappa_p, n=order).removeO()
#sol_x_eta = sol_x_eta.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
sol_x_eta = simplify(series(sol_x_eta, kappa, n=order).removeO())
sol_x_eta = simplify(series(sol_x_eta, gamma, n=order).removeO())
print("\n\n taylor ux with eta: ", sol_x_eta)
#print(simplify(sol_x_eta.subs(xi, 1)))


sol_y = sol_y.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
sol_y = series(sol_y, kappa_p, n=order).removeO()
sol_y = sol_y.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
sol_y = simplify(series(sol_y, kappa, n=order).removeO())
sol_y = simplify(series(sol_y, gamma, n=order).removeO())
print("\n\n taylor uy: ", sol_y)
#print(simplify(sol_y.subs(xi, 1)))

#gamma**2*(kappa**2*(105*kappa**2*xi**6 - 1470*kappa**2*xi**4 + 837*kappa**2*xi**2 - 16*kappa**2 + 1050*xi**4 - 3360*xi**2 + 1470) + 25200) + kappa**2*(4200*kappa**2*xi**4 - 2520*kappa**2*xi**2 + 10500*xi**2 + 2100)