from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, Gamma, kappa, kappa_p, Sc, F0, Pe, rho, rho_p, i = symbols("omega gamma Gamma, kappa kappa_p Sc F0 Pe rho rho_p i")


B, C, D, P_1 = symbols("B C D P_1")

#define known quantities
B0 = Pe*F0*tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(cosh(rho*xi)/(rho*sinh(rho)) - cosh(gamma*xi)/(gamma*sinh(gamma)))
u0  = F0*(1-cosh(gamma*xi)/cosh(gamma))/(gamma*gamma)
ux1 = (B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))

ux1 = ux1.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
ux1 = simplify(ux1.subs(C, F0*tanh(gamma)/gamma))
ux1 = ux1.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
ux1 = ux1.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))

uy1 = simplify(uy1.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma)))
P1 = P_1*cosh(kappa*xi)

RHS_y = diff(P1, xi)

my_sol_y = uy1 
my_sol_y_check = simplify(diff(my_sol_y, xi, xi) - kappa_p*kappa_p*my_sol_y - RHS_y)
my_sol_y_check = simplify(my_sol_y_check.subs(kappa_p, sqrt(kappa*kappa + gamma*gamma)))
print("y-check: ", my_sol_y_check)


RHS_x = -kappa*kappa*xi*diff(u0, xi) + 2*diff(u0, xi, xi) - P1*kappa

my_sol_x =  P_1*kappa*cosh(kappa)*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p*xi)/cosh(kappa_p))/(gamma*gamma) + F0*tanh(gamma)*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma))/gamma
my_sol_x  = simplify(my_sol_x)

my_sol_check = simplify(diff(my_sol_x, xi, xi) - kappa_p*kappa_p*my_sol_x - RHS_x)
my_sol_check = simplify(expand(my_sol_check))
my_sol_check = simplify(my_sol_check.subs(kappa_p, sqrt(kappa*kappa + gamma*gamma)))


print("\n\n x-check", simplify(expand(my_sol_check)))




