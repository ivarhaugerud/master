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
uy1 = uy1.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))

RHS = -kappa*kappa*xi*diff(B0, xi)# - u0*(Pe + kappa*cosh(kappa*xi)/sinh(kappa))
print(simplify(RHS))

my_sol1 =  kappa*kappa*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)*sinh(rho)) * (xi*sinh(rho*xi)/(rho_p*rho_p-rho*rho) + 2*rho*cosh(rho*xi)/(rho_p*rho_p-rho*rho)**2)
my_sol2 = -kappa*kappa*Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)*sinh(gamma))*(xi*sinh(gamma*xi)/(rho_p*rho_p-gamma*gamma) + 2*gamma*cosh(gamma*xi)/(rho_p*rho_p-gamma*gamma)**2)
my_sol = simplify(my_sol1+my_sol2)
my_sol_check = diff(my_sol, xi, xi) - rho_p*rho_p*my_sol - RHS
my_sol_check = simplify(expand(my_sol_check))

my_sol_check = my_sol_check.subs(gamma, sqrt(i*omega/Sc))
my_sol_check = my_sol_check.subs(rho_p, sqrt(kappa*kappa + i*omega))
my_sol_check = my_sol_check.subs(rho, sqrt(i*omega))
my_sol_check = my_sol_check.subs(i, sqrt(-1))

print("\n\n ", simplify(expand(my_sol_check)))



