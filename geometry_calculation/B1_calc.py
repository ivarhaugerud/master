from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, Gamma, kappa, kappa_p, Sc, F0, Pe, rho, rho_p = symbols("omega gamma Gamma, kappa kappa_p Sc F0 Pe rho rho_p")


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


print(simplify((kappa*xi*u0-uy1).subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))))


#-F0*(sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))/sinh(sqrt(gamma**2 + kappa**2)) - sinh(kappa*xi))*tanh(gamma)/(gamma*(1 - sqrt(gamma**2 + kappa**2)*tanh(kappa)/(kappa*tanh(sqrt(gamma**2 + kappa**2))))*cosh(kappa)) + F0*kappa*xi*(1 - cosh(gamma*xi)/cosh(gamma))/gamma**2

"""
print(simplify(series(ux1, kappa, x0=0)), "\n\n")
print(simplify(series(uy1, kappa, x0=0)))
"""
"""
print("\n-----------------------------------------------------------------")
print("CALCULATE B1")

sin_RHS = simplify(kappa*kappa*xi*diff(B0, xi) - 2*diff(B0, xi, xi) + Pe*ux1)
#cos_RHS = simplify(simplify((u0*kappa*xi - Pe*uy1*diff(B0, xi)))*diff(B0, xi))

sin_RHS = sin_RHS.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
sin_RHS = sin_RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
sin_RHS = simplify(expand(sin_RHS))
print(sin_RHS, "whole thing\n")
"""

"""
RHS_with_Pe_no_sin = simplify(RHS.subs(eta, 0))
RHS_with_Pe_no_cos = simplify(RHS.subs(eta, pi/(2*kappa)))

print(RHS_with_Pe_no_cos, " whole thing times Pe*sin(kappa*eta)\n")
print(RHS_with_Pe_no_sin, " whole thing times Pe*cos(kappa*eta)\n")
"""

#calculate RHS
"""
RHS = kappa*kappa*xi*diff(B0, xi) - 2*diff(B0, xi, xi)

#test particular solution
B1_partic = (Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(xi*sinh(rho*xi)/sinh(rho) + (kappa*kappa/(gamma*gamma-rho_p*rho_p))*(xi*sinh(gamma*xi)/sinh(gamma)) + 2*gamma*(rho*rho-gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)*tanh(gamma))*(cosh(gamma*xi)/cosh(gamma)))
#should_be_zero = simplify( -diff(B1_partic, xi, xi) + (rho*rho+kappa*kappa)*B1_partic - RHS)
#should_be_zero = simplify(expand(should_be_zero.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
#print("Test particular solution: ", should_be_zero)

#test Boundary condition
B1_homo_1 = -(Pe*F0*tanh(gamma)/(rho_p*sinh(rho_p)*gamma*gamma*gamma*(Sc-1)))*(1 + rho/tanh(rho) + (kappa*kappa/(gamma*gamma-rho_p*rho_p))*(1+gamma/tanh(gamma)) + 2*gamma*gamma*(rho*rho-gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)))*cosh(rho_p*xi)
B1_homo_2 = 0#-cos(kappa*eta)*cosh(kappa*xi)/sinh(kappa)

derivative = simplify(diff(B1_homo_1+B1_partic+B1_homo_2, xi))
boundary_value = simplify(derivative.subs(xi, 1))
boundary_value = simplify(expand(boundary_value.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
print("\n Test boundary condition: ", boundary_value)
B1 = simplify(B1_homo_2 + B1_homo_1+ B1_partic)


print("\n-----------------------------------------------------------------")
print("CALCULATE B2")
RHS = -0.5*(kappa*kappa*xi*diff(B1, xi)-2*diff(B1, xi, xi)) - (kappa*kappa*xi*diff(B0, xi) + 3*diff(B0, xi, xi) + kappa*kappa*xi*xi*diff(B0, xi, xi))
#RHS = simplify(expand(RHS.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
RHS = simplify(expand(RHS))

print("\n\n", simplify(RHS))
"""


#-F0*Pe*(Pe - 1)*((2*gamma*cosh(gamma*xi) - kappa**2*xi*sinh(gamma*xi))*sinh(rho) + (kappa**2*xi*sinh(rho*xi) - 2*rho*cosh(rho*xi))*sinh(gamma))*tanh(gamma)/(gamma**3*(Sc - 1)*sinh(gamma)*sinh(rho))  whole thing times Pe*sin(kappa*eta)
