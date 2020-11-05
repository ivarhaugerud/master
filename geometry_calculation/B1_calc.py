from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, Gamma, kappa, kappa_p, Sc, F0, Pe, rho, rho_p = symbols("omega gamma Gamma, kappa kappa_p Sc F0 Pe rho rho_p")

#define known quantities
B0 = Pe*F0*tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(cosh(rho*xi)/(rho*sinh(rho)) - cosh(gamma*xi)/(gamma*sinh(gamma)))


print("\n-----------------------------------------------------------------")
print("CALCULATE B1")

#calculate RHS
RHS = kappa*kappa*xi*diff(B0, xi) - 2*diff(B0, xi, xi)

#test particular solution
B1_partic = (Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(xi*sinh(rho*xi)/sinh(rho) + (kappa*kappa/(gamma*gamma-rho_p*rho_p))*(xi*sinh(gamma*xi)/sinh(gamma)) + 2*gamma*(rho*rho-gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)*tanh(gamma))*(cosh(gamma*xi)/cosh(gamma)))
#should_be_zero = simplify( -diff(B1_partic, xi, xi) + (rho*rho+kappa*kappa)*B1_partic - RHS)
#should_be_zero = simplify(expand(should_be_zero.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
#print("Test particular solution: ", should_be_zero)

#test Boundary condition
B1_homo_1 = -(Pe*F0*tanh(gamma)/(rho_p*sinh(rho_p)*gamma*gamma*gamma*(Sc-1)))*(1 + rho/tanh(rho) + (kappa*kappa/(gamma*gamma-rho_p*rho_p))*(1+gamma/tanh(gamma)) + 2*gamma*gamma*(rho*rho-gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)))*cosh(rho_p*xi)
B1_homo_2 = 0#-cos(kappa*eta)*cosh(kappa*xi)/sinh(kappa)
"""
derivative = simplify(diff(B1_homo_1+B1_partic+B1_homo_2, xi))
boundary_value = simplify(derivative.subs(xi, 1))
boundary_value = simplify(expand(boundary_value.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
print("\n Test boundary condition: ", boundary_value)
"""
B1 = simplify(B1_homo_2 + B1_homo_1+ B1_partic)

"""

print("\n-----------------------------------------------------------------")
print("CALCULATE B2")
RHS = -0.5*(kappa*kappa*xi*diff(B1, xi)-2*diff(B1, xi, xi)) - (kappa*kappa*xi*diff(B0, xi) + 3*diff(B0, xi, xi) + kappa*kappa*xi*xi*diff(B0, xi, xi))
#RHS = simplify(expand(RHS.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
RHS = simplify(expand(RHS))

print("\n\n", simplify(RHS))
"""
print("\n-----------------------------------------------------------------")
print("CALCULATE D_parallel")