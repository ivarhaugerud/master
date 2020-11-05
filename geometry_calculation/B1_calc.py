from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, Gamma, kappa, kappa_p, Sc, F0, Pe, rho, rho_p = symbols("omega gamma Gamma, kappa kappa_p Sc F0 Pe rho rho_p")

#define known quantities
B0 = Pe*F0*tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(cosh(rho*xi)/(rho*sinh(rho)) - cosh(gamma*xi)/(gamma*sinh(gamma)))


print("\n\n-----------------------------------------------------------------")
print("CALCULATE B1 low Pe")
RHS = kappa*kappa*xi*diff(B0, xi) - 2*diff(B0, xi, xi)

B1_extra = -cosh(kappa*xi)*cos(kappa*eta)/sinh(kappa)

B1_partic = (Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(xi*sinh(rho*xi)/sinh(rho) + (kappa*kappa/(gamma*gamma-rho_p*rho_p))*(xi*sinh(gamma*xi)/sinh(gamma)) + 2*gamma*(rho*rho-gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)*tanh(gamma))*(cosh(gamma*xi)/cosh(gamma)))
should_be_zero = simplify( -diff(B1_partic, xi, xi) + (rho*rho+kappa*kappa)*B1_partic - RHS)
should_be_zero = simplify(expand(should_be_zero.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
print("Test homogeneous + particular solution: ", should_be_zero)

derivative     = simplify(diff(B1_partic, xi))
boundary_value = simplify(derivative.subs(xi, 1))
boundary_value = simplify(expand(boundary_value.subs(rho_p, sqrt(rho*rho+kappa*kappa))))



g0 = -(Pe*F0*tanh(gamma)/(rho_p*sinh(rho_p)*gamma*gamma*gamma*(Sc-1)))*(1 + rho/tanh(rho) + (kappa*kappa/(gamma*gamma-rho_p*rho_p))*(1+gamma/tanh(gamma)) + 2*gamma*gamma*(rho*rho-gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)))
B1_homo = g0*cosh(rho_p*xi)

derivative = simplify(diff(B1_homo+B1_partic+B1_extra, xi))
boundary_value = simplify(derivative.subs(xi, 1))
boundary_value = simplify(expand(boundary_value.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
print("\n Test boundary condition: ", boundary_value)





g0 -F0*Pe*(rho*(gamma**4 - 2*gamma**2*kappa**2 - 2*gamma**2*rho**2 + kappa**4 + 2*kappa**2*rho**2 + rho**4)*tanh(gamma) - (gamma**4*sqrt(kappa**2 + rho**2)*tanh(gamma)*tanh(sqrt(kappa**2 + rho**2)) + gamma**4*tanh(gamma) - gamma**3*kappa**2 - 2*gamma**3*sqrt(kappa**2 + rho**2)*tanh(sqrt(kappa**2 + rho**2)) - gamma**2*kappa**2*sqrt(kappa**2 + rho**2)*tanh(gamma)*tanh(sqrt(kappa**2 + rho**2)) + gamma**2*kappa**2*tanh(gamma) - 2*gamma**2*rho**2*sqrt(kappa**2 + rho**2)*tanh(gamma)*tanh(sqrt(kappa**2 + rho**2)) + gamma*kappa**4 + gamma*kappa**2*rho**2 + 2*gamma*rho**2*sqrt(kappa**2 + rho**2)*tanh(sqrt(kappa**2 + rho**2)) + kappa**2*rho**2*sqrt(kappa**2 + rho**2)*tanh(gamma)*tanh(sqrt(kappa**2 + rho**2)) - kappa**2*rho**2*tanh(gamma) + rho**4*sqrt(kappa**2 + rho**2)*tanh(gamma)*tanh(sqrt(kappa**2 + rho**2)) - rho**4*tanh(gamma))*tanh(rho))/(gamma**3*rho_p*(Sc*gamma**4 - 2*Sc*gamma**2*kappa**2 - 2*Sc*gamma**2*rho**2 + Sc*kappa**4 + 2*Sc*kappa**2*rho**2 + Sc*rho**4 - gamma**4 + 2*gamma**2*kappa**2 + 2*gamma**2*rho**2 - kappa**4 - 2*kappa**2*rho**2 - rho**4)*sinh(rho_p)*tanh(rho))
