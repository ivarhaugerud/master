from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, Gamma, kappa, kappa_p, Sc, F0, Pe, rho, rho_p = symbols("omega gamma Gamma, kappa kappa_p Sc F0 Pe rho rho_p")

#define known quantities
B0 = Pe*F0*tanh(gamma)/(gamma*gamma*gamma*rho*rho) + Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*(cosh(rho*xi)/(rho*sinh(rho)) - cosh(gamma*xi)/(gamma*sinh(gamma)))


print("\n\n-----------------------------------------------------------------")
print("CALCULATE B1 low Pe")
RHS = kappa*kappa*xi*diff(B0, xi) - 2*diff(B0, xi, xi)
#my_sol = Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1))*xi*sinh(rho*xi)/sinh(rho) + (Pe*F0*tanh(gamma)/(sinh(gamma)*gamma*gamma*gamma*(Sc-1)))*(xi*kappa*kappa*sinh(gamma*xi)/(gamma*gamma-kappa*kappa-rho*rho) + 2*gamma*cosh(gamma*xi)*(rho*rho-gamma*gamma)/(gamma*gamma-rho*rho-kappa*kappa)**2)
#should_be_zero = simplify(diff(my_sol, xi, xi) - (rho*rho + kappa*kappa)*my_sol + RHS)
#print("\n\n Testing my solution:", simplify(expand(should_be_zero)))

#g0 = -(Pe*F0*tanh(gamma)/(rho_p*sinh(rho_p)*gamma*gamma*gamma*(Sc-1)))*(1+rho/tanh(rho) + kappa*kappa*gamma/(tanh(gamma)*(gamma*gamma-rho_p*rho_p)) + gamma*gamma*(rho*rho-2*gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)))
#B1_extra = -cosh(kappa*xi)*cos(kappa*eta)/sinh(kappa)
B1 = (Pe*F0*tanh(gamma)/(gamma*gamma*gamma*(Sc-1)))*(xi*sinh(rho*xi)/sinh(rho)-cosh(rho_p*xi)/cosh(rho_p) + (kappa*kappa/(gamma*gamma-rho_p*rho_p))*(xi*sinh(gamma*xi)/sinh(gamma) - cosh(rho_p*xi)/cosh(rho_p)) + 2*gamma*(rho*rho-gamma*gamma)/((gamma*gamma-rho_p*rho_p)*(gamma*gamma-rho_p*rho_p)*tanh(gamma))*(cosh(gamma*xi)/cosh(gamma)-cosh(rho_p*xi)/cosh(rho_p)))
#should_be_zero = simplify( -diff(B1, xi, xi) + (rho*rho+kappa*kappa)*B1 - RHS)
#should_be_zero = simplify(expand(should_be_zero.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
#print("Test homogeneous + particular solution: ", should_be_zero)

derivative = simplify(diff(B1, xi))
boundary_value = simplify(derivative.subs(xi, 1))
boundary_value = simplify(expand(boundary_value.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
#boundary_value = simplify(expand(boundary_value.subs(rho, gamma*sqrt(Sc))))

print("\n Test boundary condition particular solution: ", boundary_value)

derivative = simplify(diff(B1_homo, xi))
boundary_value = simplify(derivative.subs(xi, 1))
boundary_value = simplify(expand(boundary_value.subs(rho_p, sqrt(rho*rho+kappa*kappa))))
#boundary_value = simplify(expand(boundary_value.subs(rho, gamma*sqrt(Sc))))

print("\n Test boundary condition: ", boundary_value)
