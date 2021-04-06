from sympy import * 

#define variables
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i, P_1, rho_p = symbols("omega gamma kappa Sc F0 Pe rho i P_1 rho_p")
order = 2 

#define known quantities
ux0 = F0*(1-xi*xi)/4  # eta indep 
uy1 = F0*kappa*xi*(1-xi*xi)/4  #cos
ux1 = -F0*(1-xi*xi)/4 #sin
ux2 = 5*F0*(xi**2 - 1)/8
ux2_avg = -5*F0/12 
ux1_avg = -F0/6


B0 = F0*Pe*cosh(rho*xi)/(2*rho*rho*rho*sinh(rho)) - F0*Pe*xi*xi/(4*rho*rho)
B0_deriv = diff(B0, xi)
B0_deriv_deriv = diff(B0_deriv, xi)



RHS = -2*B0_deriv_deriv + ux1*Pe

beta1 = 3/(2*rho**4) + (xi*xi - 1)/(4*rho*rho) + (xi*sinh(rho*xi) - cosh(rho*xi)/tanh(rho) - 2*cosh(rho*xi)/rho )/(2*rho*rho*sinh(rho))
beta1 *= Pe*F0 

print("\n\n BETA 1")
print("Solution: ", simplify( diff(diff(beta1, xi), xi) - rho*rho*beta1 + RHS) )
print("BC: ",  simplify(diff(beta1, xi).subs(xi, 1)))
print("\n\n")

beta0_rhs = (1-xi*xi)/(rho*rho) + (2*xi*xi-1-xi**4)/2 + ((1-xi*xi)/sinh(rho))*(xi*sinh(rho*xi)  - cosh(rho*xi)*(1/tanh(rho) + 2/rho))
beta0_rhs *= Pe*Pe*F0*F0*kappa/(8*rho*rho)
beta0_rhs_num = Pe*ux0*kappa*beta1 


print("\n\n BETA 0")
beta0 = 3*xi*xi*(1-xi*xi/6)/(2*rho*rho) + xi*xi*(-1+xi*xi/3 - xi**4/15)/4 + cosh(rho*xi)*(24/(rho*rho*rho) + 6*xi*xi/rho - 2/rho)/(rho*rho*sinh(rho)) + xi*sinh(rho*xi)*(-18/(rho*rho) - xi*xi + 1)/(rho*rho*sinh(rho)) + cosh(rho*xi)*(2/rho + 1/tanh(rho))*(6/(rho*rho) + xi*xi - 1 - 4*xi*tanh(rho*xi)/rho)/(rho*rho*sinh(rho))
beta0 *= Pe*Pe*F0*F0*kappa/(8*rho*rho)

print("Solution: ", simplify(diff(diff(beta0, xi), xi) - beta0_rhs_num))
beta0 += cosh(kappa*xi)*(-1 + F0*F0*Pe*Pe/(60*rho**6)*(2*rho**4 - 30*rho**2 + 15*rho**2/(tanh(rho)**2) + 60*rho/tanh(rho) - 75) )/sinh(kappa)
print("BC: ", simplify(diff(beta0, xi).subs(xi, 1)))
print("\n\n")


beta2 = (1/4 - 2/(rho*rho) + 3/(rho**4) + xi*xi*(3/(rho*rho) - 1/2) + xi*xi*xi*xi/4) + (cosh(rho*xi)*(72/(rho*rho) + 6*xi*xi - 2)/rho + sinh(rho*xi)*(30*xi/(rho*rho) + xi**3 - xi) )/sinh(rho) + (2/rho + 1/tanh(rho))*(cosh(rho*xi)*(1-xi*xi - 10/(rho*rho)) - 4*xi*sinh(rho*xi)/rho)/sinh(rho)
beta2 *= Pe*Pe*F0*F0*kappa/(8*rho**4)

print("here: ", simplify(diff(beta2, xi).subs(xi, 1)))
beta2 -= (F0*F0*Pe*Pe*kappa/(4*sqrt(2)*sinh(sqrt(2)*rho)*rho**7))*(40 + 8*rho/tanh(rho) - 3*rho*rho/(sinh(rho)**2))*cosh(sqrt(2)*rho*xi)

beta_2_rhs = 3*(1-xi*xi)/(rho*rho) + (2*xi*xi-1-xi**4)/2 + ((1-xi*xi)/sinh(rho))*(xi*sinh(rho*xi)   - cosh(rho*xi)*(1/tanh(rho) + 2/rho))
beta_2_rhs *= Pe*Pe*F0*F0*kappa/(8*rho**2)
beta2_rhs_num = Pe*ux0*kappa*beta1 

print("the two source terms: ", simplify(beta_2_rhs-beta2_rhs_num))
print("\n\n BETA 2")
print("Solution: ", simplify(diff(diff(beta2, xi), xi) -2*rho*rho*beta2 - beta_2_rhs))
print("BC :", simplify(diff(beta2, xi).subs(xi, 1)))
print("\n\n")





rhs_beta2 = 3*diff(diff(B0, xi), xi)/2 - diff(diff(beta1, xi), xi) + Pe*(ux2 - ux1_avg/2-ux2_avg)

B2 = 5*xi*xi - 1 -3*xi*sinh(rho*xi)/sinh(rho) - 2*cosh(rho)*rho*xi*sinh(rho*xi)/(sinh(rho)*sinh(rho)) + rho*xi*xi*cosh(rho*xi)/sinh(rho) - xi*sinh(rho*xi)/sinh(rho)
B2 *= F0*Pe/(8*rho*rho)
B2 -= F0*Pe*cosh(rho*xi)*(rho*rho + 6 - 2*rho*rho/(tanh(rho)**2) -4*rho/tanh(rho) )/(8*rho*rho*rho*sinh(rho))
#B2 -= 3*F0*Pe/(4*rho*rho*rho*rho)


print(simplify( diff(diff(B2, xi), xi) - rho*rho*B2 - (-rhs_beta2)))
print(simplify(diff(B2, xi).subs(xi, 1)))





