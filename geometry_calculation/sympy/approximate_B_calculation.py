from sympy import * 

#define variables
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i, P_1, kappa_p = symbols("omega gamma kappa Sc F0 Pe rho i P_1 kappa_p")
order = 2 

#define known quantities
ux0 = F0*(1-xi*xi)/4  # eta indep 
uy1 = F0*kappa*xi*(1-xi*xi)/4  #sin
ux1 = -F0*(1-xi*xi)/4 #sin
B0 = F0*Pe*(gamma**2*rho**2*(rho**2*(21*xi**6 - 105*xi**4 + 147*xi**2 - 31) + 630*xi**4 - 1260*xi**2 + 294) + 5040*gamma**2 - 15120)/(30240*gamma**2*rho**2)
B0_deriv = F0*Pe*xi*(60*(xi*xi-1) + rho*rho*(3*xi**4 - 10*xi*xi +7))/720
B0_deriv_deriv = F0*Pe*(180*xi*xi-60+rho*rho*(15*xi**4 -30*xi**2 + 7))/720 

RHS = -2*B0_deriv_deriv + ux1*Pe
beta1 = F0*Pe*(-7 + 30*xi*xi - 15*xi**4 + (30-270*xi*xi)/(rho*rho) -540 /(rho**4))/360
beta1 += 3*F0*Pe*cosh(rho*xi)/(2*rho*rho*rho*sinh(rho))
beta1 = simplify(series(beta1, rho,   n=order+2).removeO())


print("CHECK BETA1 SOL: ", simplify(  series(diff(diff(beta1, xi), xi) - rho*rho*beta1 + RHS, rho, n=order+2).removeO()   ))
print("CHECK BOUNDARY BETA1 SOL: ", simplify(diff(beta1, xi).subs(xi, 1)), simplify(diff(beta1, xi).subs(xi, -1)))

print("\n\n")
RHS_beta0 = -Pe*Pe*F0*F0*kappa*(1-xi*xi)/(24*rho*rho)
beta0 = Pe*Pe*F0*F0*kappa*(xi**4 - 2*xi**2 + 7/15 + 8/(kappa*kappa) )/(288*rho*rho) + kappa*(1/6 - xi*xi/2 -1/(kappa*kappa))
print("TEST SOLUTION BETA0: ", series(simplify(diff(diff(beta0, xi), xi) - kappa*kappa*beta0 - RHS_beta0), kappa, n=2).removeO())
print("TEST BCS BETA0: ", simplify(diff(beta0, xi).subs(xi,  1)), simplify(diff(beta0, xi).subs(xi, -1)))
print("\n\n")



RHS_beta2 = -Pe*Pe*F0*F0*kappa*(1-xi*xi)/(24*rho*rho)
beta2 = Pe*Pe*F0*F0*kappa*(1/(2*rho*rho) - xi*xi/(2*rho*rho) - 1/(2*rho**4) )/(24*rho*rho)
beta2 += F0*F0*kappa*Pe*Pe*cosh(rho*sqrt(2)*xi)/(24*sqrt(2)*sinh(rho*sqrt(2))*rho**5)
beta2 = 
print(simplify(series(cosh(sqrt(2)*rho*xi)/(sqrt(2)*rho*sinh(sqrt(2)*rho)), rho, n=4)))
print("TEST SOLUTION BETA0: ", simplify(diff(diff(beta2, xi), xi) - 2*rho*rho*beta2 - RHS_beta2))
print("TEST BCS BETA0: ", simplify(diff(beta2, xi).subs(xi,  1)), simplify(diff(beta2, xi).subs(xi, -1)))
#print("\n\n")

(rho**2*(-30 + 90*xi**2 + rho**2*(15*xi**4 - 30*xi**2 + 7) + O(rho**4)) + 90)/(180*rho**2)

