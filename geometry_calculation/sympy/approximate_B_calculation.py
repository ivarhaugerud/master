from sympy import * 

#define variables
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i, P_1, rho_p = symbols("omega gamma kappa Sc F0 Pe rho i P_1 rho_p")
order = 2 

#define known quantities
ux0 = F0*(1-xi*xi)/4  # eta indep 
uy1 = F0*kappa*xi*(1-xi*xi)/4  #sin
ux1 = -F0*(1-xi*xi)/4 #sin
B0 = F0*Pe*(gamma**2*rho**2*(630*xi**4 - 1260*xi**2 + 294) + 5040*gamma**2 - 15120)/(30240*gamma**2*rho**2)
B0_deriv = F0*Pe*xi*(xi*xi-1)/12
B0_deriv_deriv = F0*Pe*(3*xi*xi-1)/12 



B0_grad_approx = F0*Pe*(sinh(rho*xi)/sinh(rho)-1)/(2*rho*rho)
print(simplify(series(B0_grad_approx, rho, n=2)))


RHS = -2*B0_deriv_deriv + ux1*Pe

print(simplify(RHS))
#beta1 = -rho*rho*xi**4/(24*rho_p*rho_p) + xi*xi*(rho*rho - 3 - 6*rho*rho/(rho_p*rho_p))/(12*rho_p*rho_p)  - (12*rho*rho/(rho_p**4) - 2*rho*rho/(rho_p*rho_p) + 6/(rho_p*rho_p) + 1 + 7*rho*rho/30)/(12*rho_p*rho_p)
#beta1 = -rho**2 * xi**4/2 + xi*xi*(rho*rho - 3 - 6*rho*rho/(rho_p*rho_p) )
#beta1 += -6 - 12*rho*rho/(rho_p*rho_p)
#beta1 *= Pe*F0/(12*rho_p*rho_p)
print("CHECK BETA1 SOL: ", simplify(  (diff(diff(beta1, xi), xi) - rho_p*rho_p*beta1 + RHS)   ))

"""
beta1 = F0*Pe*(-7 + 30*xi*xi - 15*xi**4 + (30-270*xi*xi)/(rho*rho) -540 /(rho**4))/360
beta1 += 3*F0*Pe*cosh(rho*xi)/(2*rho*rho*rho*sinh(rho))
#beta1 = simplify(series(beta1,  ,   n=order+2).removeO())
print("BETA1: ", simplify(RHS))

print("CHECK BETA1 SOL: ", simplify(  series(diff(diff(beta1, xi), xi) - rho*rho*beta1 + RHS, rho, n=order+2).removeO()   ))
print("CHECK BOUNDARY BETA1 SOL: ", simplify(diff(beta1, xi).subs(xi, 1)), simplify(diff(beta1, xi).subs(xi, -1)))
"""
"""
print("\n\n")
RHS_beta0 = -Pe*Pe*F0*F0*kappa*(1-xi*xi)/(24*rho*rho)
beta0 = Pe*Pe*F0*F0*kappa*(xi**4 - 2*xi**2 + 7/15 + 8/(kappa*kappa) )/(288*rho*rho) + kappa*(1/6 - xi*xi/2 -1/(kappa*kappa))
print("TEST SOLUTION BETA0: ", series(simplify(diff(diff(beta0, xi), xi) - kappa*kappa*beta0 - RHS_beta0), kappa, n=2).removeO())
print("TEST BCS BETA0: ", simplify(diff(beta0, xi).subs(xi,  1)), simplify(diff(beta0, xi).subs(xi, -1)))
print("\n\n")
"""
"""

RHS_beta2 = -Pe*Pe*F0*F0*kappa*(1-xi*xi)/(24*rho*rho)
beta2 = Pe*Pe*F0*F0*kappa*(1/2 - xi*xi/2 - 1/(2*rho**2)  + cosh(rho*xi*sqrt(2))/(rho*sqrt(2)*sinh(sqrt(2)*rho)))/(24*rho*rho*rho*rho)
#beta2 = Pe*Pe*F0*F0*kappa*(1/2 - xi*xi/2 - 1/(2*rho**2)  -1/6 + xi*xi/2 + 1/(2*rho*rho) + rho*rho*(15*xi**4 - 30*xi*xi + 7)/180 )/(24*rho**4)
#beta2 += F0*F0*kappa*Pe*Pe*cosh(rho*sqrt(2)*xi)/(24*sqrt(2)*sinh(rho*sqrt(2))*rho**5)
#beta2 = 

#print(simplify(series(cosh(sqrt(2)*rho*xi)/(sqrt(2)*rho*sinh(sqrt(2)*rho)), rho, n=8)))
print(simplify(series(beta2, rho, n=4)))

print("TEST SOLUTION BETA0: ", simplify(diff(diff(beta2, xi), xi) - 2*rho*rho*beta2 - RHS_beta2))
print("TEST BCS BETA0: ", simplify(diff(beta2, xi).subs(xi,  1)), simplify(diff(beta2, xi).subs(xi, -1)))
#print("\n\n")


#(F0**2*Pe**2*kappa*xi**4/288 - F0**2*Pe**2*kappa*xi**2/144 + 7*F0**2*Pe**2*kappa/4320)/rho**2 + rho**2*(F0**2*Pe**2*kappa*xi**8/120960 - F0**2*Pe**2*kappa*xi**6/12960 + 7*F0**2*Pe**2*kappa*xi**4/25920 - 31*F0**2*Pe**2*kappa*xi**2/90720 + 127*F0**2*Pe**2*kappa/1814400) + 0.0138888888888889*F0**2*Pe**2*kappa/rho**4 - 31*F0**2*Pe**2*kappa/90720 + 7*F0**2*Pe**2*kappa*xi**2/4320 - F0**2*Pe**2*kappa*xi**4/864 + F0**2*Pe**2*kappa*xi**6/4320 + O(rho**4)
"""
"""
rho_p = sqrt(kappa*kappa + 2*rho*rho)
beta2 = Pe*Pe*F0*F0*kappa*( (1-xi*xi)/(rho_p*rho_p) - 2/(rho_p**4) + 2*cosh(rho_p*xi)/(rho_p*rho_p*rho_p*sinh(rho_p)))/(24*rho*rho)
RHS = -Pe*Pe*F0*F0*kappa*(1-xi*xi)/(24*rho*rho)

test = diff(diff(beta2, xi), xi) - rho_p*rho_p*beta2 -RHS

test = series(beta2, rho_p, n=4).removeO()
test = series(test, rho, n=4).removeO()
print(simplify(test))
"""