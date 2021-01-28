from sympy import * 

#define variables
xi = symbols("xi")
omega, gamma, kappa, Sc, F0, Pe, rho, i, P_1, kappa_p = symbols("omega gamma kappa Sc F0 Pe rho i P_1 kappa_p")

#define known quantities
ux0  = F0*(1-cosh(gamma*xi)/cosh(gamma))/(2*gamma*gamma)
ux1 = (((P_1*kappa*cosh(kappa)/(gamma*gamma))*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + (F0*tanh(gamma)/gamma)*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma))))/2
uy1 = ((kappa*P_1*sinh(kappa)/(gamma*gamma))*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa)))/2
P1  = (P_1*cosh(kappa*xi))/2


order = 2
B0_deriv = Pe*F0*tanh(gamma)*(sinh(rho*xi)/sinh(rho)-sinh(gamma*xi)/sinh(gamma))/(2*gamma*(rho*rho-gamma*gamma))
B0_deriv = series(B0_deriv, kappa, n=order).removeO()
B0_deriv = simplify(series(B0_deriv, gamma, n=order).removeO())
B0_deriv = simplify(series(B0_deriv, rho,   n=order+2).removeO())

print("B0_deriv_taylor: ", simplify(expand(simplify(B0_deriv))))
print("B0_deriv_deriv_taylor: ", simplify(expand(simplify(diff(B0_deriv, xi)))))


ux0 = simplify(series(ux0, gamma, n=order).removeO())
print("\n ux0 (1): ", ux0)


uy1 = uy1.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
uy1 = series(uy1, kappa_p, n=order).removeO()
uy1 = uy1.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
uy1 = series(uy1, kappa, n=order).removeO()
uy1 = simplify(series(uy1, gamma, n=order).removeO())
print("\n uy (1): ", uy1)

ux1 = ux1.subs(P_1, (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
ux1 = series(ux1, kappa_p, n=order).removeO()
ux1 = ux1.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
ux1 = series(ux1, kappa, n=order).removeO()
ux1 = simplify(series(ux1, gamma, n=order).removeO())
print("\n ux (1): ", ux1)

RHS = simplify(expand(2*diff(B0_deriv, xi) - Pe*ux1))
sol_B1 = Pe*F0*(1-18/(rho*rho) - 9*xi*xi + rho*rho*xi*xi -7*rho*rho/30 - rho*rho*xi*xi*xi*xi/2)/(12*rho*rho)
difference = diff(diff(sol_B1, xi), xi) -rho*rho*sol_B1 - RHS
print("Check B1 particular solution: ", simplify(difference))

homo_sol = 3*Pe*F0*(simplify(series(cosh(rho*xi)/(rho*rho*rho*sinh(rho)), rho, n=3).removeO()))/2
#print(homo_sol)
#print(diff(F0*Pe*(rho**4*(630*xi**4 - 1260*xi**2 + 294+2*xi**6 -31) + 2520*rho**2*(3*xi**2 - 1) + 15120)/(10080*rho**4)+sol_B1, xi).subs(xi, 1))
#sol_B1 += F0*Pe*(rho**6*(-105*xi*xi*xi*xi+147*xi*xi+2*xi**6-31) +rho**4*(630*xi**4 - 1260*xi**2 + 294) + 2520*rho**2*(3*xi**2 - 1) + 15120)/(10080*rho**4)
#new_B = F0*Pe*(-2 - rho*rho*xi*xi/2 + 7*rho*rho/60 + xi**4 * rho*rho/4)/(12*rho*rho)
new_B = F0*Pe*(-2 - rho*rho*xi*xi/2 + 7*rho*rho/60 + xi**4 * rho*rho/4 - xi**4 * rho**4 /8 + 7*xi*xi*rho**4/40 - 31*rho**4/840 + rho**4*xi**6/420)/(12*rho*rho)

print("test boundary: ", diff(new_B, xi).subs(xi, 1))
difference = diff(diff(new_B, xi), xi) -rho*rho*new_B - RHS
print("test solution: ", simplify((difference)))

