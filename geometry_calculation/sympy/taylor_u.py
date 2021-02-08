from sympy import * 
order = 4

#define variables
xi = symbols("xi")
omega, gamma, kappa, kappa_p, Sc, rho = symbols("omega gamma kappa kappa_p Sc rho", real=True)
P_1, F0, D = symbols("P_1, F0 D", real=True)

#define known quantities
u0  = F0*(1-cosh(gamma*xi)/cosh(gamma))/(2*gamma*gamma)
ux1 = ((P_1*kappa*cosh(kappa)/(gamma*gamma))*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + (F0*tanh(gamma)/gamma)*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))/2
uy1 = ((kappa*P_1*sinh(kappa)/(gamma*gamma))*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa)))/2
P1  = P_1*cosh(kappa*xi)/2

P = (F0*gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p)))
P = series(P, kappa_p, n=order).removeO()
P = P.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
P = series(P, kappa, n=order).removeO()
P = simplify(series(P, gamma, n=order).removeO())
print("P (1): ", P)

order = 4
A = gamma*cosh(gamma*xi)/sinh(gamma)
A = simplify(series(A, gamma, n=order).removeO())
print(A)
order = 2


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
ux1 = simplify(ux1.subs(P_1, P))
ux1 = simplify(series(ux1, kappa, n=order).removeO())
print("\n ux (1): ", ux1)

Pe = symbols("Pe", real=True)
#B0 = Pe*F0*tanh(gamma)/(gamma**3) + Pe*F0*tanh(gamma)/(gamma**3 * (Sc-1))*(sinh(xi*sqrt(i*omega))/(sqrt(i*omega)*sinh(sqrt(i*omega))) - sinh(gamma*xi)/(gamma*sinh(gamma)))
#B0 = simplify(series(B0, gamma, n=order).removeO())
#print("\n B0 (0): ", B0)


order = 3
B0_grad = Pe*F0*tanh(gamma)/(2*gamma*(rho*rho-gamma*gamma))*(sinh(xi*rho)/sinh(rho) - sinh(gamma*xi)/sinh(gamma))
B0_grad = simplify(series(B0_grad, gamma, n=order).removeO())
B0_grad = simplify(series(B0_grad, rho, n=order).removeO())
B0_grad = B0_grad.subs(gamma, sqrt(I*omega/Sc))
B0_grad = B0_grad.subs(rho,   sqrt(I*omega/D))
print("\n B0_grad (0): ", B0_grad)


print(simplify(expand(simplify(integrate(B0_grad*(re(B0_grad)-I*im(B0_grad)), (xi, -1, 1))/2))))


F0**2*Pe**2*omega*(Sc**2*(13650*I*omega - 135135) + Sc*(5932*omega**2 + 117390*I*omega - 585585) + 9*omega*(2833*omega + 28210*I))/(2554051500*Sc**2)


(25497*D*omega**2 + 13650*I*Sc**2*omega + 5932*Sc*omega**2)