from sympy import * 
order = 2

#define variables
xi = symbols("xi")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
P_1, F0 = symbols("P_1, F0")

#define known quantities
u0  = F0*(1-cosh(gamma*xi)/cosh(gamma))/(gamma*gamma)
ux1 = ((P_1*kappa*cosh(kappa)/(gamma*gamma))*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + (F0*tanh(gamma)/gamma)*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = (kappa*P_1*sinh(kappa)/(gamma*gamma))*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)

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

i, Pe = symbols("i, Pe")
B0 = Pe*F0*tanh(gamma)/(gamma**3) + Pe*F0*tanh(gamma)/(gamma**3 * (Sc-1))*(sinh(xi*sqrt(i*omega))/(sqrt(i*omega)*sinh(sqrt(i*omega))) - sinh(gamma*xi)/(gamma*sinh(gamma)))
B0 = simplify(series(B0, gamma, n=order).removeO())
print("\n B0 (0): ", B0)


B0_grad = Pe*F0*tanh(gamma)/(gamma**3 * (Sc-1))*(sinh(xi*sqrt(i*omega))/sinh(sqrt(i*omega)) - sinh(gamma*xi)/sinh(gamma))
B0_grad = simplify(series(B0_grad, gamma, n=order).removeO())
print("\n B0_grad (0): ", B0_grad)

u0 = simplify(series(u0, gamma, n=order).removeO())
print("\n u (0): ", u0)

