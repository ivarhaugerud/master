from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, kappa, kappa_p, kappa_pp, Sc = symbols("omega gamma kappa kappa_p kappa_pp Sc")
B, C, D, P_1, F0, i, P_2 = symbols("B C D P_1 F0 i P_2")
gamma_p = symbols("gamma_p")
from sympy import I 

#define known quantities
u0  = F0*(1-cosh(gamma*xi)/cosh(gamma))/(gamma*gamma)
ux1 = sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)*cos(kappa*eta)

sol_x_no_eta = P_1*sinh(kappa)*(kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa) - kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p) + gamma*gamma*cosh(gamma*xi)/cosh(gamma))/(2*gamma*gamma) + cosh(gamma*xi)*F0*(1-xi*xi)/(4*cosh(gamma))
sol_x_eta = P_1*sinh(kappa)*(kappa_p*kappa_p*xi*sinh(kappa_p*xi)/sinh(kappa_p) -kappa*kappa*xi*sinh(kappa*xi)/sinh(kappa) - gamma*gamma*cosh(kappa_pp*xi)/cosh(kappa_pp))/(2*gamma*gamma) + F0*(xi*xi*cosh(gamma*xi)/cosh(gamma) - cosh(kappa_pp*xi)/cosh(kappa_pp))/4 - P_2*cosh(2*kappa)*(cosh(2*kappa*xi)/cosh(2*kappa) - cosh(kappa_pp*xi)/cosh(kappa_pp))
sol_y = (P_1*kappa*sinh(kappa)/(2*gamma*gamma*tanh(kappa_p)))*(kappa_p*xi*cosh(kappa_p*xi)/cosh(kappa_p) - kappa_p*sinh(kappa_pp*xi)/sinh(kappa_pp) -kappa*(tanh(kappa_p)/tanh(kappa))*(xi*cosh(kappa*xi)/cosh(kappa) - sinh(kappa_pp*xi)/sinh(kappa_pp))) - P_2*sinh(2*kappa)*(sinh(2*kappa*xi)/sinh(2*kappa) - sinh(kappa_pp*xi)/sinh(kappa_pp))

u0_prime = u0.subs(gamma, gamma_p)

test = integrate(u0*u0_prime, (xi, -1, 1))
print(simplify(test))