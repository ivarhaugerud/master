from sympy import * 

#define variables
xi, eta, t = symbols("xi eta t")
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

#insert complex values 

#zeroth order
"""
u0 = u0.subs(gamma, sqrt(I*omega/Sc))
u0 = u0*exp(I*omega*t)
u0 = (u0 + u0.subs(I, -I))/2

zeroth_order_int = simplify(integrate(u0*u0, (xi, -1, 1), conds="none"))
print("0'th: ", zeroth_order_int)
#zeroth_order_int = integrate(zeroth_order_int.subs(I, sqrt(-1)), (t, 0, 2*pi/omega), conds="none")

#first order
ux1 = ux1.subs(kappa_p, sqrt(kappa*kappa+I*omega/Sc))
ux1 = ux1.subs(gamma, sqrt(I*omega/Sc))
ux1 = ux1*exp(I*omega*t)
"""
uy1 = uy1.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
uy1 = uy1.subs(P_1, (gamma*F0*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
uy1 = uy1.subs(kappa_p, sqrt(kappa*kappa+I*omega/Sc))
uy1 = uy1*exp(I*omega*t)

#ux1 = (ux1 + ux1.subs(I, -I))/2
uy1 = (uy1 + uy1.subs(I, -I))/2
print("added together")
#u1x_int = simplify(integrate(ux1*ux1, (xi, -1, 1), conds="none"))
u1y_int = simplify(integrate(uy1*uy1, (xi, -1, 1), conds="none"))
print("1'th: ", (u1y_int))


#second order
"""
sol_x_no_eta = sol_x_no_eta.subs(P_1,      (gamma*F0*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
sol_x_no_eta = sol_x_no_eta.subs(kappa_pp, sqrt(4*kappa*kappa+I*omega/Sc))
sol_x_no_eta = sol_x_no_eta.subs(kappa_p,  sqrt(kappa*kappa+I*omega/Sc))
sol_x_no_eta = sol_x_no_eta.subs(gamma,    sqrt(I*omega/Sc))

sol_x_no_eta = re(sol_x_no_eta*exp(I*omega*t))
"""