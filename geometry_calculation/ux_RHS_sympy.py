from sympy import * 

#define variables
xi, eta = symbols("xi eta")
omega, gamma, kappa, kappa_p, Sc = symbols("omega gamma kappa kappa_p Sc")
A, B, C, D, P_1 = symbols("A B C D P_1")
i = symbols("i")

#define known quantities
u0  = A*(1-cosh(gamma*xi)/cosh(gamma))
ux1 = sin(kappa*eta)*(B*(cosh(kappa*xi)/cosh(kappa) - cosh(kappa_p *xi)/cosh(kappa_p)) + C*(cosh(kappa_p*xi)/cosh(kappa_p) - xi*sinh(gamma*xi)/sinh(gamma)))
uy1 = D*cos(kappa*eta)*(sinh(kappa_p*xi)/sinh(kappa_p) - sinh(kappa*xi)/sinh(kappa))
P1  = P_1*cosh(kappa*xi)*cos(kappa*eta)
P2_partic = (P_1*sin(2*kappa*eta)/3)*(kappa*xi*sinh(kappa*xi)/2 + cosh(kappa*xi)/3)

#define operators acting on known quantities
laplas1_u1x = -2*kappa*xi*cos(kappa*eta)*diff(ux1, xi, eta) + kappa*kappa*xi*sin(kappa*eta)*diff(ux1, xi) - 2*sin(kappa*eta)*diff(ux1, xi, xi)
grad1x_P1 = -kappa*xi*cos(kappa*eta)*diff(P1, xi)

#define all the terms
term1 = diff(u0, xi, xi) + diff(u0, eta, eta)
term2 = laplas1_u1x
term3 = grad1x_P1
term4 = diff(P2_partic, eta)

#add together, the substitute back and simplify
RHS = term1 + term2 + term3 + term4
RHS = RHS.subs(A, 1/(gamma*gamma))
RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
#RHS = RHS.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
#RHS = RHS.subs(gamma, sqrt(i*omega/Sc))

"""
print("added them together \n")
print("RHS: ", RHS)
expanded_RHS = expand(RHS)
print("\n \n Expanded RHS: ", expanded_RHS)
print("\n \n Simplifty RHS: ", simplify(RHS))
SIMP = simplify(expanded_RHS)
print("\n \n Simplifty expanded RHS: ", expand(SIMP))
print("\n \n Simplifty expanded RHS: ", latex(expand(SIMP)))
"""

#TEST IF RHS IS CORRECT
print("TEST IF THE RHS I HAVE USED IS CORRECT")
old_output = RHS 

my_RHS_term1 = -P_1*kappa*kappa*(kappa_p*kappa_p*xi*sinh(kappa*xi)/2 + kappa*cosh(kappa*xi))/(gamma*gamma)
my_RHS_term2 = P_1*sinh(kappa)*kappa_p*kappa_p/(sinh(kappa_p)*gamma*gamma)*(kappa*kappa*xi*sinh(kappa_p*xi)/2 + kappa_p*cosh(kappa_p*xi))
my_RHS_term3 = ((1+kappa*kappa*xi*xi/2)*cosh(gamma*xi)+ xi*sinh(gamma*xi)*(gamma*gamma+kappa*kappa/2)/gamma)/cosh(gamma)
my_RHS = my_RHS_term1+my_RHS_term2+my_RHS_term3

old_output = old_output.subs(eta, pi/(4*kappa))
old_output = old_output.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))

sol = old_output-my_RHS
sol = sol.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
sol = sol.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
sol = simplify(sol)
print("Eta independent terms add to: ", sol)





### CHECK IF SOLUTION IS CORRECT

print("\n\n\n THE FOLLOWING TEST SOLUTION OF H(xi), ALL SHOULD GIVE ZERO\n")
#kappa*xi term:
sol = (P_1*kappa*kappa/(2*gamma*gamma))* (( xi*sinh(kappa*xi)*(kappa**4 - gamma**4) - 4*gamma*gamma*kappa*cosh(kappa*xi))/((gamma*gamma-kappa*kappa)**2))
should_be_zero = ( diff(sol, xi, xi) - gamma*gamma*sol +my_RHS_term1)
should_be_zero = should_be_zero.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
print("first term gives ", simplify(should_be_zero))

#kappa_p term
sol = -P_1*sinh(kappa)*xi*sinh(kappa_p*xi)*kappa_p*kappa_p/(2*gamma*gamma*sinh(kappa_p))
should_be_zero = diff(sol, xi,xi) - gamma*gamma*sol + my_RHS_term2
should_be_zero = should_be_zero.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
should_be_zero = should_be_zero.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
print("second term gives ", simplify(expand(should_be_zero)))

#gamma*xi term
sol = (xi*sinh(gamma*xi)*(1+kappa*kappa*xi*xi/3) + gamma*xi*xi*cosh(gamma*xi))/(-4*gamma*cosh(gamma))
should_be_zero = diff(sol, xi,xi) - gamma*gamma*sol + my_RHS_term3
should_be_zero = should_be_zero.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
print("third term gives ", simplify(should_be_zero))


full_H = (P_1*kappa*kappa/(2*gamma*gamma))* (( xi*sinh(kappa*xi)*(kappa**4 - gamma**4) - 4*gamma*gamma*kappa*cosh(kappa*xi))/((gamma*gamma-kappa*kappa)**2)) -(xi*sinh(kappa_p*xi)/(2))*kappa_p*kappa_p*tanh(kappa)*tanh(gamma)/(gamma*cosh(kappa_p)*(kappa*tanh(kappa_p)-kappa_p*tanh(kappa))) + (xi*sinh(gamma*xi)*(1+kappa*kappa*xi*xi/3) + gamma*xi*xi*cosh(gamma*xi))/(-4*gamma*cosh(gamma))
answer = diff(full_H, xi, xi) - gamma*gamma*full_H + my_RHS
answer = answer.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
answer = answer.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
print("Test full H solution: ", simplify(answer))






