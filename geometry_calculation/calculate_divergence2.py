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

#calculate RHS, and insert known values
RHS = -sin(kappa*eta)*cos(kappa*eta)*kappa*xi*diff(u0, xi) + kappa*xi*cos(kappa*eta)*diff(ux1, xi) + sin(kappa*eta)*diff(uy1, xi)
RHS = RHS.subs(A, 1/(gamma*gamma))
RHS = RHS.subs(B, P_1*kappa*cosh(kappa)/(gamma*gamma))
RHS = RHS.subs(C, tanh(gamma)/gamma)
RHS = RHS.subs(D, kappa*P_1*sinh(kappa)/(gamma*gamma))
RHS = RHS.subs(kappa_p, sqrt(kappa*kappa+gamma*gamma))
RHS = simplify(RHS/(sin(2*kappa*eta)))
print("\n\n RHS: ", RHS)


#define second order velocity field with homogeneous pressure part, and without homogeneous part
Fx_1 = -P_1*kappa*((2*kappa*kappa + gamma*gamma/6 + 9*kappa*kappa*kappa*kappa/(2*gamma*gamma))*kappa*xi*sinh(kappa*xi) - 2*(2*kappa*kappa + gamma*gamma/3)*cosh(kappa*xi)/3)/((gamma*gamma+3*kappa*kappa)**2)
Fx_2 = P_1*(gamma*gamma+kappa*kappa)*sinh(kappa)*xi*sinh(sqrt(gamma*gamma+kappa*kappa)*xi)/(2*gamma*gamma*sinh(sqrt(kappa*kappa+gamma*gamma)))
Fx_3 = 3*(cosh(gamma*xi)*(2*xi*xi*kappa**4 - 2*kappa*kappa/3 + gamma*gamma/3)/(8*kappa**4) + gamma*xi*sinh(gamma*xi)*(kappa*kappa/(gamma*gamma)+1/3)/(4*kappa*kappa))/(2*cosh(gamma))
Fx = Fx_1+Fx_2+Fx_3

#Fy_1 = 2*sinh(kappa*xi)*kappa*kappa*gamma*gamma*(1-gamma*gamma/(3*kappa*kappa))/(3*(gamma*gamma+3*kappa*kappa)**2)
#Fy_2 = kappa*xi*cosh(kappa*xi)*( (gamma**4)/6 - 9*(kappa**4)/2 - gamma*gamma*kappa*kappa )/((gamma*gamma+3*kappa*kappa))
#Fy_3 = sinh(kappa)*kappa_p*xi*cosh(kappa_p*xi)/(2*sinh(kappa_p))
Fy = 2*sinh(kappa*xi)*kappa*kappa*gamma*gamma*(1-gamma*gamma/(3*kappa*kappa))/(3*(gamma*gamma+3*kappa*kappa)**2) + kappa*xi*cosh(kappa*xi)*(gamma*gamma*gamma*gamma/6 - 9*kappa*kappa*kappa*kappa/2 - gamma*gamma*kappa*kappa)/((gamma*gamma+3*kappa*kappa)**2) + sinh(kappa)*kappa_p*xi*cosh(kappa_p*xi)/(2*sinh(kappa_p))
Fy = P_1*kappa*Fy/(gamma*gamma)
Fy = Fy.subs(kappa_p, sqrt(gamma*gamma+kappa*kappa))

LHS = simplify((diff(Fy, xi) - 2*kappa*Fx))
print("\n\n LHS: ", LHS)

difference = RHS-LHS
difference = difference.subs(P_1, (gamma*tanh(gamma)/(kappa*cosh(kappa)))/(1-kappa_p*tanh(kappa)/(kappa*tanh(kappa_p))))
difference = difference.subs(kappa_p, sqrt(gamma*gamma+kappa*kappa))
#difference = difference.subs(gamma, i*omega/Sc)
#difference = difference.subs(i, sqrt(-1))

print("\n\n Difference: ", simplify(expand(difference)))

1.0*(0.125*gamma**3*cosh(gamma*xi) + 0.25*gamma**2*kappa**2*xi*sinh(gamma*xi) + 0.25*gamma*kappa**4*xi**2*cosh(gamma*xi) - 0.25*gamma*kappa**2*cosh(gamma*xi) + 0.75*kappa**4*xi*sinh(gamma*xi))/(gamma*kappa**3*cosh(gamma))

 #Difference:  -(-0.125*gamma**3*kappa*cosh(kappa)*cosh(gamma*xi)*tanh(kappa_p) + 0.125*gamma**3*kappa_p*sinh(kappa)*cosh(gamma*xi) + 1.5*gamma**2*kappa**4*xi*sinh(gamma)*sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))*tanh(kappa_p)/sinh(sqrt(gamma**2 + kappa**2)) - 0.25*gamma**2*kappa**3*xi*sinh(gamma*xi)*cosh(kappa)*tanh(kappa_p) + 0.25*gamma**2*kappa**2*kappa_p*xi*sinh(kappa)*sinh(gamma*xi) - 0.25*gamma*kappa**5*xi**2*cosh(kappa)*cosh(gamma*xi)*tanh(kappa_p) + 0.25*gamma*kappa**4*kappa_p*xi**2*sinh(kappa)*cosh(gamma*xi) + 0.25*gamma*kappa**3*cosh(kappa)*cosh(gamma*xi)*tanh(kappa_p) - 0.25*gamma*kappa**2*kappa_p*sinh(kappa)*cosh(gamma*xi) + 1.5*kappa**6*xi*sinh(gamma)*sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))*tanh(kappa_p)/sinh(sqrt(gamma**2 + kappa**2)) - 0.75*kappa**5*xi*sinh(gamma*xi)*cosh(kappa)*tanh(kappa_p) + 0.5*kappa**4*kappa_p*xi*sqrt(gamma**2 + kappa**2)*sinh(gamma)*sinh(kappa)*sinh(xi*sqrt(gamma**2 + kappa**2))/cosh(sqrt(gamma**2 + kappa**2)) + 0.75*kappa**4*kappa_p*xi*sinh(kappa)*sinh(gamma*xi))/(gamma*kappa**3*(kappa*cosh(kappa)*tanh(kappa_p) - kappa_p*sinh(kappa))*cosh(gamma))
