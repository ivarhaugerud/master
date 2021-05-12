from sympy import * 

#define variables
xi = symbols("xi", real=True)
kappa, F0, Pe, r = symbols("kappa F0 Pe r", real=True)
rho = symbols("rho")


#rho = sqrt(I*r)

B0 =  F0*Pe*cosh(rho*xi)/(2*rho*rho*rho*sinh(rho)) #- F0*Pe*xi*xi/(4*rho*rho) #
B0_deriv = diff(B0, xi)
B0_deriv_deriv = diff(B0_deriv, xi)




beta1 = 3/(2*rho**4) + (xi*xi - 1)/(4*rho*rho) + (xi*sinh(rho*xi) - cosh(rho*xi)/tanh(rho) - 2*cosh(rho*xi)/rho )/(2*rho*rho*sinh(rho))
beta1 *= Pe*F0 

beta0 = 3*xi*xi*(1-xi*xi/6)/(2*rho*rho) + xi*xi*(-1+xi*xi/3 - xi**4/15)/4 + cosh(rho*xi)*(24/(rho*rho*rho) + 6*xi*xi/rho - 2/rho)/(rho*rho*sinh(rho)) + xi*sinh(rho*xi)*(-18/(rho*rho) - xi*xi + 1)/(rho*rho*sinh(rho)) + cosh(rho*xi)*(2/rho + 1/tanh(rho))*(6/(rho*rho) + xi*xi - 1 - 4*xi*tanh(rho*xi)/rho)/(rho*rho*sinh(rho))
beta0 *= Pe*Pe*F0*F0*kappa/(8*rho*rho)
beta0 += cosh(kappa*xi)*(-1 + F0*F0*Pe*Pe/(60*rho**6)*(2*rho**4 - 30*rho**2 + 15*rho**2/(tanh(rho)**2) + 60*rho/tanh(rho) - 75) )/sinh(kappa)


beta2 = (1/4 - 2/(rho*rho) + 3/(rho**4) + xi*xi*(3/(rho*rho) - 1/2) + xi*xi*xi*xi/4) + (cosh(rho*xi)*(72/(rho*rho) + 6*xi*xi - 2)/rho + sinh(rho*xi)*(30*xi/(rho*rho) + xi**3 - xi) )/sinh(rho) + (2/rho + 1/tanh(rho))*(cosh(rho*xi)*(1-xi*xi - 10/(rho*rho)) - 4*xi*sinh(rho*xi)/rho)/sinh(rho)
beta2 *= Pe*Pe*F0*F0*kappa/(8*rho**4)
beta2 -= (F0*F0*Pe*Pe*kappa/(4*sqrt(2)*sinh(sqrt(2)*rho)*rho**7))*(40 + 8*rho/tanh(rho) - 3*rho*rho/(sinh(rho)**2))*cosh(sqrt(2)*rho*xi)


B2 = 5*xi*xi - 1 -3*xi*sinh(rho*xi)/sinh(rho) - 2*cosh(rho)*rho*xi*sinh(rho*xi)/(sinh(rho)*sinh(rho)) + rho*xi*xi*cosh(rho*xi)/sinh(rho) - xi*sinh(rho*xi)/sinh(rho)
B2 *= F0*Pe/(8*rho*rho)
B2 -= F0*Pe*cosh(rho*xi)*(rho*rho + 6 - 2*rho*rho/(tanh(rho)**2) -4*rho/tanh(rho) )/(8*rho*rho*rho*sinh(rho))
"""

#term1 = kappa* ( beta0 + conjugate(beta0) ) #might be kappa^2
term1 = 0
term2 = 0.5*diff(beta1, xi)*diff(conjugate(beta1), xi)
term3 = conjugate(B0_deriv) * ( simplify(2*diff(beta2, xi) - diff(beta1, xi) - kappa*xi*xi*beta1))
term3 = (simplify(term3))
#term3 += 0#conjugate(term3)
term4 = 0.5*B0_deriv*conjugate(B0_deriv)
"""


#term3_2 = simplify(-(re(B0_deriv)-im(B0_deriv)) * diff(beta1, xi))
#print("\n\n", (integrate(term3_2, (xi, -1, 1))))
#term3_2_int =  I*F0**2*Pe**2*(Integral(xi*sqrt(I*r)*sinh(xi*sqrt(I*r)), (xi, -1, 1)) + Integral(xi*sinh(xi*sqrt(I*r))*tanh(sqrt(I*r)), (xi, -1, 1)) + Integral(-xi**2*sinh(sqrt(I*r))*tanh(sqrt(I*r)), (xi, -1, 1)) + Integral(-sqrt(I*r)*sinh(xi*sqrt(I*r))*re(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)) + Integral(-sqrt(I*r)*sinh(xi*sqrt(I*r))*im(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)) + Integral(-sinh(xi*sqrt(I*r))*tanh(sqrt(I*r))*re(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)) + Integral(-sinh(xi*sqrt(I*r))*tanh(sqrt(I*r))*im(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)) + Integral(xi*sinh(sqrt(I*r))*tanh(sqrt(I*r))*re(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)) + Integral(xi*sinh(sqrt(I*r))*tanh(sqrt(I*r))*im(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)) + Integral(-xi**2*sqrt(I*r)*cosh(xi*sqrt(I*r))*tanh(sqrt(I*r)), (xi, -1, 1)) + Integral(xi*sqrt(I*r)*cosh(xi*sqrt(I*r))*tanh(sqrt(I*r))*re(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)) + Integral(xi*sqrt(I*r)*cosh(xi*sqrt(I*r))*tanh(sqrt(I*r))*im(sinh(xi*sqrt(I*r))/sinh(sqrt(I*r))), (xi, -1, 1)))/(4*r**2*sinh(sqrt(I*r))*tanh(sqrt(I*r)))
"""
beta1 = (xi*sinh(rho*xi) - cosh(rho*xi)/tanh(rho) - 2*cosh(rho*xi)/rho )/(2*rho*rho*sinh(rho))
beta1 *= Pe*F0 
term3_3 = simplify( (conjugate(B0_deriv) )* ( (- kappa*xi*beta1)))
print("\n\n", (integrate(term3_3, (xi, -1, 1))))

beta1 = 3/(2*rho**4) + (xi*xi - 1)/(4*rho*rho)
beta1 *= Pe*F0 
term3_3 = simplify( (conjugate(B0_deriv) )* ( (- kappa*xi*beta1)))
print("\n\n", (integrate(term3_3, (xi, -1, 1))))
"""
#term3_3_integrated = F0**2*Pe**2*kappa*(-sinh(rho)/tanh(rho) + cosh(rho) - 5*sinh(rho)/rho + 2*cosh(rho)/(rho*tanh(rho)) - 2*sinh(rho)/(rho**2*tanh(rho)) + 10*cosh(rho)/rho**2 - 10*sinh(rho)/rho**3)/(4*rho**3*sinh(rho)*conjugate(rho)**2) - F0**2*Pe**2*kappa*(sinh(rho)/tanh(rho) - cosh(rho) + 5*sinh(rho)/rho - 2*cosh(rho)/(rho*tanh(rho)) + 2*sinh(rho)/(rho**2*tanh(rho)) - 10*cosh(rho)/rho**2 + 10*sinh(rho)/rho**3)/(4*rho**3*sinh(rho)*conjugate(rho)**2)
#term3_3_integrated +=  -F0**2*Pe**2*kappa/(30*rho**2*conjugate(rho)**2) + F0**2*Pe**2*kappa/(2*rho**4*conjugate(rho)**2)
#term3_3_integrated += -F0**2*Pe**2*kappa*(-2*rho**2*sinh(conjugate(rho))/conjugate(rho)**2 - 6*rho**2*sinh(conjugate(rho))/conjugate(rho)**4 + 6*rho**2*cosh(conjugate(rho))/conjugate(rho)**3 - 6*sinh(conjugate(rho))/conjugate(rho)**2 + 6*cosh(conjugate(rho))/conjugate(rho))/(8*rho**4*sinh(conjugate(rho))*conjugate(rho)**2) + F0**2*Pe**2*kappa*(2*rho**2*sinh(conjugate(rho))/conjugate(rho)**2 + 6*rho**2*sinh(conjugate(rho))/conjugate(rho)**4 - 6*rho**2*cosh(conjugate(rho))/conjugate(rho)**3 + 6*sinh(conjugate(rho))/conjugate(rho)**2 - 6*cosh(conjugate(rho))/conjugate(rho))/(8*rho**4*sinh(conjugate(rho))*conjugate(rho)**2)

term3_1 = simplify(conjugate(B0_deriv) * ( (2*diff(beta2, xi))))

print("\n\n", (integrate(term3_1, (xi, -1, 1))))
term3_1_integrated = 

"""
#print("\n\n\n term 2 : ", term2)
#print("\n\n\n term 3: ", term3 + conjugate(term3))
print("\n\n\n term 4: ", term4)

#full = term1 + term2 + term3 + term4 
#print(full)
print("\n\n\n\n")
#full = 
#full = (F0**2*Pe**2*(2.0*rho**4*(xi*sinh(rho) - sinh(rho*xi))*(xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*sinh(rho)*sinh(sqrt(2)*rho)*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*tanh(rho)*tanh(conjugate(rho))*conjugate(rho)**4 + rho**4*(xi*sinh(rho) - sinh(rho*xi))*(2*F0*Pe*kappa*(40*sinh(conjugate(rho))**2*tanh(conjugate(rho)) + 8*sinh(conjugate(rho))**2*conjugate(rho) - 3*tanh(conjugate(rho))*conjugate(rho)**2)*sinh(sqrt(2)*xi*conjugate(rho)) - F0*Pe*kappa*(xi**3*sinh(conjugate(rho))*tanh(conjugate(rho))*conjugate(rho)**2 - 2*xi*(0.5*conjugate(rho)**2 - 3)*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*((xi**2 - 1)*conjugate(rho)**2 + 42)*cosh(xi*conjugate(rho))*conjugate(rho) + 3*((3*xi**2 - 1)*conjugate(rho)**2 + 34)*sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - (2*tanh(conjugate(rho)) + conjugate(rho))*(6*xi*cosh(xi*conjugate(rho))*conjugate(rho) + ((xi**2 - 1)*conjugate(rho)**2 + 10)*sinh(xi*conjugate(rho)) + 4*sinh(xi*conjugate(rho))))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho)) + kappa*xi**2*((xi**2 - 1)*sinh(conjugate(rho))*tanh(conjugate(rho))*conjugate(rho)**2 - 2*(-xi*sinh(xi*conjugate(rho))*tanh(conjugate(rho))*conjugate(rho) + 2*cosh(xi*conjugate(rho))*tanh(conjugate(rho)) + cosh(xi*conjugate(rho))*conjugate(rho))*conjugate(rho) + 6*sinh(conjugate(rho))*tanh(conjugate(rho)))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*conjugate(rho)**2 + 2*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*conjugate(rho)**4)*sinh(rho)*sinh(sqrt(2)*rho)*tanh(rho) + rho**4*(-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(rho)*sinh(sqrt(2)*rho)*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*conjugate(rho)**4 + (xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*(2*F0*Pe*kappa*(-3*rho**2*tanh(rho) + 8*rho*sinh(rho)**2 + 40*sinh(rho)**2*tanh(rho))*sinh(sqrt(2)*rho*xi) - F0*Pe*kappa*(rho**2*xi**3*sinh(rho)*tanh(rho) - 2*xi*(0.5*rho**2 - 3)*sinh(rho)*tanh(rho) - (rho + 2*tanh(rho))*(6*rho*xi*cosh(rho*xi) + (rho**2*(xi**2 - 1) + 10)*sinh(rho*xi) + 4*sinh(rho*xi)) + (rho*xi*(rho**2*(xi**2 - 1) + 42)*cosh(rho*xi) + 3*(rho**2*(3*xi**2 - 1) + 34)*sinh(rho*xi))*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho) + kappa*rho**2*xi**2*(rho**2*(xi**2 - 1)*sinh(rho)*tanh(rho) - 2*rho*(-rho*xi*sinh(rho*xi)*tanh(rho) + rho*cosh(rho*xi) + 2*cosh(rho*xi)*tanh(rho)) + 6*sinh(rho)*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho) + 2*rho**4*(-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*tanh(conjugate(rho))*conjugate(rho)**4)/(8*rho**6*sinh(rho)**2*sinh(sqrt(2)*rho)*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))**2*tanh(rho)*tanh(conjugate(rho))*conjugate(rho)**6))
answ = integrate(term3, (xi, -1, 1))
print(answ)
print("\n\n\n\n")
print(simplify(answ))
#print(term3)
#print("second term")
#print(simplify(term2 + term3 + conjugate(term3)))

term2_integrated = F0**2*Pe**2*( (rho**8*conjugate(rho)**2 )/3 + 0.75*rho**8 - 0.5*rho**8*conjugate(rho)/tanh(conjugate(rho)) - 0.25*rho**8*conjugate(rho)**2/tanh(conjugate(rho))**2 - 1.5*rho**6*conjugate(rho)**2 + 1.5*rho**6*conjugate(rho)**3/tanh(conjugate(rho)) - 0.25*rho**6*conjugate(rho)**4/tanh(rho)**2 - 1.0*rho**5*conjugate(rho)**4/tanh(rho) + rho**4*conjugate(rho)**5/tanh(conjugate(rho)) + 0.25*rho**4*conjugate(rho)**6/tanh(conjugate(rho))**2 - 1.5*rho**3*conjugate(rho)**6/tanh(rho) - 0.333333333333333*rho**2*conjugate(rho)**8 + 1.5*rho**2*conjugate(rho)**6 + 0.25*rho**2*conjugate(rho)**8/tanh(rho)**2 + 0.5*rho*conjugate(rho)**8/tanh(rho) - 0.75*conjugate(rho)**8)/(rho**4*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*conjugate(rho)**4)
term4_integrated = F0**2*Pe**2*(0.0833333333333333*rho**4*conjugate(rho)**2 + 0.25*rho**4 - 0.25*rho**4*conjugate(rho)/tanh(conjugate(rho)) - 0.0833333333333333*rho**2*conjugate(rho)**4 + 0.25*rho*conjugate(rho)**4/tanh(rho) - 0.25*conjugate(rho)**4)/(rho**4*(rho**2 - conjugate(rho)**2)*conjugate(rho)**4)
"""

#print(simplify(term2_integrated + term4_integrated))
