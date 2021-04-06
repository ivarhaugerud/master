from numpy import *
import matplotlib.pyplot as plt 
import scipy.integrate as scp 
F0 = 6
tau = 3.0
D = 0.1
Pe = 1/D
omega = 2*pi/tau
xi = linspace(-1, 1, int(1e4))



kappas = linspace(0.01, 0.4, 10)
omega = logspace(-3, 3, 60)
D0 = zeros(len(omega))

D_eff = zeros((len(kappas), len(omega)))

for j in range(len(kappas)):
	for i in range(len(omega)):
		rho = sqrt(1j*omega[i]/D)
		kappa = kappas[j]
		gamma = sqrt(1j*omega[i]/1e4)
		D0[i] = 1 + Pe*Pe*F0*F0*tanh(gamma)*tanh(conjugate(gamma))/(4*gamma*conjugate(gamma)*(gamma**4 - rho**4))*(1/(gamma*gamma)*(gamma/tanh(gamma) - conjugate(gamma/tanh(gamma))) - 1/(rho*rho)*(rho/tanh(rho) - conjugate(rho/tanh(rho))))

		"""

		B0_full   = F0*Pe/(rho*rho)*(cosh(rho*xi)/(rho*sinh(rho)) - cosh(gamma*xi)/(gamma*sinh(gamma)))
		B0        = F0*Pe*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) - F0*Pe*xi*xi/(2*rho*rho)
		B0_approx = F0*Pe*cosh(rho*xi)/(rho*rho*rho*sinh(rho)) - F0*Pe*xi*xi/(2*rho*rho)

		beta1  = 3/(2*rho**4) + (xi*xi - 1)/(4*rho*rho) + (xi*sinh(rho*xi) - cosh(rho*xi)/tanh(rho) - 2*cosh(rho*xi)/rho )/(2*rho*rho*sinh(rho))
		beta1 *= Pe*F0 


		beta0 = 3*xi*xi*(1-xi*xi/6)/(2*rho*rho) + xi*xi*(-1+xi*xi/3 - xi**4/15)/4 + cosh(rho*xi)*(24/(rho*rho*rho) + 6*xi*xi/rho - 2/rho)/(rho*rho*sinh(rho)) + xi*sinh(rho*xi)*(-18/(rho*rho) - xi*xi + 1)/(rho*rho*sinh(rho)) + cosh(rho*xi)*(2/rho + 1/tanh(rho))*(6/(rho*rho) + xi*xi - 1 - 4*xi*tanh(rho*xi)/rho)/(rho*rho*sinh(rho))
		beta0 *= Pe*Pe*F0*F0*kappa/(8*rho*rho)
		beta0 += (-1 + F0*F0*Pe*Pe/(60*rho**6)*(2*rho**4 - 20*rho**2 + 15*rho**2/(tanh(rho)**2) + 60*rho/tanh(rho) - 75) )*cosh(kappa*xi)/sinh(kappa)#*(1/kappa + kappa*(xi*xi-1/3)/2)#*cosh(kappa*xi)/sinh(kappa)
		
		

		beta2 = (1/4 - 2/(rho*rho) + 3/(rho**4) + xi*xi*(2/(rho*rho) - 1/2) + xi*xi*xi*xi/4) + (cosh(rho*xi)*(72/(rho*rho) + 6*xi*xi - 2)/rho + sinh(rho*xi)*(30*xi/(rho*rho) + xi**3 - xi) )/sinh(rho) + (2/rho + 1/tanh(rho))*(cosh(rho*xi)*(1-xi*xi - 10/(rho*rho)) - 4*xi*sinh(rho*xi)/rho)/sinh(rho)
		beta2 *= Pe*Pe*F0*F0*kappa/(8*rho**4)
		beta2 -= (F0*F0*Pe*Pe*kappa/(4*sqrt(2)*sinh(sqrt(2)*rho)*rho**7))*(40 + 8*rho/tanh(rho) - 3*rho*rho/(sinh(rho)**2))*cosh(sqrt(2)*rho*xi)
		
		B2 = 5*xi*xi - 1 -3*xi*sinh(rho*xi)/sinh(rho) - 2*cosh(rho)*rho*xi*sinh(rho*xi)/(sinh(rho)*sinh(rho)) + rho*xi*xi*cosh(rho*xi)/sinh(rho) - xi*sinh(rho*xi)/sinh(rho)
		B2 *= F0*Pe/(8*rho*rho)
		B2 -= F0*Pe*cosh(rho*xi)*(rho*rho + 6 - 2*rho*rho/(tanh(rho)**2) -4*rho/tanh(rho) )/(8*rho*rho*rho*sinh(rho))

		plt.title(str(omega[i]))
		plt.plot(xi, gradient(B0_approx, xi))
		plt.plot(xi, gradient(beta0, xi), "k")
		plt.plot(xi, gradient(beta1, xi), "r")
		plt.plot(xi, gradient(beta2, xi), "g")
		plt.plot(xi, gradient(B2, xi), "b")
		plt.legend(loc="best")
		plt.show()
		"""
		#D_para_num_int = F0**2*Pe**2*(rho**4*(xi*sinh(rho) - sinh(rho*xi))*(2*F0*Pe*kappa*(40*sinh(conjugate(rho))**2*tanh(conjugate(rho)) + 8*sinh(conjugate(rho))**2*conjugate(rho) - 3*tanh(conjugate(rho))*conjugate(rho)**2)*sinh(sqrt(2)*xi*conjugate(rho)) - F0*Pe*kappa*(xi**3*sinh(conjugate(rho))*tanh(conjugate(rho))*conjugate(rho)**2 - 2*xi*(0.5*conjugate(rho)**2 - 3)*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*((xi**2 - 1)*conjugate(rho)**2 + 42)*cosh(xi*conjugate(rho))*conjugate(rho) + 3*((3*xi**2 - 1)*conjugate(rho)**2 + 34)*sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - (2*tanh(conjugate(rho)) + conjugate(rho))*(6*xi*cosh(xi*conjugate(rho))*conjugate(rho) + ((xi**2 - 1)*conjugate(rho)**2 + 10)*sinh(xi*conjugate(rho)) + 4*sinh(xi*conjugate(rho))))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho)) + kappa*xi**2*((xi**2 - 1)*sinh(conjugate(rho))*tanh(conjugate(rho))*conjugate(rho)**2 - 2*(-xi*sinh(xi*conjugate(rho))*tanh(conjugate(rho))*conjugate(rho) + 2*cosh(xi*conjugate(rho))*tanh(conjugate(rho)) + cosh(xi*conjugate(rho))*conjugate(rho))*conjugate(rho) + 6*sinh(conjugate(rho))*tanh(conjugate(rho)))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*conjugate(rho)**2 + 2*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*conjugate(rho)**4)*sinh(rho)*sinh(sqrt(2)*rho)*tanh(rho) + rho**4*(-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*(xi*sinh(conjugate(rho))*tanh(conjugate(rho)) + (xi*cosh(xi*conjugate(rho))*conjugate(rho) - sinh(xi*conjugate(rho)))*tanh(conjugate(rho)) - sinh(xi*conjugate(rho))*conjugate(rho))*sinh(rho)*sinh(sqrt(2)*rho)*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*conjugate(rho)**4 + (xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*(2*F0*Pe*kappa*(-3*rho**2*tanh(rho) + 8*rho*sinh(rho)**2 + 40*sinh(rho)**2*tanh(rho))*sinh(sqrt(2)*rho*xi) - F0*Pe*kappa*(rho**2*xi**3*sinh(rho)*tanh(rho) - 2*xi*(0.5*rho**2 - 3)*sinh(rho)*tanh(rho) - (rho + 2*tanh(rho))*(6*rho*xi*cosh(rho*xi) + (rho**2*(xi**2 - 1) + 10)*sinh(rho*xi) + 4*sinh(rho*xi)) + (rho*xi*(rho**2*(xi**2 - 1) + 42)*cosh(rho*xi) + 3*(rho**2*(3*xi**2 - 1) + 34)*sinh(rho*xi))*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho) + kappa*rho**2*xi**2*(rho**2*(xi**2 - 1)*sinh(rho)*tanh(rho) - 2*rho*(-rho*xi*sinh(rho*xi)*tanh(rho) + rho*cosh(rho*xi) + 2*cosh(rho*xi)*tanh(rho)) + 6*sinh(rho)*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho) + 2*rho**4*(-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho))*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))*tanh(conjugate(rho))*conjugate(rho)**4)/(8*rho**6*sinh(rho)**2*sinh(sqrt(2)*rho)*sinh(sqrt(2)*conjugate(rho))*sinh(conjugate(rho))**2*tanh(rho)*tanh(conjugate(rho))*conjugate(rho)**6)
		#D_para_num_int += F0**2*Pe**2*(xi*sinh(rho) - sinh(rho*xi))*(xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))/(8*rho**2*sinh(rho)*sinh(conjugate(rho))*conjugate(rho)**2)

		#D_para_num_int = (-0.25*F0*Pe*xi/rho**2 + 0.25*F0*Pe*sinh(rho*xi)/(rho**2*sinh(rho)))*(-F0*Pe*xi/(2*conjugate(rho)**2) + F0*Pe*sinh(xi*conjugate(rho))/(2*sinh(conjugate(rho))*conjugate(rho)**2))
		#D_para_num_int += (-0.25*F0*Pe*xi/rho**2 + 0.25*F0*Pe*sinh(rho*xi)/(rho**2*sinh(rho)))*(-F0*Pe*xi/(2*conjugate(rho)**2) + F0*Pe*sinh(xi*conjugate(rho))/(2*sinh(conjugate(rho))*conjugate(rho)**2))
		
		D_para = F0**2*Pe**2*(xi*sinh(conjugate(rho)) - sinh(xi*conjugate(rho)))*(2*F0*Pe*kappa*(-3*rho**2*tanh(rho) + 8*rho*sinh(rho)**2 + 40*sinh(rho)**2*tanh(rho))*sinh(sqrt(2)*rho*xi) - F0*Pe*kappa*(rho**2*xi**3*sinh(rho)*tanh(rho) - 2*xi*(0.5*rho**2 - 3)*sinh(rho)*tanh(rho) - (rho + 2*tanh(rho))*(6*rho*xi*cosh(rho*xi) + (rho**2*(xi**2 - 1) + 10)*sinh(rho*xi) + 4*sinh(rho*xi)) + (rho*xi*(rho**2*(xi**2 - 1) + 42)*cosh(rho*xi) + 3*(rho**2*(3*xi**2 - 1) + 34)*sinh(rho*xi))*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho) + kappa*rho**2*xi**2*(rho**2*(xi**2 - 1)*sinh(rho)*tanh(rho) - 2*rho*(-rho*xi*sinh(rho*xi)*tanh(rho) + rho*cosh(rho*xi) + 2*cosh(rho*xi)*tanh(rho)) + 6*sinh(rho)*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho) + 2*rho**4*(-rho*sinh(rho*xi) + xi*sinh(rho)*tanh(rho) + (rho*xi*cosh(rho*xi) - sinh(rho*xi))*tanh(rho))*sinh(rho)*sinh(sqrt(2)*rho))/(8*rho**6*sinh(rho)**2*sinh(sqrt(2)*rho)*sinh(conjugate(rho))*tanh(rho)*conjugate(rho)**2)
		#D_para_num_int += D_para
		analytic_integration = -F0**2*Pe**2*((rho**2 - conjugate(rho)**2)*(0.75*rho**4*(rho**4 - conjugate(rho)**4)*tanh(rho)**2*conjugate(rho)**2 - rho**4*(-1.5*rho**4 + 4.5*rho**2*conjugate(rho)**2 + 3*conjugate(rho)**4)*tanh(rho)**2*tanh(conjugate(rho))*conjugate(rho) + 0.75*rho**2*(rho**4 - conjugate(rho)**4)*tanh(conjugate(rho))**2*conjugate(rho)**4 + rho*(3.0*rho**4 + 4.5*rho**2*conjugate(rho)**2 - 1.5*conjugate(rho)**4)*tanh(rho)*tanh(conjugate(rho))**2*conjugate(rho)**4 - (rho**8*conjugate(rho)**2 + 2.25*rho**8 - 4.5*rho**6*conjugate(rho)**2 - 0.999999999999999*rho**2*conjugate(rho)**8 + 4.5*rho**2*conjugate(rho)**6 - 2.25*conjugate(rho)**8)*tanh(rho)**2*tanh(conjugate(rho))**2) - 3*(-0.25*rho**4*tanh(rho)*conjugate(rho) + 0.25*rho*tanh(conjugate(rho))*conjugate(rho)**4 + (0.0833333333333333*rho**4*conjugate(rho)**2 + 0.25*rho**4 - 0.0833333333333333*rho**2*conjugate(rho)**4 - 0.25*conjugate(rho)**4)*tanh(rho)*tanh(conjugate(rho)))*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*tanh(rho)*tanh(conjugate(rho)))/(3*rho**4*(rho**2 - conjugate(rho)**2)*(rho**6 - 3*rho**4*conjugate(rho)**2 + 3*rho**2*conjugate(rho)**4 - conjugate(rho)**6)*tanh(rho)**2*tanh(conjugate(rho))**2*conjugate(rho)**4)
		D_eff[j, i] = scp.trapz(D_para, xi) + analytic_integration
		#print(real(scp.trapz(D_para_num_int, xi) - D_eff[j, i]))
		#plt.plot(xi, D_para)
		#plt.show()

a_o = logspace(-2.5, 2.5, 10)
a = load("../finite_element/data/D_parallels_kappa_D1.npy")

plt.plot(a_o, a, "o")
for i in range(len(kappas)):
	plt.plot(omega, D_eff[i, :])
plt.xscale("log")
plt.show()

epsilon = 0.3
for i in range(len(kappas)):
	plt.plot(omega, D0[:] +epsilon*epsilon* D_eff[i, :])
plt.xscale("log")
plt.show()


