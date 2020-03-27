#from dolfin import *
import numpy as np
import dolfin as dlf

l = 10
Nx = 32
Ny = 32

p0 = dlf.Point(np.array([0.0, 0.0]))
p1 = dlf.Point(np.array([l, l]))
#creates points for the upper and lower corner 

mesh = dlf.RectangleMesh(p0, p1, Nx, Ny)
#creates a 2D mesh between the two corners

mark = {"generic": 0,
		"wall": 1,
		"lid": 2}

subdomains = dlf.MeshFunction("size_t", mesh, 1)
#lager en klasse av all punktene i mesh-et vårt

subdomains.set_all(mark["generic"])
#setter alle elementene i meshet til å være væskeelement = generic

class Wall(dlf.SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and (dlf.near(x[1], 0) or dlf.near(x[1], l) or dlf.near(x[0], 0))
		#if on_boundary is False then it will always return False
		#dlf.near checks if the x-position is on the boundary
		#both must be true for it to return true

class Lid(dlf.SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and dlf.near(x[1], l)
		#if on_boundary is False then it will always return False
		#dlf.near checks if the x-position is on the boundary
		#both must be true for it to return true

wall = Wall()
wall.mark(subdomains, mark["wall"])

lid = Lid()
lid.mark(subdomains, mark["lid"])

#does not work
#dlf.plot(subdomains, title="Subdomains")
#np.interactive()


#CG =  continuous-Galerkin
#We use dlf.triangle since our domain is discretized into subdomains, or elements, which are triangles 
#last argument is the order of the CG function
#we use second order on velocity and first order on pressure
#so velocity is more precise 
V = dlf.VectorElement("CG", dlf.triangle, 2)
P = dlf.FiniteElement("CG", dlf.triangle, 1)
W = dlf.FunctionSpace(mesh, V*P)

dx = dlf.Measure("dx", domain=mesh, subdomain_data=subdomains) # Volume integration
ds = dlf.Measure("ds", domain=mesh, subdomain_data=subdomains) # Surface integration

# Surface normal
n = dlf.FacetNormal(mesh)

#apply boundary conditions on velocity = V = W.sub(0)
#since the velocity is zero at the boundary noslip = (0,0)
noslip  = dlf.Constant((0.0, 0.0))
bc_wall = dlf.DirichletBC(W.sub(0), noslip, subdomains, mark["wall"])

#apply boundary conditions on velocity = V = W.sub(0) on lid
#since the velocity is constant on the lid
u0 = 1
slip  = dlf.Constant((u0, 0.0))
bc_lid = dlf.DirichletBC(W.sub(0), slip, subdomains, mark["lid"])

#apply boundary conditions on pressure = P = w.sub(1)
bc_p = dlf.DirichletBC(W.sub(1), 0, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")


#combine bounadry conditions
bcs = [bc_wall, bc_lid, bc_p]



#define som parameters
Reynolds = 1.0
nu = 1./Reynolds

# Solution vectors
w = dlf.Function(W)
v, q = dlf.TestFunctions(W)
u, p = (dlf.as_vector((w[0], w[1])), w[2])

F = dlf.inner(dlf.grad(u)*u, v)*dx \
		  + nu*dlf.inner(dlf.grad(u), dlf.grad(v))*dx \
		  - p*dlf.div(v)*dx \
		  - q*dlf.div(u)*dx

J = dlf.derivative(F, w) # Jacobian
dlf.solve(F == 0, w, bcs, J=J)


# Plot solution
#plot(u, title="Velocity")
#plot(p, title="Pressure")
#interactive()

#plot(dlf.sqrt(u**2), "Speed")
#interavtive()


