#from dolfin import *
import numpy as np
import dolfin as dlf

delta = 0.5
l = 10
D = 4

Nx = int(l/delta)
Ny = int(D/delta)

p0 = dlf.Point(np.array([0.0, 0.0]))
p1 = dlf.Point(np.array([l, D]))
#creates points for the upper and lower corner 

mesh = dlf.RectangleMesh(p0, p1, Nx, Ny)
#creates a 2D mesh between the two corners

mark = {"generic": 0,
		"wall": 1,
		"left": 2,
		"right": 3 }

subdomains = dlf.MeshFunction("size_t", mesh, 1)
#lager en klasse av all punktene i mesh-et vårt

subdomains.set_all(mark["generic"])
#setter alle elementene i meshet til å være væskeelement = generic

#lager klasse for venstre vegg
class Left(dlf.SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and dlf.near(x[0], 0)
		#if on_boundary is False then it will always return False
		#dlf.near checks if the x-position is on the boundary
		#both must be true for it to return true

class Right(dlf.SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and dlf.near(x[0], l)
		#if on_boundary is False then it will always return False
		#dlf.near checks if the x-position is on the boundary
		#both must be true for it to return true

class Wall(dlf.SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and dlf.near(x[1], 0) or dlf.near(x[1], D)
		#if on_boundary is False then it will always return False
		#dlf.near checks if the x-position is on the boundary
		#both must be true for it to return true

left = Left()
left.mark(subdomains, mark["left"])

right = Right()
right.mark(subdomains, mark["right"])

wall = Wall()
wall.mark(subdomains, mark["wall"])

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

# Define variational problem
(u, p) = dlf.TrialFunctions(W)
(v, q) = dlf.TestFunctions(W)

dx = dlf.Measure("dx", domain=mesh, subdomain_data=subdomains) # Volume integration
ds = dlf.Measure("ds", domain=mesh, subdomain_data=subdomains) # Surface integration

# Surface normal
n = dlf.FacetNormal(mesh)

# Pressures. First define the numbers (for later use):
P_left = 1.0
P_right = 0.0

# ...and then the DOLFIN constants:
pressure_left  = dlf.Constant(P_left)
pressure_right = dlf.Constant(P_right)


# Body force:
force = dlf.Constant((0.0, 0.0))

#give dlf. the variational forms of our equations and boundary conditions
#dx gives the volume integration, and ds gives the surface integration

a = dlf.inner(dlf.grad(u), dlf.grad(v))*dx - p*dlf.div(v)*dx + q*dlf.div(u)*dx
L = dlf.inner(force, v) * dx \
		- pressure_left  * dlf.inner(n, v) * ds(mark["left"]) \
		- pressure_right * dlf.inner(n, v) * ds(mark["right"])

#apply boundary conditions on velocity = V = W.sub(0)
#since the velocity is zero at the boundary noslip = (0,0)
noslip  = dlf.Constant((0.0, 0.0))
bc_wall = dlf.DirichletBC(W.sub(0), noslip, subdomains, mark["wall"])

#apply boundary conditions on pressure = P = w.sub(1)
bc_right = dlf.DirichletBC(W.sub(1), P_right, subdomains, mark["right"])
bc_left  = dlf.DirichletBC(W.sub(1), P_left, subdomains,  mark["left"])

#combine bounadry conditions
bcs = [bc_wall, bc_left, bc_right]

# Compute solution
w = dlf.Function(W)
dlf.solve(a == L, w, bcs)
(u, p) = w.split(True)

# Plot solution
#plot(u, title="Velocity")
#plot(p, title="Pressure")
#interactive()

#plot(dlf.sqrt(u**2), "Speed")
#interavtive()


