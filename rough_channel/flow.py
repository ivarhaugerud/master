import dolfin as df
import argparse
import matplotlib.pyplot as plt
import numpy as np
from mesh_box import square_rough_mesh
import os
from bcs import PeriodicBC, Wall, NotWall
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
folder = "data_square"


# Form compiler options
df.parameters["form_compiler"]["optimize"] = True
df.parameters["form_compiler"]["cpp_optimize"] = True

# Command line arguments
parser = argparse.ArgumentParser(description="Taylor dispersion")
parser.add_argument("-res", type=int, default=32, help="Resolution")
parser.add_argument("-R", type=float, default=0.0,
                    help="Imposed Reynolds number")
parser.add_argument("-logPe_min", type=float, default=0.0, help="Min Pe")
parser.add_argument("-logPe_max", type=float, default=5.0, help="Max Pe")
parser.add_argument("-logPe_N", type=int, default=10, help="Num Pe")
parser.add_argument("-b", type=float, default=0.5, help="Roughness b")
parser.add_argument("--onlyflow", action="store_true", help="only flow")
args = parser.parse_args()

# Calculate resolution
N = args.res
R = df.Constant(args.R)
Ly = 2.0
dx = Ly/N

#define mesh and coordinates
mesh = square_rough_mesh(args.b, dx)
coords = mesh.coordinates()[:]
Lx = df.MPI.max(df.MPI.comm_world, coords[:, 0].max())

# define lagrange elements for velocity (Eu), pressure (Ep) and Brenner field (Echi)
Eu = df.VectorElement("Lagrange", mesh.ufl_cell(), 3)
Ep = df.FiniteElement("Lagrange", mesh.ufl_cell(), 2)
Echi = df.FiniteElement("Lagrange", mesh.ufl_cell(), 2)

#create periodic boundary conditions and define a boundary-function
pbc = PeriodicBC(Lx, Ly)
wall = Wall(Lx, Ly)
notwall = NotWall(Lx, Ly)


#set wall points to 1, and fluid points to 0
subd = df.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
subd.set_all(0)
wall.mark(subd, 1)
notwall.mark(subd, 0)

#create folder for saving
if rank == 0 and not os.path.exists(folder):
    os.makedirs(folder)

with df.XDMFFile(mesh.mpi_comm(), "{}/subd_b{}.xdmf".format(
        folder, args.b)) as xdmff:
    xdmff.write(subd)

#define problem for u and P
W = df.FunctionSpace(mesh, df.MixedElement([Eu, Ep]),
                     constrained_domain=pbc)
#define problem for Brenner
S = df.FunctionSpace(mesh, Echi, constrained_domain=pbc)

#define trial and test functions
w_ = df.Function(W)
u_, p_ = df.split(w_)
w = df.TrialFunction(W)
u, p = df.split(w)
v, q = df.TestFunctions(W)

#define external force
f = df.Constant((1., 0.))

#NS on variational form, including continuity equation
F = (
    R*df.inner(df.grad(u_)*u_, v)*df.dx
    + df.inner(df.grad(u_), df.grad(v))*df.dx
    - df.div(v)*p_*df.dx - df.div(u_)*q*df.dx
    - df.dot(f, v)*df.dx
)

#for parallelization
x0, y0 = coords[0, 0], coords[0, 1]
x0 = comm.bcast(x0, root=0)
y0 = comm.bcast(y0, root=0)

#calculate jacobian
J = df.derivative(F, w_, du=w)

#define no-slip boundary conditions for fluid
bcu = df.DirichletBC(W.sub(0), df.Constant((0., 0.)),
                     subd, 1)

#boundary condition for pressure
bcp = df.DirichletBC(W.sub(1), df.Constant(0.),
                     ("abs(x[0]-({x0})) < DOLFIN_EPS && "
                      "abs(x[1]-({y0})) < DOLFIN_EPS").format(x0=x0, y0=y0),
                     "pointwise")

#total BC is BC for velocity field and pressure
bcs = [bcu, bcp]

#define problem for velocity field
problem = df.NonlinearVariationalProblem(F, w_, bcs=bcs, J=J)
solver = df.NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-12
solver.parameters["newton_solver"]["linear_solver"] = "mumps"
solver.parameters["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1e-14
solver.parameters["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True

#solve problem
solver.solve()

#define a function with a value of 1 over the fluid domain, and calculate the domain volume from integrating this function
one = df.interpolate(df.Constant(1.), S)
V_Omega = df.assemble(one*df.dx)

#calculate mean velocity value, y-comp will not contribute
ux_mean = df.assemble(u_[0]*df.dx)/V_Omega

#print mean velocity
if rank == 0:
    print("ux_mean = {}".format(ux_mean))

#calculate the aposteriori Reynolds number
Re = args.R*ux_mean
if rank == 0:
    print("Re = {}".format(Re))


###
# NOW THE VELOCITY FIELD IS KNOWN, AND THE BRENNER FIELD CAN BE CALCULATED
###
# Part 2: chi

#rename velocity field and pressure, then normalize velocity field to have an average value of 1
U_, P_ = w_.split(deepcopy=True)
U_.rename("u", "tmp")
U_.vector()[:] /= ux_mean

#set up for writing of data
with df.XDMFFile(mesh.mpi_comm(), "{}/U_Re{}_b{}.xdmf".format(
        folder, Re, args.b)) as xdmff:
    xdmff.write(U_)

with df.HDF5File(mesh.mpi_comm(), "{}/flow_Re{}_b{}.h5".format(
        folder, Re, args.b), "w") as h5f:
    h5f.write(mesh, "mesh")
    h5f.write(U_, "U")
    h5f.write(subd, "subd")

#if the Brenner field is of no importance, stop here
if args.onlyflow:
    exit("Only flow!")

#calculate normal vector of boundary
n = df.FacetNormal(mesh)

#calculate Peclet number
Pe = df.Constant(1.0)

#define trial and test functions using mesh for B-field (S)
chi  = df.TrialFunction(S)
chi_ = df.Function(S)
psi  = df.TestFunction(S)

#define surface integral
ds = df.Measure("ds", domain=mesh, subdomain_data=subd)

#variational problem for Brenner field, where the Neumann BCs are included
F_chi = (n[0]*psi*ds(1)
         + df.inner(df.grad(chi), df.grad(psi))*df.dx
         + Pe*psi*df.dot(U_, df.grad(chi))*df.dx
         + Pe*(U_[0] - df.Constant(1.))*psi*df.dx)

#define left and right hand side
a_chi, L_chi = df.lhs(F_chi), df.rhs(F_chi)

#define problem and solver
problem_chi2 = df.LinearVariationalProblem(a_chi, L_chi, chi_, bcs=[])
solver_chi2  = df.LinearVariationalSolver(problem_chi2)
solver_chi2.parameters["krylov_solver"]["absolute_tolerance"] = 1e-15

#array of Peclet numbers to investigate
logPe_arr = np.linspace(args.logPe_min, args.logPe_max, args.logPe_N)

if rank == 0:
    data = np.zeros((len(logPe_arr), 3))

#solve the Problem for each Peclet number and print
for iPe, logPe in enumerate(logPe_arr):
    Pe_loc = 10**logPe
    Pe.assign(Pe_loc)

    solver_chi2.solve()

    integral = (2*df.assemble(chi_.dx(0)*df.dx)/V_Omega
                + df.assemble(df.inner(df.grad(chi_),
                                       df.grad(chi_))*df.dx)/V_Omega)

    if rank == 0:
        print("Pe = {}, D_eff/D = {}".format(Pe_loc, 1+integral))

        data[iPe, 0] = Pe_loc
        data[iPe, 1] = 1+integral
        data[iPe, 2] = integral/Pe_loc**2

#save data
if rank == 0:
    np.savetxt("{}/Re{}_b{}_res{}.dat".format(folder, str(Re)[:5], args.b, args.res), data)
