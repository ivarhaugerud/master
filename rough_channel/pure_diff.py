import dolfin as df
from mesh_box import square_rough_mesh
from bcs import PeriodicBC, Wall, NotWall
import numpy as np
import argparse
import os
from mpi4py import MPI

### CODE IS SIMILAR TO flow.py, but only calculates for Pe=0, so that the velocity field is irrelevant
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
folder = "data_square"


# Form compiler options
df.parameters["form_compiler"]["optimize"] = True
df.parameters["form_compiler"]["cpp_optimize"] = True

#command line arguments
parser = argparse.ArgumentParser(description="Pure diffusion")
parser.add_argument("-res", type=int, default=32, help="Resolution")
parser.add_argument("-b_min", type=float, default=0.0,
                    help="Roughness b min")
parser.add_argument("-b_max", type=float, default=1.5,
                    help="Roughness b max")
parser.add_argument("-b_N", type=int, default=15,
                    help="Roughness b N")
args = parser.parse_args()

#implicitly defined resolution of mesh
N = args.res
Ly = 2.0
dx = Ly/N

#data depending on number of roughness values investigated
data = np.zeros((args.b_N, 2))

#loop over b-values
for ib, b in enumerate(np.linspace(args.b_min, args.b_max, args.b_N)):
    #get mesh and coordinates for for b
    mesh = square_rough_mesh(b, dx)
    coords = mesh.coordinates()[:]

    #get system length
    Lx = df.MPI.max(df.MPI.comm_world, coords[:, 0].max())

    #define lagrange elements for Brenner field
    Echi = df.FiniteElement("Lagrange", mesh.ufl_cell(), 2)

    #define periodic BCs, and if a point is a wall point or not
    pbc = PeriodicBC(Lx, Ly)
    wall = Wall(Lx, Ly)
    notwall = NotWall(Lx, Ly)

    #define a function on the discrete mesh points with a value of 1 at wall, and a value of 0 at the fluid nodes
    subd = df.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    subd.set_all(0)
    wall.mark(subd, 1)
    notwall.mark(subd, 0)

    #define periodic function of the Brenner field
    S = df.FunctionSpace(mesh, Echi, constrained_domain=pbc)

    #calculate fluid volume
    one = df.interpolate(df.Constant(1.), S)
    V_Omega = df.assemble(one*df.dx)

    #calculate normal vector to the edges of the mesh
    n = df.FacetNormal(mesh)

    #define trial and test function
    chi  = df.TrialFunction(S)
    chi_ = df.Function(S, name="chi")
    psi  = df.TestFunction(S)

    #define surface integral
    ds = df.Measure("ds", domain=mesh, subdomain_data=subd)

    #Brenner equation with Pe=0 on variational form, Neumann boundary conditions are included
    F_chi = (n[0]*psi*ds(1)
             + df.inner(df.grad(chi), df.grad(psi))*df.dx)

    #LHS and RHS of above expression
    a_chi, L_chi = df.lhs(F_chi), df.rhs(F_chi)

    #define problem
    problem_chi2 = df.LinearVariationalProblem(a_chi, L_chi, chi_, bcs=[])
    solver_chi2 = df.LinearVariationalSolver(problem_chi2)
    solver_chi2.parameters["krylov_solver"]["absolute_tolerance"] = 1e-15

    #solve problem
    solver_chi2.solve()

    #save data
    with df.XDMFFile(mesh.mpi_comm(),
                     "chi_Pe0_b{}.xdmf".format(b)) as xdmff:
        xdmff.write(chi_)
    
    #calculate integral of solution, which gives dispersion tensor
    integral = (2*df.assemble(chi_.dx(0)*df.dx)/V_Omega
                + df.assemble(df.inner(df.grad(chi_),
                                       df.grad(chi_))*df.dx)/V_Omega)
    #print effective diffusion
    if rank == 0:
        print("b = {}, D_eff/D = {}".format(b, 1+integral))

    #format data file
    data[ib, 0] = b
    data[ib, 1] = 1+integral

#save data to file
if rank == 0:
    np.savetxt("{}/Pe0_square_rough_mesh_res{}_bmax{}.dat".format(folder, N, args.b_max), data)




