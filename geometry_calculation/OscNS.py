import argparse
import os

import dolfin as df
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
basefolder = "results_oscwavychannel"

# Form compiler options
df.parameters["form_compiler"]["optimize"] = True
df.parameters["form_compiler"]["cpp_optimize"] = True


class PBC(df.SubDomain):
    def __init__(self, Lx):
        self.Lx = Lx
        df.SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two slave edges
        return bool(df.near(x[0], 0) and on_boundary)

    def map(self, x, y):
        if df.near(x[0], self.Lx):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
        else:  # near(x[2], Lz/2.):
            y[0] = x[0]
            y[1] = x[1]


class Walls(df.SubDomain):
    def __init__(self, Lx, Ly, epsilon):
        self.Lx = Lx
        self.Ly = Ly
        self.eps = epsilon
        self.tol = 1e-3
        df.SubDomain.__init__(self)

    def inside(self, x, on_bnd):
        return (
            x[1] < -self.Ly / 2. *
            (1.0 + self.eps * np.cos(2 * np.pi * x[0] / self.Lx)) + self.tol
            or x[1] > self.Ly / 2. *
            (1.0 + self.eps * np.cos(2 * np.pi * x[0] / self.Lx)) - self.tol)


def inlet(x, on_bnd):
    return on_bnd and df.near(x[0], 0)


parser = argparse.ArgumentParser(
    description="Time-dependent Taylor dispersion")
parser.add_argument("-Lx",   type=float, default=2.0,   help="Length of cell")
parser.add_argument("-res",  type=int,   default=32,    help="Resolution")
parser.add_argument("-dt",   type=float, default=0.02,  help="Time step")
parser.add_argument("-T",    type=float, default=100.0, help="Total time")
parser.add_argument("-tau",  type=float, default=5.0,   help="Period")
parser.add_argument("-nu",   type=float, default=1.0,   help="Kinematic viscosity")
parser.add_argument("-D",    type=float, default=0.1,   help="Molecular diffusivity")
parser.add_argument("-f0",   type=float, default=0.0,   help="Base force")
parser.add_argument("-f1",   type=float, default=1.0,   help="Oscillatory force")

parser.add_argument("-tol", "--tolerance", type=float, default=1e-5, help="Convergence tolerance")
parser.add_argument("-eps", "--epsilon",   type=float, default=0.2,  help="Roughness amplitude")

parser.add_argument("--disable_inertia", action="store_true", help="Disable inertia")
parser.add_argument("--onlyflow",        action="store_true", help="only flow")
args = parser.parse_args()

Ny = args.res
D = df.Constant(args.D)
nu = df.Constant(args.nu)
enable_inertia = not args.disable_inertia
tau = args.tau
f0 = args.f0
f1 = args.f1
T = args.T

Lx = args.Lx
Ly = 2.0
eps = args.epsilon

dt = args.dt
dx = Ly / Ny
Nx = 4*args.res#int(Lx / dx)
tol = args.tolerance

folder = os.path.join(
    basefolder, "Lx{}_tau{}_eps{}_nu{}_D{}_fzero{}_fone{}_res{}_dt{}".format(
        Lx, tau, eps, args.nu, args.D, f0, f1, Ny, dt))
timestampsfolder = os.path.join(folder, "timestamps")

if not os.path.exists(folder):
    os.makedirs(folder)
if not os.path.exists(timestampsfolder):
    os.makedirs(timestampsfolder)

mesh = df.RectangleMesh(df.Point(0., -Ly / 2), df.Point(Lx, Ly / 2), Nx, Ny)
x = mesh.coordinates()[:]
x[:, 1]  = (np.arctan(2. * np.pi * x[:, 1] / Ly) / np.arctan(1. * np.pi))
x[:, 1] *= (1 + eps * np.cos(2 * np.pi * x[:, 0] / Lx))

Eu   = df.VectorElement("Lagrange", mesh.ufl_cell(), 2)
Ep   = df.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
Echi = df.FiniteElement("Lagrange", mesh.ufl_cell(), 1)

pbc   = PBC(Lx)
walls = Walls(Lx, Ly, eps)
inlet = df.AutoSubDomain(inlet)

subd = df.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subd.set_all(0)
walls.mark(subd, 1)
inlet.mark(subd, 2)

if rank == 0 and not os.path.exists(folder):
    os.makedirs(folder)

with df.XDMFFile(mesh.mpi_comm(), "{}/subd.xdmf".format(folder)) as xdmff:
    xdmff.write(subd)

W = df.FunctionSpace(mesh, df.MixedElement([Eu, Ep]), constrained_domain=pbc)
S = df.FunctionSpace(mesh, Echi, constrained_domain=pbc)

ds = df.Measure("ds", domain=mesh, subdomain_data=subd)

w_ = df.Function(W)
w_1 = df.Function(W)
w_2 = df.Function(W)
u_, p_ = df.split(w_)
u_1, p_1 = df.split(w_1)
u_2, p_2 = df.split(w_2)
w = df.TrialFunction(W)
u, p = df.split(w)
v, q = df.TestFunctions(W)

f = df.Expression(("f_0 + f_1*cos(2*M_PI*t/tau)", "0."),
                  f_0=f0,
                  f_1=f1,
                  t=0.,
                  tau=tau,
                  degree=2)

u_CN = 0.5 * (u + u_1)
u_CN_ = 0.5 * (u_ + u_1)
u_AB_ = 1.5 * u_1 + 0.5 * u_2

F = (df.dot(u - u_1, v) / dt * df.dx +
     nu * df.inner(df.grad(u_CN), df.grad(v)) * df.dx - df.div(v) * p * df.dx -
     df.div(u) * q * df.dx - df.dot(f, v) * df.dx)
if enable_inertia:
    F += df.inner(df.grad(u_CN) * u_AB_, v) * df.dx

bcu = df.DirichletBC(W.sub(0), df.Constant((0., 0.)), subd, 1)

x0, y0 = x[0, 0], x[0, 1]
x0 = comm.bcast(x0, root=0)
y0 = comm.bcast(y0, root=0)
# distribute!

bcp = df.DirichletBC(W.sub(1), df.Constant(0.),
                     ("abs(x[0]-({x0})) < DOLFIN_EPS && "
                      "abs(x[1]-({y0})) < DOLFIN_EPS").format(x0=x0, y0=y0),
                     "pointwise")

bcs = [bcu, bcp]

a = df.lhs(F)
L = df.rhs(F)

n = df.FacetNormal(mesh)
B = df.TrialFunction(S)
B_ = df.Function(S, name="B")
B_1 = df.Function(S)
psi = df.TestFunction(S)

B_CN = 0.5 * (B + B_1)

ux_CN_mean = df.Expression("ux", degree=1, ux=0.)

F_chi = ((B - B_1) * psi / dt * df.dx +
         D * df.inner(df.grad(B_CN), df.grad(psi)) * df.dx +
         psi * df.dot(u_CN_, df.grad(B_CN)) * df.dx -
         (u_CN_[0] - ux_CN_mean) * psi * df.dx - n[0] * psi * ds(1))

a_chi, L_chi = df.lhs(F_chi), df.rhs(F_chi)

problem = df.LinearVariationalProblem(a, L, w_, bcs=bcs)
solver = df.LinearVariationalSolver(problem)
#solver.parameters["absolute_tolerance"] = 1e-14
#solver.parameters["krylov_solver"]["absolute_tolerance"] = 1e-14
#solver.parameters["krylov_solver"]["nonzero_initial_guess"] = True

problem_B = df.LinearVariationalProblem(a_chi, L_chi, B_, bcs=[])
solver_B = df.LinearVariationalSolver(problem_B)
solver_B.parameters["krylov_solver"]["absolute_tolerance"] = 1e-14

t = 0

one = df.interpolate(df.Constant(1.), S)
V_Omega = df.assemble(one * df.dx)
r = df.SpatialCoordinate(mesh)

xdmf_u = df.XDMFFile(mesh.mpi_comm(), "{}/u.xdmf".format(folder))
xdmf_B = df.XDMFFile(mesh.mpi_comm(), "{}/B.xdmf".format(folder))

ux_mean = 0.

Ntau = int(tau / dt)
normu_arr = np.zeros(Ntau)
it = 0

if rank == 0:
    ofile = open("{}/tdata.dat".format(folder), "w")

normu_max = 0.
nt_conv = 0
it_conv = -1

eq_data = np.zeros((Ntau, 8))
if rank == 0:
    timestampsfile = open(os.path.join(timestampsfolder, "timestamps.dat"),
                          "w")

while t < T:
    t += dt
    it += 1
    # Assign
    w_2.assign(w_1)
    w_1.assign(w_)
    B_1.assign(B_)
    f.t = t - 0.5 * dt

    solver.solve()
    ux_mean_prev = ux_mean
    ux_mean = df.assemble(u_[0] * df.dx) / V_Omega
    uy_mean = df.assemble(u_[1] * df.dx) / V_Omega
    u2_mean = df.assemble(df.dot(u_, u_) * df.dx) / V_Omega

    ux_CN_mean.ux = 0.5 * (ux_mean + ux_mean_prev)
    solver_B.solve()

    C2 = df.assemble(((B_ - r[0])**2 - (B_1 - r[0])**2) * df.dx) / V_Omega
    B2 = df.assemble((B_**2 - B_1**2) * df.dx) / V_Omega
    Bmean = df.assemble(B_ * df.dx) / V_Omega

    integral_2 = (
        1 - 2 * df.assemble(B_.dx(0) * df.dx) / V_Omega +
        df.assemble(df.inner(df.grad(B_), df.grad(B_)) * df.dx) / V_Omega)

    ires = it % Ntau
    normu = df.norm(w_.vector())
    normu_max = max(normu_max, abs(normu))
    res = abs(normu - normu_arr[ires]) / normu_max
    normu_arr[ires] = normu
    if res < tol:
        nt_conv += 1
    else:
        nt_conv = 0

    if nt_conv == Ntau:
        it_conv = it

    f.t = t
    fxt = df.assemble(f[0] * df.dx(domain=mesh)) / V_Omega
    if rank == 0:
        items = (t, fxt, ux_mean, uy_mean, u2_mean, C2, B2, Bmean, integral_2,
                 normu, res)
        print("time =", t, "\t ux_mean =", ux_mean, "\t res =", res, "vs.",
              tol)
        ofile.write(("\t".join(["{}" for _ in items]) + "\n").format(*items))

    U_, P_ = w_.split(deepcopy=True)
    U_.rename("u", "tmp")
    P_.rename("p", "tmp")
    xdmf_u.write(U_, t)
    xdmf_B.write(B_, t)
    if it_conv > 0:
        if it > it_conv + Ntau:
            break
        with df.HDF5File(mesh.mpi_comm(),
                         "{}/up_{}.h5".format(timestampsfolder,
                                              ires), "w") as h5f:
            h5f.write(U_, "u")
            h5f.write(P_, "p")
        eq_data[ires, 0] = ires * dt
        eq_data[ires, 1] = ux_mean
        eq_data[ires, 2] = uy_mean
        eq_data[ires, 3] = u2_mean
        eq_data[ires, 4] = C2
        eq_data[ires, 5] = integral_2
        if rank == 0:
            timestampsfile.write("{} up_{}.h5\n".format(ires * dt, ires))

if rank == 0:
    print("Time-and-space averaged:")
    print("  u_mean_avg =", np.mean(eq_data[:, 1]))
    u_rms = np.sqrt(np.mean(eq_data[:, 3]))
    print("  u_rms =", u_rms)
    print("A posteriori estimates:")
    print(
        "  Re =", u_rms / args.nu, "(inertia was enabled)"
        if enable_inertia else "(but inertia was disabled)")
    print("  Pe =", u_rms / args.D)
    print("  Sc =", args.nu / args.D)
    print("Finally:")
    print("  D_eff/D =", np.mean(eq_data[:, 5]))

    dC2dt = (np.roll(eq_data[:, 4], -1) - np.roll(eq_data[:, 4], 1)) / (2 * dt)
    eq_data[:, 6] = dC2dt
    eq_data[:, 7] = 0.5 * dC2dt + args.D * eq_data[:, 5]
    np.savetxt(os.path.join(folder, "eqdata.dat"), eq_data)
    timestampsfile.close()

if rank == 0:
    ofile.close()
