import dolfin as df

#define periodic boundaries
class PeriodicBC(df.SubDomain):
    def __init__(self, Lx, Ly):
        self.Lx = Lx
        self.Ly = Ly
        df.SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        return bool(x[0] < df.DOLFIN_EPS_LARGE and
                    x[0] > -df.DOLFIN_EPS_LARGE and
                    on_boundary)

    def map(self, x, y):
        y[0] = x[0] - self.Lx
        y[1] = x[1]

#to check if on boundary node
class Wall(df.SubDomain):
    def __init__(self, Lx, Ly):
        self.Lx = Lx
        self.Ly = Ly
        df.SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        return on_boundary

#to check if one fluid node
class NotWall(df.SubDomain):
    def __init__(self, Lx, Ly):
        self.Lx = Lx
        self.Ly = Ly
        df.SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        return bool(x[0] < df.DOLFIN_EPS_LARGE or
                    x[0] > self.Lx - df.DOLFIN_EPS_LARGE) and on_boundary
