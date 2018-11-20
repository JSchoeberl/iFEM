from ngsolve import *
from netgen.geom2d import unit_square
mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))


b = CoefficientFunction( (y-0.5, 0.5-x) )
Draw (b, mesh, "wind")

fes = L2(mesh, order=3)

u = fes.TrialFunction()
v = fes.TestFunction()
a = BilinearForm(fes, nonassemble=True)

a += SymbolicBFI(-b*u*grad(v))

# the upwind-term:
n = specialcf.normal(2)
uup = IfPos(b*n, u, u.Other(bnd=0))
a += SymbolicBFI(b*n*uup*v, element_boundary=True)

f = LinearForm(fes)
f.Assemble()

gfu = GridFunction(fes)
gfu.Set(exp(-10**2*((x-0.5)*(x-0.5) +(y-0.75)*(y-0.75))))

Draw(gfu, min=0, max=1, autoscale=False)


tau = 0.001
tend = 10
t = 0

# we need a help vector
w = gfu.vec.CreateVector()

with TaskManager():
    while t < tend:
        # apply the transport operator 
        a.Apply (gfu.vec, w)
        
        # use an efficient, matrix-free technique 
        # to solve with the mass matrix
        fes.SolveM (rho=CoefficientFunction(1), vec=w)
        gfu.vec.data -= tau * w
        t += tau
        Redraw()
