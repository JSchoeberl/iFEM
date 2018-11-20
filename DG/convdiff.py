import netgen.gui
from ngsolve import *
from netgen.geom2d import unit_square
mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

order = 3
fes1 = L2(mesh, order=order)
fes2 = FacetFESpace(mesh, order=order, dirichlet=".*", highest_order_dc=True)
fes = FESpace( [fes1,fes2])

b = CoefficientFunction( (y-0.5, 0.5-x) )

tau = 1e-3
eps = 1e-3
tend = 10

u,uhat = fes.TrialFunction()
v,vhat = fes.TestFunction()

h = specialcf.mesh_size
n = specialcf.normal(2)
alpha = 2

mstar = BilinearForm(fes)
mstar += SymbolicBFI(u*v)

adiff = BilinearForm(fes)
adiff += SymbolicBFI(eps*grad(u)*grad(v))
adiff += SymbolicBFI(eps*(-n*grad(u)*(v-vhat)-n*grad(v)*(u-uhat)),
                 element_boundary=True)
adiff += SymbolicBFI(eps*alpha*(order+1)**2/h*(u-uhat)*(v-vhat),element_boundary=True)


aconv = BilinearForm(fes, nonassemble=True)
aconv += SymbolicBFI(-b*u*grad(v))
uup = IfPos(b*n, u, u.Other(bnd=0))
aconv += SymbolicBFI(b*n*uup*v, element_boundary=True)


f = LinearForm(fes)

mstar.Assemble()
adiff.Assemble()
mstar.mat.AsVector().data += tau * adiff.mat.AsVector()
f.Assemble()


gfu = GridFunction(fes)
gfu.components[0].Set(exp(-10**2*((x-0.5)*(x-0.5) +(y-0.75)*(y-0.75))))
Draw(gfu.components[0], min=0, max=1, autoscale=False)

convu = gfu.vec.CreateVector()
w = gfu.vec.CreateVector()
r = gfu.vec.CreateVector()

inv = mstar.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
t = 0
with TaskManager():
  while t < tend:
    t += tau
    aconv.Apply(gfu.vec, convu)
    r.data = f.vec - convu - adiff.mat * gfu.vec
    w.data = inv * r
    gfu.vec.data += tau*w
    Redraw()
    
