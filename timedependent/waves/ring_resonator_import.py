import math
from ngsolve import *
import netgen.geom2d as gm

SetHeapSize(50*1000*1000)
# geometrical and material properties
geo = gm.SplineGeometry()

xneg  =-0.43
xpos  = 0.43
yneg  =-0.48
ypos  = 0.48
wslab = 0.04
cringx = 0.0
cringy = 0.0
rring = 0.4
gap   = 0.005

pntx = [xneg,xpos]
pnty = [yneg,-rring-gap-wslab,-rring-gap,rring+gap,rring+gap+wslab,ypos]


pts = []
for yi in pnty:
    for xi in pntx:
        pts.append (geo.AddPoint(xi,yi))

#### parameters for source position


### geometry #######################################################
#inner rects
geo.Append (["line", pts[0], pts[1]], leftdomain=1, rightdomain=0)
geo.Append (["line", pts[1], pts[3]], leftdomain=1, rightdomain=0)
geo.Append (["line", pts[3], pts[2]], leftdomain=1, rightdomain=2)
geo.Append (["line", pts[2], pts[0]], leftdomain=1, rightdomain=0)

geo.Append (["line", pts[3], pts[5]], leftdomain=2, rightdomain=0,bc="normal_wg_rightbottom")
geo.Append (["line", pts[5], pts[4]], leftdomain=2, rightdomain=3)
geo.Append (["line", pts[4], pts[2]], leftdomain=2, rightdomain=0,bc="normal_wg_leftbottom")

geo.Append (["line", pts[5], pts[7]], leftdomain=3, rightdomain=0)
geo.Append (["line", pts[7], pts[6]], leftdomain=3, rightdomain=4)
geo.Append (["line", pts[6], pts[4]], leftdomain=3, rightdomain=0)

geo.Append (["line", pts[7], pts[9]], leftdomain=4, rightdomain=0,bc="normal_wg_righttop")
geo.Append (["line", pts[9], pts[8]], leftdomain=4, rightdomain=5)
geo.Append (["line", pts[8], pts[6]], leftdomain=4, rightdomain=0,bc="normal_wg_lefttop")

geo.Append (["line", pts[9], pts[11]], leftdomain=5, rightdomain=0)
geo.Append (["line", pts[11], pts[10]], leftdomain=5, rightdomain=0)
geo.Append (["line", pts[10], pts[8]], leftdomain=5, rightdomain=0)

geo.AddCircle(c=(cringx,cringy), r=rring, leftdomain=6, rightdomain=3)
geo.AddCircle(c=(cringx,cringy), r=rring-wslab, leftdomain=7, rightdomain=6)

for i in (1,3,5,7):
    geo.SetMaterial(i, "air")
for i in (2,4,6):
    geo.SetMaterial(i, "eps_nine")
data = geo.CreatePML(0.05)
normals = data["normals"]


mesh = Mesh(geo.GenerateMesh(maxh=0.05))
mesh.Curve(3)
Draw(mesh)

eps_r = {"air" : 1,
         "eps_nine" : 3**3}

order = 3

for mat in mesh.GetMaterials():
    if mat.startswith("pml_normal_wg"):
        eps_r[mat] = eps_r["eps_nine"]

for mat in mesh.GetMaterials():
    if mat not in eps_r:
        eps_r[mat] = eps_r["air"]

# print ('materials:',[mat for mat in mesh.GetMaterials()])

def gaussp(pt):
    return CoefficientFunction(sin( (math.pi/wslab)*(y-pt[1])))
    return CoefficientFunction(exp ( -3200 * ( (y-pt[1])**2 )) )

# ~ center_source = (pntx[0],0.5*pnty[3]+0.5*pnty[4])
center_source = (pntx[0],pnty[3])


### Parameters for Source field ##########################################################################
wavelength = 1.542
fcen       = 5/wavelength
df         = 0.1
tpeak      = 1

fes_facet = FacetFESpace(mesh, order=order+1)
gfsource = GridFunction(fes_facet)
source_cf = gaussp(center_source)
gfsource.Set(source_cf,definedon=mesh.Boundaries("normal_wg_lefttop"))

fes_test = H1(mesh)
gf_test = GridFunction(fes_test)
gf_test.Set(source_cf)
Draw(gf_test)

##########################################################################################################


fes_u = VectorL2(mesh, order=order, piola=True, order_policy=ORDER_POLICY.CONSTANT)
fes_p = L2(mesh, order=order+1, all_dofs_together=True, order_policy=ORDER_POLICY.CONSTANT)
fes_tr = FacetFESpace(mesh, order=order+1)
fes_hdiv = HDiv(mesh, order=order+1, orderinner=1)
fes = FESpace( [fes_p,fes_p,fes_u, fes_u] )

p,phat, u,uhat = fes.TrialFunction()
q,qhat, v,vhat = fes.TestFunction()

n = specialcf.normal(2) 
h = specialcf.mesh_size

Bel = BilinearForm(trialspace=fes_p, testspace=fes_u, geom_free=True)
Bel += SymbolicBFI ( grad(fes_p.TrialFunction())*fes_u.TestFunction())
Bel += SymbolicBFI ( -fes_p.TrialFunction()*(fes_u.TestFunction()*n), element_boundary=True)
Bel.Assemble()

Btr = BilinearForm(trialspace=fes_tr, testspace=fes_u, geom_free=True)
Btr += SymbolicBFI (0.5 * fes_tr.TrialFunction() * (n*fes_u.TestFunction()), element_boundary=True)
Btr.Assemble()

Bstab = BilinearForm(trialspace=fes_p, testspace=fes_hdiv, geom_free=True)
Bstab += SymbolicBFI(fes_p.TrialFunction()*(fes_hdiv.TestFunction()*n), element_boundary=True)
Bstab.Assemble()

Mstab = BilinearForm(fes_hdiv)
Mstab += SymbolicBFI(fes_hdiv.TrialFunction() * fes_hdiv.TestFunction())
Mstab.Assemble()
# Mstabinv = Mstab.mat.Inverse(inverse="sparsecholesky")
Mstabinv = Mstab.mat.CreateSmoother()


### LinearForm for the source ##############################################
Lsrc  = LinearForm(fes)
# ~ Lsrc  += SymbolicLFI( gfsource*(n*v),element_boundary=True)
# ~ Lsrc  += SymbolicLFI( eps*SrcField*q)
Lsrc  += SymbolicLFI( gfsource*q,element_boundary=True)
Lsrc.Assemble()
############################################################################


nvec = { mat : ((normals[mat][0], normals[mat][1]) if mat in normals else (0,0)) for mat in mesh.GetMaterials() }

cfn = CoefficientFunction( [CoefficientFunction(nvec[mat]) for mat in mesh.GetMaterials()])
cft = CoefficientFunction( ( cfn[1], -cfn[0] ) )

pml1d = mesh.Materials("pml_default.*|pml_normal.*")

print(pml1d.Mask())

eps = CoefficientFunction([eps_r[mat] for mat in mesh.GetMaterials()])

Draw(eps, mesh, "eps")


# damping matrices
sigma = 10   # pml damping parameter
dampp1 = fes_p.Mass (eps, definedon=pml1d)
dampp2 = fes_p.Mass (eps, definedon=mesh.Materials("pml_corner"))
dampu1 = fes_u.Mass (OuterProduct(cfn,cfn), definedon=pml1d)
dampu2 = fes_u.Mass (OuterProduct(cft,cft), definedon=pml1d)

gfu = GridFunction(fes)
gfu.vec[:] = 0
# gfu.components[0].Set (gaussp(prism1_center) + gaussp(mirrory(prism1_center)))

scalval = 1e-1
Draw(gfu.components[0], mesh, "p")

# For the time stepping
tau = 2e-4
tend = 100
t = 0

pdofs = fes.Range(0)
phatdofs = fes.Range(1)
udofs = fes.Range(2)
uhatdofs = fes.Range(3)

w = gfu.vec.CreateVector()
hv = gfu.vec.CreateVector()
hvtrace = BaseVector(fes_tr.ndof)
hvu = BaseVector(fes_u.ndof)
hvu2 = BaseVector(fes_u.ndof)

emb_u = Embedding(fes.ndof, udofs)
emb_uhat = Embedding(fes.ndof, uhatdofs)
emb_p = Embedding(fes.ndof, pdofs)
emb_phat = Embedding(fes.ndof, phatdofs)

traceop = fes_p.TraceOperator(fes_tr, False)
fullB = emb_u @ (Bel.mat + Btr.mat @ traceop) @ emb_p.T

dampingu = emb_u @ dampu1 @ emb_u.T + (-emb_u + emb_uhat) @ dampu2 @ (emb_u.T + emb_uhat.T)
dampingp = emb_p @ dampp1 @ emb_p.T + emb_p @ dampp2 @ (2*emb_p.T-emb_phat.T) + emb_phat @ dampp2 @ emb_p.T

invmassp = fes_p.Mass(eps).Inverse()
invmassu = fes_u.Mass(Id(mesh.dim)).Inverse()
invp = emb_p @ invmassp @ emb_p.T + emb_phat @ invmassp @ emb_phat.T
invu = emb_u @ invmassu @ emb_u.T + emb_uhat @ invmassu @ emb_uhat.T


gfstab = GridFunction(fes_hdiv)
hvstab = gfstab.vec.CreateVector()

from ngsolve.internal import visoptions, VideoStart, VideoAddFrame, VideoFinalize, SnapShot
visoptions.mmaxval = scalval
visoptions.mminval = -scalval
visoptions.autoscale    = 0
visoptions.subdivisions = 4

# ~ SetNumThreads(16)
# SetNumThreads(1)
i = 0
n = 5
# VideoStart("test.mp4")
with TaskManager(): # pajetrace=100*1000*1000):
    # while t < tend:
    if False:
        # print(t)
        w.data = -fullB.T * gfu.vec
        ### time envelope for the src ################################################################
        # ~ t_envelope = sin(2*pi*fcen*t)*exp ( -( 0.5*(t-tpeak)**2 )/(df**2))
        if abs((t-tpeak)/tpeak) < 1:
           t_envelope = (2*exp(1)/sqrt(math.pi))*sin(2*math.pi*fcen*t)*exp (-1/(1-((t-tpeak)/tpeak)**2))
        else:
           t_envelope = 0

        # ~ t_envelope = exp ( ( 0.5*(t-tpeak)**2 )/(df**2))
        w.data += t_envelope*Lsrc.vec
        ###############################################################################################
        w.data -= sigma * dampingp * gfu.vec
        w.data -= emb_p @ Bstab.mat.T * gfstab.vec

        # fes_p.SolveM (vec=w[pdofs], rho=eps)
        # fes_p.SolveM (vec=w[phatdofs], rho=eps)
        # gfu.vec.data += tau * w
        gfu.vec.data += tau * invp * w
        
        w.data = fullB * gfu.vec
        hvstab.data = Bstab.mat @ emb_p.T * gfu.vec
        w.data -= sigma * dampingu * gfu.vec
        
        # fes_u.SolveM (vec=w[udofs])
        # fes_u.SolveM (vec=w[uhatdofs])
        # gfu.vec.data += tau * w
        gfu.vec.data += tau * invu * w        
        gfstab.vec.data += tau * Mstabinv * hvstab

        t += tau
        i+=1
        if i%n == 0:
            Redraw(blocking=False)
            # SnapShot("test_" + str(i//n+1000) + ".png")
            # VideoAddFrame()
# VideoFinalize()
