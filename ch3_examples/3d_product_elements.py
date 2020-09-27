# Computations from subsection 3.3.2.

from firedrake import *

N = 10
base_mesh = UnitSquareMesh(N, N)
mesh = ExtrudedMesh(base_mesh, layers=N, layer_height=1./N)

P2t = FiniteElement("CG", triangle, 2)
P2i = FiniteElement("CG", interval, 2)
H1_element = TensorProductElement(P2t, P2i)
H1 = FunctionSpace(mesh, H1_element)

dP1t = FiniteElement("DG", triangle, 1)
dP1i = FiniteElement("DG", interval, 1)
L2_element = TensorProductElement(dP1t, dP1i)
L2 = FunctionSpace(mesh, L2_element)

N2_1 = FiniteElement("N2curl", triangle, 1)
P2i = FiniteElement("CG", interval, 2)
Hcurl_h = HCurl(TensorProductElement(N2_1, P2i))
Hcurl_v = HCurl(TensorProductElement(P2t, dP1i))
Hcurl_element = Hcurl_h + Hcurl_v
Hcurl = FunctionSpace(mesh, Hcurl_element)

RT2 = FiniteElement("RT", triangle, 2)
Hdiv_h = HDiv(TensorProductElement(RT2, dP1i))
Hdiv_v = HDiv(TensorProductElement(dP1t, P2i))
Hdiv_element = Hdiv_h + Hdiv_v
Hdiv = FunctionSpace(mesh, Hdiv_element)

# visualize functions from these spaces using Paraview
x, y, z = SpatialCoordinate(mesh)
uH1 = Function(H1).interpolate(x*x + y*y + z*z)
uH1.rename('uH1')
uL2 = Function(L2).interpolate(x*x + y*y + z*z)
uL2.rename('uL2')
uHcurl = Function(Hcurl).project(as_vector([z*z,y*y,x*x]))
uHcurl.rename('uHcurl')
uHdiv = Function(Hdiv).project(as_vector([z*z,y*y,x*x]))
uHdiv.rename('uHdiv')
File("unitcube.pvd").write(uH1,uL2,uHcurl,uHdiv)

