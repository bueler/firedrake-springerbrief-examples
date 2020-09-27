# Computations from section 3.2.  Use Paraview to generate something like Fig. 3.2.

from firedrake import *

base_mesh = CubedSphereMesh(radius=6400.e6, refinement_level=5,
                            degree=2)
mesh = ExtrudedMesh(base_mesh, layers=20, layer_height=10,
                    extrusion_type="radial")

# visualize a function from this space with Paraview
V = FunctionSpace(mesh, "CG", 1)
u = Function(V)
File("cubedsphere.pvd").write(u)

