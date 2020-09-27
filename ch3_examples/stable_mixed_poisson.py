# Computations from section 3.1.  Generate plot similar to Fig. 3.1 (right).

from firedrake import *

mesh = UnitSquareMesh(16, 16)

U = FunctionSpace(mesh, "RT", 1)
V = FunctionSpace(mesh, "DG", 0)
W = U * V

sigma, u = TrialFunctions(W)
tau, v = TestFunctions(W)
f = Function(V)

x, y = SpatialCoordinate(mesh)
f.interpolate(-2*(x-1)*x - 2*(y-1)*y)

a = (dot(tau, sigma) - u*div(tau) + div(sigma)*v)*dx
L = f*v*dx

w = Function(W, name="Solution")
solve(a == L, w)

# visualize
import matplotlib.pyplot as plt
sigma_h, u_h = w.split()
pic = tripcolor(u_h, cmap='coolwarm')
plt.axis('equal')
plt.colorbar(pic)
plt.show()

