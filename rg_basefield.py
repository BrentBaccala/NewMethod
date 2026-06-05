# RG experiment with the constants moved INTO the coefficient field
# (BaseFieldExtension over Q(E,v1,...,c1)), to see if it tames the
# bad_remainder_irreducible_factorwise / baz_Yun factorization blowup.
import time
T0 = time.time()
def stamp(msg):
    print(f"[{time.time()-T0:7.2f}s] {msg}", flush=True)

import sympy
import DifferentialAlgebra
from DifferentialAlgebra import BaseFieldExtension

x, y, z = sympy.var('x,y,z')
E = sympy.var('E')
v1, v2, v3, v4 = sympy.var('v1,v2,v3,v4')
a0, a1, b0, b1, c0, c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]

Psi, DPsi, DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')

DiffRing = DifferentialAlgebra.DifferentialRing(
    derivations=[x, y, z],
    blocks=[[DDPsi, DPsi, Psi, v, r], constants],
    parameters=constants,
    notation='jet')

# Move the constants into the base field: Q(E,v1,...,c1)
F = BaseFieldExtension(generators=constants, ring=DiffRing)
stamp("base field extension built over the constants")

ansatz = [Psi[x] - DPsi * v[x],
          Psi[y] - DPsi * v[y],
          Psi[z] - DPsi * v[z],
          DPsi[x] - DDPsi * v[x],
          DPsi[y] - DDPsi * v[y],
          DPsi[z] - DDPsi * v[z],
          (a0 + a1 * v) * DDPsi + (b0 + b1 * v) * DPsi + (c0 + c1 * v) * Psi,
          v - (v1 * x + v2 * y + v3 * z + v4 * r),
          r**2 - x**2 - y**2 - z**2]
ansatz = list(map(sympy.expand, ansatz))
stamp("setup done; RosenfeldGroebner(ansatz, basefield=F) START")

rg = DiffRing.RosenfeldGroebner(ansatz, basefield=F)
stamp(f"RosenfeldGroebner DONE: {len(rg)} component(s)")
for i, C in enumerate(rg):
    stamp(f"  component {i}: {len(C.equations())} equations")
