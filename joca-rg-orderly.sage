# -*- mode: python -*-
#
# Discriminator: same RosenfeldGroebner(ansatz) experiment as joca-rg.sage,
# but with an ORDERLY ranking (all dependent variables + constants in a single
# block) instead of joca.sage's block-ELIMINATION ranking
# ([[DDPsi,DPsi,Psi,v,r], constants]). If RG now returns quickly, the cliff is
# the elimination ranking (a BLAD computational pathology), not the marshalling.

import time
T0 = time.time()
def stamp(msg):
    print(f"[{time.time()-T0:7.2f}s] {msg}", flush=True)

import sympy
import DifferentialAlgebra

x,y,z = sympy.var('x,y,z')
E = sympy.var('E')
v1,v2,v3,v4 = sympy.var('v1,v2,v3,v4')
a0,a1,b0,b1,c0,c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E,v1,v2,v3,v4,a0,a1,b0,b1,c0,c1]

Psi,DPsi,DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')

# ORDERLY: one block containing everything (no >> elimination between blocks)
DiffRing = DifferentialAlgebra.DifferentialRing (derivations = [x,y,z],
                                                 blocks = [[DDPsi,DPsi,Psi,v,r] + constants],
                                                 parameters = constants,
                                                 notation = 'jet')

PDE = -int(1)/int(2)*(Psi[x,x] + Psi[y,y] + Psi[z,z])*r - Psi - E*r*Psi

ansatz = [Psi[x] - DPsi * v[x],
          Psi[y] - DPsi * v[y],
          Psi[z] - DPsi * v[z],
          DPsi[x] - DDPsi * v[x],
          DPsi[y] - DDPsi * v[y],
          DPsi[z] - DDPsi * v[z],
          (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
          v - (v1*x + v2*y + v3*z + v4*r),
          r**2 - x**2 - y**2 - z**2]
ansatz = list(map(sympy.expand, ansatz))
stamp("setup done (ORDERLY ranking)")

stamp("RosenfeldGroebner(ansatz) START")
rg = DiffRing.RosenfeldGroebner(ansatz)
stamp(f"RosenfeldGroebner(ansatz) DONE: {len(rg)} component(s)")
for i, C in enumerate(rg):
    stamp(f"  component {i}: {len(C.equations())} equations")
