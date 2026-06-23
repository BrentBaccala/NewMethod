# -*- mode: python -*-
#
# Verification harness for the joca.sage (NewMethod, generic) vs joca-thomas.sage
# (differential-Thomas) discrepancy.  Run with:
#
#     ~/miniforge3/envs/sage/bin/sage -python verify_discrepancy.py
#
# It reproduces and PROVES, for each of the 5+5 reported primes, its theoretical
# role: matched (good locus) / spurious-on-{h=0} / genuine-degenerate-stratum.
# See ~/project/reports/joca-thomas-vs-newmethod-discrepancy.md for the writeup.

import sympy
import DifferentialAlgebra
from sage.all import PolynomialRing, QQ

# ---------------------------------------------------------------------------
# 1. The explicit Ritt multiplier h of the GENERIC reduction (joca.sage).
# ---------------------------------------------------------------------------
x, y, z = sympy.var('x,y,z')
E = sympy.var('E')
v1, v2, v3, v4 = sympy.var('v1,v2,v3,v4')
a0, a1, b0, b1, c0, c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]
Psi, DPsi, DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')
DiffRing = DifferentialAlgebra.DifferentialRing(
    derivations=[x, y, z], blocks=[[DDPsi, DPsi, Psi, v, r], constants],
    parameters=constants, notation='jet')
PDE = -sympy.Rational(1, 2) * (Psi[x, x] + Psi[y, y] + Psi[z, z]) * r - Psi - E * r * Psi
ansatz0 = [Psi[x] - DPsi*v[x], Psi[y] - DPsi*v[y], Psi[z] - DPsi*v[z],
           DPsi[x] - DDPsi*v[x], DPsi[y] - DDPsi*v[y], DPsi[z] - DDPsi*v[z],
           (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
           v - (v1*x + v2*y + v3*z + v4*r), r**2 - x**2 - y**2 - z**2]
ansatz0 = list(map(sympy.expand, ansatz0))
h, rem = DiffRing.differential_prem(PDE, ansatz0)
print("GENERIC multiplier  h =", sympy.factor(h))
print("   => h = 64 * (a0 + a1*v) * r^2 ;  initial factor (a0+a1*v), separant factor r^2.\n")

# ---------------------------------------------------------------------------
# 2. Classify each prime by whether the ANSATZ-substituted initial (a0+a1*v),
#    v = v1*x+v2*y+v3*z+v4*r, vanishes IDENTICALLY in space (=> on {h=0}).
#    The initial polynomial's coeffs in 1,x,y,z,r are a0, a1*v1, a1*v2, a1*v3, a1*v4.
# ---------------------------------------------------------------------------
V = ['x', 'y', 'z', 'DDPsi', 'DPsi', 'Psi', 'v', 'r', 'E',
     'v1', 'v2', 'v3', 'v4', 'a0', 'a1', 'b0', 'b1', 'c0', 'c1']
R = PolynomialRing(QQ, V)
def Idl(*s): return R.ideal([R(e) for e in s])
NM = {'NM1': Idl('a0', 'v4', 'v3', 'v2', 'v1'),
      'NM2': Idl('c1', 'c0', 'b1', 'b0', 'a1', 'a0'),
      'NM3': Idl('v1^2+v2^2+v3^2', 'a1', 'a0', 'v4'),
      'NM4': Idl('v4*c0-b0', 'E*c0-v4*c1', '-v4^2*c1+E*b0', 'b1', 'a1-1/2*b0', 'a0', 'v3', 'v2', 'v1'),
      'NM5': Idl('v4*c0-b0', 'v1^2+v2^2+v3^2-v4^2', 'c1', 'b1', 'a1-b0', 'a0', 'E')}
TH = {'TH_A': Idl('-b1*c0+b0*c1', 'v4*c1-b1', 'v4*c0-b0', 'a1', 'a0', 'v3', 'v2', 'v1', 'E+1/2'),
      'TH_B': Idl('v4*c0-b0', 'E*c0-v4*c1', '-v4^2*c1+E*b0', 'b1', 'a1-1/2*b0', 'a0', 'v3', 'v2', 'v1'),
      'TH_C': Idl('c0^2+4*b0*c1', '4*v4*c1+c0', 'v4*c0-b0', 'b1+1/2*c0', 'a1', 'a0', 'v3', 'v2', 'v1', 'E+1/8'),
      'TH_D': Idl('v4*c0-b0', 'v1^2+v2^2+v3^2-v4^2', 'c1', 'b1', 'a1-b0', 'a0', 'E'),
      'TH_E': Idl('b0', 'a1', 'a0', 'v4', 'v3', 'v2', 'v1')}
init_coeffs = ['a0', 'a1*v1', 'a1*v2', 'a1*v3', 'a1*v4']
print("On-{h=0} test (initial a0+a1*v vanishes identically on the stratum):")
for nm, J in list(NM.items()) + list(TH.items()):
    on = all(R(p) in J for p in init_coeffs)
    print(f"   {nm}: {on}")
print()

# ---------------------------------------------------------------------------
# 3. The BRIDGE: saturate the generic ideal (= intersection of NM_i) by the
#    multiplier's parameter shadow.  Removes every {h=0} component.
# ---------------------------------------------------------------------------
J = NM['NM1']
for K in [NM['NM2'], NM['NM3'], NM['NM4'], NM['NM5']]:
    J = J.intersection(K)
Hinit = Idl('a0', 'a1*v1', 'a1*v2', 'a1*v3', 'a1*v4')
Jsat = J.saturation(Hinit)[0]
print("Generic ideal saturated by (initial)^inf leaves:")
for P in Jsat.minimal_associated_primes():
    print("   ", list(P.gens()))
print("   => exactly NM4=TH_B and NM5=TH_D (the GOOD-locus solutions).\n")

# ---------------------------------------------------------------------------
# 4. The degenerate radial stratum (cell 26: a0=a1=v1=v2=v3=0) re-reduction,
#    which the generic pass multiplies away.  Produces TH_A (E=-1/2), TH_C (E=-1/8).
# ---------------------------------------------------------------------------
sub = {a0: 0, a1: 0, v1: 0, v2: 0, v3: 0}
spec = [e for e in (sympy.expand(a.subs(sub)) for a in ansatz0) if e != 0]
hd, remd = DiffRing.differential_prem(PDE, spec)
print("Degenerate radial stratum multiplier h =", sympy.factor(hd),
      "  (leader dropped Psi'' -> Psi'; new initial b0+b1*v)")
PolyRing = PolynomialRing(QQ, names=[str(i) for i in DiffRing.indets(selection='all')])
PolyRing_constants = list(map(PolyRing, constants))
def build_system(eqn, consts):
    ring = eqn.parent(); sysd = {}
    ncs = tuple(1 if ring.gen(n) in consts else ring.gen(n) for n in range(ring.ngens()))
    for coeff, mon in eqn:
        ncp = mon(ncs); cp = coeff*mon // ncp
        sysd[ncp] = sysd.get(ncp, 0) + cp
    return tuple(set(sysd.values()))
eqns = build_system(PolyRing(remd), PolyRing_constants)
gens = list(eqns) + [PolyRing(c) for c in [a0, a1, v1, v2, v3]]
primesd = PolyRing.ideal(gens).minimal_associated_primes()
print("Degenerate stratum minimal primes (the strata joca.sage drops):")
for P in primesd:
    Eg = [g for g in P.gens() if 'E' in str(g)]
    print("   E-cond", Eg, "|", list(P.gens()))
