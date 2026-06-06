# -*- mode: python -*-
#
# joca-normalform.sage
#
# Companion to joca.sage / joca-rg.sage.  Makes the Rosenfeld–Gröbner
# SATURATION explicit.
#
# joca-rg.sage discharges hypothesis (3) by regularizing the ansatz with
# RosenfeldGroebner, then projects against the RAW ansatz.  The reason it must
# NOT instead reduce against the regularized chain's bare .equations() is the
# subject of this script: a regular differential chain C represents the
# *saturated* ideal  [C] : H_C^∞  (BLOP, "Computing representations for radicals
# of finitely generated differential ideals", Thm 28), where H_C is the product
# of the chain's initials and separants.  Reducing by the bare [C] drops the
# ":H_C^∞" and re-admits the bad locus H_C = 0 as the spurious prime (a0,a1).
#
# Here we call the chain's *saturation-aware* reduction, normal_form(), which
# returns the PDE's normal form as a FRACTION  NUM / DEN.  The denominator DEN
# is the chain's separant  S̄ = (a0+a1·v)²  reduced mod r² = x²+y²+z² — i.e. the
# ":S̄^∞" made explicit.  DEN vanishes exactly on a0 = a1 = 0, so the bad locus
# is a POLE of the saturated reduction: that is precisely why the regularized
# chain cannot see the two genuine strata #2, #3 that live inside a0 = a1 = 0.
#
# Run with:  ~/miniforge3/envs/sage/bin/sage joca-normalform.sage
# See the writeup  rg-saturation-and-the-bad-locus.pdf  for the full discussion.

import sympy
from sympy import expand, fraction, factor, together
import DifferentialAlgebra

x, y, z = sympy.var('x,y,z')
E = sympy.var('E')
v1, v2, v3, v4 = sympy.var('v1,v2,v3,v4')
a0, a1, b0, b1, c0, c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]

Psi, DPsi, DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')

DiffRing = DifferentialAlgebra.DifferentialRing(derivations=[x, y, z],
            blocks=[[DDPsi, DPsi, Psi, v, r], constants],
            parameters=constants, notation='jet')

# 2x the Schrödinger PDE (cleared of the 1/2 so normal_form's parser stays
# exact -- a constant scale does not change the redundancy locus)
PDE = -(Psi[x, x] + Psi[y, y] + Psi[z, z])*r - int(2)*Psi - int(2)*E*r*Psi

ansatz = [Psi[x] - DPsi*v[x], Psi[y] - DPsi*v[y], Psi[z] - DPsi*v[z],
          DPsi[x] - DDPsi*v[x], DPsi[y] - DDPsi*v[y], DPsi[z] - DDPsi*v[z],
          (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
          v - (v1*x + v2*y + v3*z + v4*r), r**2 - x**2 - y**2 - z**2]
ansatz = list(map(sympy.expand, ansatz))

F = DifferentialAlgebra.BaseFieldExtension(generators=constants, ring=DiffRing)
C = DiffRing.RosenfeldGroebner(ansatz, basefield=F)[0]

# ---------------------------------------------------------------------------
# Saturation-aware reduction.  normal_form returns NUM / DEN.
# ---------------------------------------------------------------------------
nf = C.normal_form(PDE)
num, den = fraction(together(nf))
num, den = expand(num), expand(den)

print("=== normal_form(PDE) = NUM / DEN ===\n")
print("DEN  (the chain separant S̄, = (a0+a1·v)² reduced mod r²=x²+y²+z²):")
print("   ", den)
print("\n   factored:", factor(den))
print("   DEN at a0=a1=0:", den.subs({a0: 0, a1: 0}),
      "  <-- the bad locus a0=a1=0 is a POLE of the saturated reduction")

# ---------------------------------------------------------------------------
# The redundancy locus is  NUM = 0  with DEN != 0.  Collect coefficients of NUM
# (a differential polynomial in Psi, DPsi, x, y, z, r) into a constant ideal.
# ---------------------------------------------------------------------------
PolyRing = PolynomialRing(QQ, names=[str(i) for i in DiffRing.indets(selection='all')])
PolyRing_constants = list(map(PolyRing, constants))

def build_system_of_equations(eqn, consts):
    ring = eqn.parent()
    d = dict()
    nc = tuple(1 if ring.gen(n) in consts else ring.gen(n) for n in range(ring.ngens()))
    for coeff, mon in eqn:
        ncp = mon(nc); cp = coeff*mon // ncp
        d[ncp] = d.get(ncp, 0) + cp
    return tuple(set(d.values()))

eqns = build_system_of_equations(PolyRing(num), PolyRing_constants)
I = ideal(eqns)

primes_num = I.minimal_associated_primes()
primes_num.sort(key=lambda p: str(p))
print("\n=== {NUM = 0}  (raw numerator ideal) -> %d primes ===" % len(primes_num))
for p in primes_num:
    print("   ", p.gens())
print("   The (a0,a1) here is exactly DEN = 0: it is the saturated reduction's")
print("   pole, NOT a redundancy cell (verify: v4*b1 is in the true redundancy")
print("   ideal and does not vanish on generic a0=a1=0).")

# Saturate away the pole DEN=0 (the bad locus) -> the generic strata only.
sat = I.saturation(ideal(PolyRing(a0), PolyRing(a1)))[0]
primes_gen = sat.minimal_associated_primes()
primes_gen.sort(key=lambda p: str(p))
print("\n=== {NUM = 0, DEN != 0}  (saturate by (a0,a1)) -> %d GENERIC strata ===" % len(primes_gen))
for p in primes_gen:
    print("   ", p.gens())
print("\nThe two strata lost to the pole (#2, #3 of joca.sage) live inside")
print("a0=a1=0, where the regularized ODE (a0+a1v)·(ODE) collapses to 0=0.  The")
print("raw projection in joca.sage keeps them because the RAW ODE degenerates to")
print("a usable 1st-order rule there instead of vanishing.  See the writeup.")
