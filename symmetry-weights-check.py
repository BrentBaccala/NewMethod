#!/usr/bin/env python3
"""Verify the Tier-1 (scaling) weight lattice of the hydrogen / ansatz 5 system.

A weight is a vector in the free lattice with coordinates (Theta, mu, P):
Theta = the ODE-multiplier direction, mu = the scaling of the inner variable v,
P = the scaling of Psi.  A generator is admissible iff all of its monomials
carry the same weight.  Run:  python3 symmetry-weights-check.py
"""
from itertools import combinations
from fractions import Fraction

W = {                       # (Theta, mu, P)
    'x': (0,0,0), 'y': (0,0,0), 'z': (0,0,0), 'r': (0,0,0),
    'Psi': (0,0,1), 'Psi_x': (0,0,1), 'Psi_xx': (0,0,1),
    'DPsi': (0,-1,1), 'DPsi_x': (0,-1,1),
    'DDPsi': (0,-2,1),
    'v': (0,1,0), 'v_x': (0,1,0),
    'v1': (0,1,0), 'v2': (0,1,0), 'v3': (0,1,0), 'v4': (0,1,0),
    'a0': (1,2,0), 'a1': (1,1,0),
    'b0': (1,1,0), 'b1': (1,0,0),
    'c0': (1,0,0), 'c1': (1,-1,0),
    'E':  (0,0,0),
}
PARAMS = ['v1','v2','v3','v4','a0','a1','b0','b1','c0','c1','E']

def wt(mon):
    t = [0,0,0]
    for s in mon:
        for i in range(3):
            t[i] += W[s][i]
    return tuple(t)

# each generator as a list of monomials (a monomial = list of symbols)
GENS = {
 'ansatz 1  Psi[x] - DPsi*v[x]':      [['Psi_x'], ['DPsi','v_x']],
 'ansatz 4  DPsi[x] - DDPsi*v[x]':    [['DPsi_x'], ['DDPsi','v_x']],
 'ansatz 7  (a0+a1 v)DDPsi + (b0+b1 v)DPsi + (c0+c1 v)Psi':
     [['a0','DDPsi'], ['a1','v','DDPsi'],
      ['b0','DPsi'],  ['b1','v','DPsi'],
      ['c0','Psi'],   ['c1','v','Psi']],
 'ansatz 8  v - (v1 x + v2 y + v3 z + v4 r)':
     [['v'], ['v1','x'], ['v2','y'], ['v3','z'], ['v4','r']],
 'ansatz 9  r^2 - (x^2+y^2+z^2)':
     [['r','r'], ['x','x'], ['y','y'], ['z','z']],
 'PDE  -Psi[x,x] r - 2 Psi r E - 2 Psi':
     [['Psi_xx','r'], ['Psi','r','E'], ['Psi']],
 'inequation v':   [['v']],
 'inequation Psi': [['Psi']],
}

print("=== generator-by-generator quasi-homogeneity ===")
ok = True
for name, mons in GENS.items():
    ws = {wt(m) for m in mons}
    good = len(ws) == 1
    ok &= good
    print(f"  [{'OK ' if good else 'BAD'}] {name}\n         weights: {sorted(ws)}")

print(f"\nAll generators quasi-homogeneous: {ok}")

print("\n=== parameter weight table (Theta, mu) ===")
for p in PARAMS:
    print(f"  {p:>3}: {W[p][:2]}")

def det2(u, v):
    return u[0]*v[1] - u[1]*v[0]

print("\n=== which 2-parameter pins are clean? (|det| == 1) ===")
clean, cover, degen = [], [], []
for p, q in combinations(PARAMS, 2):
    d = det2(W[p][:2], W[q][:2])
    if d == 0:      degen.append((p,q))
    elif abs(d)==1: clean.append((p,q))
    else:           cover.append((p,q,abs(d)))
print(f"  clean (unimodular): {len(clean)}, e.g. {clean[:6]}")
print(f"  DEFECTIVE d-fold covers ({len(cover)}):")
for p,q,d in cover:
    print(f"      {p} = {q} = 1   ->  {d}-fold cover")
print(f"  degenerate (rank<2, not a valid pin): {len(degen)}")
print(f"      involving E: {[c for c in degen if 'E' in c]}")
