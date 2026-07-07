# -*- mode: python -*-
#
# joca-prolong2.sage -- the bounded-prolongation ("adjoin-then-eliminate",
# existence) computation for the hydrogen example.
#
# Author: Brent Baccala
# Date: July 6, 2026
#
# Tested on Ubuntu 24 with Sage 10.7.  Pure Sage: unlike joca.sage and
# joca-rg.sage, NO differential-algebra package is used -- the only
# differential operation is a mechanical prolongation on jet variables.
#
# Implements the algorithm of bounded-prolongation-theorem.md (SS2bis/SS2ter,
# single-ideal route):
#
#   1. PROLONG.  Differentiate every equation of the combined system A u {P}
#      (ansatz + PDE together, on equal footing) up to total order N
#      (--order, default N = 2, the conjectured involutive order N0).
#   2. ALGEBRAIZE.  Treat every jet as an independent polynomial
#      indeterminate, giving an ideal J in K[jets, x,y,z, constants].
#      The monic-linear equations (chain rules, the v-relation, and their
#      prolongations) are substituted away first -- an exact triangular
#      elimination that inverts no initials, so the bad locus is untouched.
#      Where two different prolongations define the same jet (e.g. Psi_xy
#      via d/dy of the x-chain-rule or d/dx of the y-chain-rule), the
#      difference is kept as a generator: these are exactly the
#      integrability (Delta-polynomial) conditions, surfaced mechanically.
#      NO saturation by the initials/separants H_A is performed.
#   3. NONTRIVIALITY.  Every generator is linear homogeneous in the jet
#      block (Psi, DPsi, DDPsi and their jets), so a fiber point with
#      Psi != 0 can be rescaled to Psi = 1: dehomogenizing by Psi = 1
#      effects the quotient J : Psi^oo exactly, without the GB-hostile
#      Rabinowitsch generator t*Psi - 1.  This deletes only the zero
#      section Psi == 0 -- without it the elimination below is vacuous
#      (J n K[c] = (0), see "Two saturations" in the note).  It is NOT
#      the H_A saturation.
#   4. ELIMINATE + DECOMPOSE.  J_c = (J : Psi^oo) n K[c] by block-order
#      elimination of all non-constant variables; the minimal associated
#      primes of J_c are the irreducible components of (the Zariski closure
#      of) the EXISTENCE locus -- the constants at which the ansatz family
#      contains a nonzero solution of the PDE -- bad-locus strata included.
#
# Contrast joca.sage (reduce-then-project): that computes the MEMBERSHIP
# locus {c : P in [A(c)]} off the denominator locus h != 0, and so cannot
# certify the a0 = a1 = 0 strata.  This script's existence locus keeps them:
# expect, in addition to the classical and E = 0 branches, the
# first-order-ODE strata at E = -1/2 (ground state e^-r) and E = -1/8 (the
# 2s state (1 - r/2) e^{-r/2}) -- compare cells 26/27 of the differential
# Thomas run, joca-thomas-openmaple.out.
#
# PROLONG-UNTIL-STABLE.  N = 2 is the involutive order of the ANSATZ alone
# (mixed partials of v commute); it is NOT the involutive order of the
# combined system A u {P}.  The PDE and the ODE share the leader Psi'' and
# form a critical pair whose resolution is the Ritt remainder rho; the
# conditions that cut the constants are rho and its total derivatives
# D^k rho, and D^k rho enters the prolonged ideal only through D^k P,
# i.e. at N = 2 + k.  At N = 2 only rho itself (one pointwise condition)
# is present and the existence locus is far too big.  Under-prolonging
# errs safe -- a superset: spurious constants possible, solutions never
# lost.  Run increasing --order and watch J_c grow; stabilization pins
# N0(A u {P}) empirically (Noetherian termination, see the note SS2bis).
#
# Options:
#   --order N    prolongation order (default 2; must be >= 2 = ord P)
#   --char P     compute over GF(P) instead of QQ (fast reconnaissance,
#                e.g. --char 32003; primes valid with high probability)

import sys
import time

def opt(flag, default, cast):
    if flag in sys.argv:
        return cast(sys.argv[sys.argv.index(flag) + 1])
    return default

N = opt('--order', 2, int)
char = opt('--char', 0, int)
assert N >= 2, "--order must be >= 2 (the PDE has order 2)"
K = QQ if char == 0 else GF(char)

t0 = time.time()
def stage(msg):
    print(f"[{time.time()-t0:8.1f}s] {msg}", flush=True)

print(f"bounded-prolongation existence computation: order N = {N}, field = {K}")

# ----------------------------------------------------------------- jet ring
#
# Dependent symbols and their jets to order N.  Psi', Psi'' are the separate
# indeterminates DPsi, DDPsi exactly as in joca.sage; each of the five
# symbols gets the full set of spatial jets (DDPsi has no chain rule, so its
# jets are genuine unknowns -- their only constraints come from prolonging
# the ODE, which is where the bad-locus branching lives).

deps = ['Psi', 'DPsi', 'DDPsi', 'v', 'r']
const_names = ['E', 'v1', 'v2', 'v3', 'v4', 'a0', 'a1', 'b0', 'b1', 'c0', 'c1']

def multi_indices(n):
    """All (i,j,k) with i+j+k <= n, ascending total order."""
    return [(i, j, t - i - j)
            for t in range(n + 1)
            for i in range(t + 1)
            for j in range(t - i + 1)]

midx = multi_indices(N)

def jetname(w, m):
    s = 'x' * m[0] + 'y' * m[1] + 'z' * m[2]
    return w if not s else w + '_' + s

jet_list = [(w, m) for w in deps for m in midx]
names = ([jetname(w, m) for (w, m) in jet_list]
         + ['x', 'y', 'z']
         + const_names)
n_nonconst = len(jet_list) + 3
R = PolynomialRing(K, names,
                   order=TermOrder('degrevlex', n_nonconst)
                       + TermOrder('degrevlex', len(const_names)))
g = dict(zip(names, R.gens()))
jet = {(w, m): g[jetname(w, m)] for (w, m) in jet_list}

X, Y, Z = g['x'], g['y'], g['z']
E = g['E']
v1, v2, v3, v4 = g['v1'], g['v2'], g['v3'], g['v4']
a0, a1, b0, b1, c0, c1 = (g[s] for s in ['a0', 'a1', 'b0', 'b1', 'c0', 'c1'])
Psi, DPsi, DDPsi, v, r = (g[s] for s in deps)

# Total-derivative tables: succ[w][d] = D_d(w) for every ring generator.

succ = {}
for (w, m) in jet_list:
    succ[jet[(w, m)]] = [jet.get((w, tuple(m[i] + (1 if i == d else 0)
                                           for i in range(3))))
                         for d in range(3)]
succ[X] = [R.one(), R.zero(), R.zero()]
succ[Y] = [R.zero(), R.one(), R.zero()]
succ[Z] = [R.zero(), R.zero(), R.one()]
for s in const_names:
    succ[g[s]] = [R.zero()] * 3

def Dtot(poly, d):
    """Total derivative of a jet-ring polynomial in direction d."""
    out = R.zero()
    for w in poly.variables():
        s = succ[w][d]
        if s is None:
            raise ValueError(f"jet order overflow differentiating {w}")
        if s:
            out += poly.derivative(w) * s
    return out

def prolongations(eq, maxord):
    """{multi-index m: D^m eq} for |m| <= maxord (derivations commute)."""
    memo = {(0, 0, 0): eq}
    for m in multi_indices(maxord):
        if sum(m) == 0:
            continue
        d = next(i for i in range(3) if m[i] > 0)
        m0 = tuple(m[i] - (1 if i == d else 0) for i in range(3))
        memo[m] = Dtot(memo[m0], d)
    return memo

# ------------------------------------------- the combined system A u {P}

def e(w, *m):
    return jet[(w, tuple(m))]

# Monic-linear elements: (equation, leader).  These and their prolongations
# become exact substitution rules, not ideal generators.
solved = [
    (e('Psi', 1, 0, 0) - DPsi * e('v', 1, 0, 0), ('Psi', (1, 0, 0)), 1),
    (e('Psi', 0, 1, 0) - DPsi * e('v', 0, 1, 0), ('Psi', (0, 1, 0)), 1),
    (e('Psi', 0, 0, 1) - DPsi * e('v', 0, 0, 1), ('Psi', (0, 0, 1)), 1),
    (e('DPsi', 1, 0, 0) - DDPsi * e('v', 1, 0, 0), ('DPsi', (1, 0, 0)), 1),
    (e('DPsi', 0, 1, 0) - DDPsi * e('v', 0, 1, 0), ('DPsi', (0, 1, 0)), 1),
    (e('DPsi', 0, 0, 1) - DDPsi * e('v', 0, 0, 1), ('DPsi', (0, 0, 1)), 1),
    (v - (v1*X + v2*Y + v3*Z + v4*r), ('v', (0, 0, 0)), 0),
]

# Kept elements: (equation, differentiation order already present).  The
# PDE is the paper's, cleared of the 1/2 by multiplying through by -2 (a
# unit -- same ideal): -1/2 (Psi_xx+Psi_yy+Psi_zz) r - Psi - E r Psi = 0.
ODE = (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi
r_rel = r**2 - X**2 - Y**2 - Z**2
PDE = (e('Psi', 2, 0, 0) + e('Psi', 0, 2, 0) + e('Psi', 0, 0, 2))*r \
      + 2*Psi + 2*E*r*Psi
kept = [(ODE, 0), (r_rel, 0), (PDE, 2)]

# --------------------------------- step 1+2: prolong, collect rules + gens

rules = {}
gens = []
integrability = 0

for eq, (w, m0), o in solved:
    for m, peq in prolongations(eq, N - o).items():
        L = jet[(w, tuple(m0[i] + m[i] for i in range(3)))]
        rhs = L - peq                       # peq = L - F, so rhs = F
        if L in rules:
            # Same jet reached by two prolongation paths: keep the
            # difference -- an integrability (Delta-polynomial) condition.
            diff = rules[L] - rhs
            if diff:
                gens.append(diff)
                integrability += 1
        else:
            rules[L] = rhs

for eq, o in kept:
    gens.extend(prolongations(eq, N - o).values())

stage(f"prolonged: {len(rules)} substitution rules, {len(gens)} generators "
      f"({integrability} nonzero integrability differences)")

# Resolve the rules against each other (the dependency graph is acyclic:
# Psi-jets -> DPsi-jets -> DDPsi/v-jets -> r-jets), then substitute into
# the kept generators.  This is an exact elimination of monic-linear
# variables: no initial is inverted, the bad locus is untouched.

keys = set(rules)
changed = True
while changed:
    changed = False
    for L in rules:
        pending = [u for u in rules[L].variables() if u in keys]
        if pending:
            rules[L] = rules[L].subs({u: rules[u] for u in pending})
            changed = True

def substitute(q):
    pending = [u for u in q.variables() if u in keys]
    return q.subs({u: rules[u] for u in pending}) if pending else q

gens = [substitute(q) for q in gens]
gens = [q for q in gens if q]
assert all(u not in keys for q in gens for u in q.variables())

used_vars = sorted({str(u) for q in gens for u in q.variables()})
stage(f"substituted: {len(gens)} generators in {len(used_vars)} variables")

# ------------------------- step 3+4: nontriviality, eliminate, decompose
#
# Each generator is linear homogeneous in the jet block (the r-relation and
# its prolongations have jet-degree 0), so scaling a Psi != 0 fiber point by
# 1/Psi normalizes Psi = 1: dehomogenize instead of saturating.

jetblock = set(jet.values())
for q in gens:
    degs = {sum(mo.degree(u) for u in mo.variables() if u in jetblock)
            for mo in q.monomials()}
    assert len(degs) == 1, f"not jet-homogeneous: {q}"

gens = [q.subs({Psi: R.one()}) for q in gens]

J = R.ideal(gens)
elim = [R.gen(i) for i in range(n_nonconst)]
stage(f"eliminating {n_nonconst} non-constant variables ...")
Jc = J.elimination_ideal(elim)
stage(f"elimination ideal J_c: {len(Jc.gens())} generators")

Rc = PolynomialRing(K, const_names, order='degrevlex')
Jc = Rc.ideal([Rc(str(q)) for q in Jc.gens()])

gens_str = [str(q) for q in Jc.gens()]
if sum(len(s) for s in gens_str) < 4000:
    print("\nJ_c generators:", *gens_str, sep='\n')
else:
    print(f"\nJ_c generators: {len(gens_str)} (large; printing degrees only)")
    print(sorted(q.degree() for q in Jc.gens()))

stage("computing minimal associated primes of J_c ...")
primes = Jc.minimal_associated_primes()
primes.sort(key=str)

print(f"\nMinimal associated primes of J_c "
      f"(components of the existence locus, order {N}):")
for P in primes:
    print(tuple(P.gens()))

print("\n(Existence semantics: constants where the ansatz family contains a"
      "\nnonzero PDE solution.  Compare joca.sage's five membership primes and"
      "\nthe nontrivial solution varieties of joca-thomas-openmaple.out; rerun"
      f"\nwith --order {N+1} to check prolong-until-stable.)")
