# -*- mode: python -*-
#
# joca-thomas_combined.sage --- the JOINT-VARIETY variant of joca-thomas.sage.
#
# joca-thomas.sage decomposes the ansatz A alone and then REDUCES the PDE P
# modulo each cell (Ritt) -> the MEMBERSHIP locus  V = {c : P in [A(c)]}.
#
# This script instead consumes a decomposition of the COMBINED system A u {P}
# (produced by ex4_hydrogen_combined.mpl, which adjoins the Schrodinger PDE), and
# runs NO reduction at the end.  Each cell already encodes Sol(A) cap Sol(P); we
# simply PROJECT it onto the constants -- the cell's parameter relations are the
# constant stratum -- prune by the cell's parameter inequations, and union.  By
# membership-vs-variety-partial-strata.tex this yields
#       pi_c(Sol(A) cap Sol(P))  =  V  u  {partial-solution strata}  >=  V,
# i.e. it also surfaces parameter values where a proper, lower-dimensional
# sub-family of the ansatz solves P (e.g. a normalizable bound-state mode).
#
#   sage joca-thomas_combined.sage [--log PATH] [--cells 1,2,3]

import os, sys
import sympy
try:
    import DifferentialAlgebra
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(ex.msg + "\nInstall it with '%pip install DifferentialAlgebra'")

def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default

LOG = _argval('--log', os.path.expanduser('~/thomas-experiments/hydrogen_combined.log'))
CELLS_ARG = _argval('--cells')
WANT_CELLS = set(int(n) for n in CELLS_ARG.split(',')) if CELLS_ARG else None

_here = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv and sys.argv[0] else ''
for _d in (_here, os.path.dirname(LOG)):
    if _d and _d not in sys.path:
        sys.path.insert(0, _d)
import thomas_cells as tc

# The combined system carries the energy as a parameter named EE (Maple) -> E.
tc.PARAM_NAMES = dict(tc.PARAM_NAMES); tc.PARAM_NAMES['EE'] = 'E'

# --- constants / ring (constants now include the energy E) ----------------
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
PolyRing = PolynomialRing(QQ, names=[str(indet) for indet in DiffRing.indets(selection='all')])

jetbases = {'DDPsi': DDPsi, 'DPsi': DPsi, 'Psi': Psi, 'v': v, 'rho_': r}
const_syms = {'v1': v1, 'v2': v2, 'v3': v3, 'v4': v4, 'a0': a0, 'a1': a1,
              'b0': b0, 'b1': b1, 'c0': c0, 'c1': c1, 'E': E}
parse = tc.make_translator(jetbases, const_syms, {'x': x, 'y': y, 'z': z})

with open(LOG) as fh:
    cells = tc.split_cells(fh.read())
print(f"Loaded {len(cells)} combined (A u P) Thomas cells from {LOG}\n" + "=" * 72)

# --- helpers (shared with joca-thomas.sage) -------------------------------
def prime_key(P):
    return tuple(sorted(str(g) for g in P.gens()))

def _ideal_subset(J, I):
    return all(g in I for g in J.gens())

def keep_minimal(d):
    items = list(d.items()); kept, dropped = {}, []
    for i, (ki, (Pi, ci)) in enumerate(items):
        absorber = None
        for j, (kj, (Pj, cj)) in enumerate(items):
            if i == j:
                continue
            sub_ji = _ideal_subset(Pj, Pi)
            if sub_ji and not _ideal_subset(Pi, Pj):
                absorber = kj; break
            if sub_ji and _ideal_subset(Pi, Pj) and j < i:
                absorber = kj; break
        if absorber is None:
            kept[ki] = (Pi, ci)
        else:
            dropped.append((ki, Pi, ci, absorber))
    return kept, dropped

_PSI_MAPLE_BASE = next(k for k, vv in tc.JET_NAMES.items() if vv == 'Psi')
def cell_forces_psi_zero(cell):
    target = f"{_PSI_MAPLE_BASE}(x,y,z)"
    for eq in tc._maple_list(cell['eqs_str']):
        if '=' not in eq:
            continue
        lhs, rhs = (s.replace(' ', '') for s in eq.split('=', 1))
        if rhs == '0' and lhs == target:
            return True
    return False

# --- project each combined cell onto its constant stratum (NO reduction) --
union_primes, trivial_primes = {}, {}
for cell in cells:
    num = cell['num']
    if WANT_CELLS is not None and num not in WANT_CELLS:
        continue
    cp = tc.cell_polys(cell, parse, const_syms)
    Z = sorted((p for p in cp['param_eqs']), key=str)          # constant relations of the cell
    Zideal = ideal([PolyRing(p) for p in Z]) if Z else ideal(PolyRing.zero())
    v_zero = all(PolyRing(vv) in Zideal for vv in (v1, v2, v3, v4))
    trivial = cell_forces_psi_zero(cell) or v_zero

    gens = [PolyRing(p) for p in Z]
    I = ideal(gens) if gens else ideal(PolyRing.zero())
    primes = I.minimal_associated_primes()

    ineq_polys = [PolyRing(f) for f in cp['param_ineqs']]
    survivors = []
    for P in primes:
        if P.is_one():
            continue
        if any(g in P for g in ineq_polys):
            continue
        survivors.append(P)

    tag = "TRIVIAL (Psi==0)" if trivial else "stratum"
    print(f"\n--- cell {num}: zero {{{', '.join(map(str, Z)) or '(none, generic)'}}}; "
          f"{len(cp['param_ineqs'])} param-ineqs, {len(cp['jet_ineqs'])} jet-ineqs; {tag} ---")
    if not survivors:
        print("  surviving strata: NONE (pruned by inequations / empty)")
    bucket = trivial_primes if trivial else union_primes
    for P in survivors:
        print("   V:", P, "" if not trivial else "   [trivial Psi==0]")
        bucket.setdefault(prime_key(P), (P, []))[1].append(num)

# --- union over cells -----------------------------------------------------
def dump_union(title, d):
    kept, dropped = keep_minimal(d)
    note = f"; {len(dropped)} redundant superset ideal(s) dropped" if dropped else ""
    print(f"\n{title} ({len(kept)} maximal varieties{note}):\n")
    for key, (P, cells_for) in sorted(kept.items(), key=lambda kv: str(kv[0])):
        print("  V:", P)
        print("       (from cells:", ", ".join(map(str, sorted(set(cells_for)))) + ")")
    if dropped:
        print("\n  dropped (superset of a kept ideal => smaller, redundant variety):")
        kbygens = {k: P for k, (P, _) in d.items()}
        for key, P, cells_for, absorber in sorted(dropped, key=lambda t: str(t[0])):
            print("    V:", P, "  (from cells:",
                  ", ".join(map(str, sorted(set(cells_for)))) + ")")
            print("        ⊇", kbygens[absorber])

print("\n" + "=" * 72)
scope = 'selected' if WANT_CELLS else 'all'
dump_union(f"JOINT-VARIETY solution strata over {scope} cells (membership ∪ partial)", union_primes)
print("\n" + "-" * 72)
dump_union(f"TRIVIAL (Psi==0) strata over {scope} cells", trivial_primes)
