# -*- mode: python -*-
#
# Sage script for the paper in the Journal of Computational Algebra.
#
# This is the differential-Thomas variant of joca-rg.sage.  Instead of running
# the paper's algorithm once against the single (generic) ansatz, it loads the
# disjoint simple systems ("cells") of the differential Thomas decomposition of
# the hydrogen ansatz -- computed by Maple's DifferentialThomas package and
# logged in ~/thomas-experiments/hydrogen_thomas.log -- and runs the algorithm
# once per cell:
#
#     for each Thomas cell:
#         specialize the ansatz to the cell's parameter stratum
#         reduce the PDE modulo that ansatz        (Ritt full reduction)
#         gather the remainder's coefficients in the parameters   (projection)
#         compute minimal associated primes        (GTZ)
#         drop primes excluded by the cell's inequations          (saturation)
#     union the surviving primes over all cells
#
# Each cell's parameter stratum (verified empirically across all 73 cells) is a
# subset of the constants {a0,a1,b0,b1,c0,c1,v1,v2,v3,v4} set to zero; so the
# cell-specialized ansatz is just the original ansatz with those constants set
# to 0 (the chain-rule relations survive intact, so Ritt reduction reduces
# fully).  Because the cells are *disjoint* and each is a squarefree regular
# simple system, the per-cell saturation is valid pointwise -- and the union
# over the cells is the complete decomposition.  This is the paper's
# completeness-by-composition thesis made executable: differential Thomas as
# the completeness oracle, the basic algorithm run cell-by-cell.
#
# Author: Brent Baccala
# Date: June 2026
# Tested on Ubuntu 24 with Sage 10.7 and DifferentialAlgebra 1.14.
#
#   sage joca-thomas.sage [--log PATH] [--cells 1,45,46] [--verbose-remainder]

import os, sys

import sympy

try:
    import DifferentialAlgebra
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(ex.msg + "\nInstall it with '%pip install DifferentialAlgebra'")

# --- command-line options -------------------------------------------------
def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default

LOG = _argval('--log', os.path.expanduser('~/thomas-experiments/hydrogen_thomas.log'))
CELLS_ARG = _argval('--cells')                       # e.g. "1,45,46"; default all
WANT_CELLS = set(int(n) for n in CELLS_ARG.split(',')) if CELLS_ARG else None
VERBOSE_REM = '--verbose-remainder' in sys.argv

# the cell parser (thomas_cells.py) ships next to this script; also look in the
# log's directory so the parser can live with the log if preferred
_here = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv and sys.argv[0] else ''
for _d in (_here, os.path.dirname(LOG)):
    if _d and _d not in sys.path:
        sys.path.insert(0, _d)
import thomas_cells as tc

# Claude Sonnet 4's solution to printing "DPsi" as "\Psi'" in LaTeX
def patch_latex_varify():
    from sage.misc.latex import latex_varify
    import sage.misc.latex
    original = latex_varify
    def custom(a, is_fname=False):
        return r"\Psi'" if a == "DPsi" else original(a, is_fname=is_fname)
    sage.misc.latex.latex_varify = custom
patch_latex_varify()

# --- ring / PDE / ansatz (identical to joca.sage) -------------------------
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
                                                parameters=constants,
                                                notation='jet')

PDE = -int(1)/int(2)*(Psi[x, x] + Psi[y, y] + Psi[z, z])*r - Psi - E*r*Psi
print("PDE:", PDE)

ansatz0 = [Psi[x] - DPsi*v[x],
           Psi[y] - DPsi*v[y],
           Psi[z] - DPsi*v[z],
           DPsi[x] - DDPsi*v[x],
           DPsi[y] - DDPsi*v[y],
           DPsi[z] - DDPsi*v[z],
           (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
           v - (v1*x + v2*y + v3*z + v4*r),
           r**2 - x**2 - y**2 - z**2]
ansatz0 = list(map(sympy.expand, ansatz0))

# --- Sage polynomial ring (identical to joca.sage) ------------------------
PolyRing = PolynomialRing(QQ, names=[str(indet) for indet in DiffRing.indets(selection='all')])
PolyRing_constants = list(map(PolyRing, constants))


def build_system_of_equations(eqn, constants):
    """Factor each term into constant * non-constant parts, gather like
    non-constant monomials; return the (constant-only) coefficient polynomials.
    Identical to joca.sage / joca-rg.sage."""
    ring = eqn.parent()
    system = dict()
    non_constant_sub = tuple(1 if ring.gen(n) in constants else ring.gen(n)
                             for n in range(ring.ngens()))
    for coeff, monomial in eqn:
        non_constant_part = monomial(non_constant_sub)
        constant_part = coeff * monomial // non_constant_part
        if non_constant_part in system:
            system[non_constant_part] += constant_part
        else:
            system[non_constant_part] = constant_part
    return tuple(set(system.values()))


# --- the differential-Thomas cells ----------------------------------------
jetbases = {'DDPsi': DDPsi, 'DPsi': DPsi, 'Psi': Psi, 'v': v, 'rho_': r}
const_syms = {'v1': v1, 'v2': v2, 'v3': v3, 'v4': v4,
              'a0': a0, 'a1': a1, 'b0': b0, 'b1': b1, 'c0': c0, 'c1': c1}
parse = tc.make_translator(jetbases, const_syms, {'x': x, 'y': y, 'z': z})

with open(LOG) as fh:
    cells = tc.split_cells(fh.read())
print(f"\nLoaded {len(cells)} Thomas cells from {LOG}\n" + "=" * 72)


def prime_key(P):
    """A canonical, sortable string for a prime ideal (for dedup / union)."""
    return tuple(sorted(str(g) for g in P.gens()))


def _ideal_subset(J, I):
    """True iff ideal J ⊆ ideal I (every generator of J lies in I)."""
    return all(g in I for g in J.gens())


def keep_minimal(d):
    """Keep only the ⊆-minimal ideals in a {key: (ideal, cells)} dict.

    We report the largest varieties, i.e. the smallest ideals: if I_a ⊊ I_b then
    V(I_a) ⊋ V(I_b), so V(I_b) is already contained in V(I_a) and is redundant in
    the union.  Drop every ideal that strictly contains another in the list (the
    superset ideals); the survivors are each presented in their simplest form.
    Returns (kept_dict, dropped) with dropped = [(key, P, cells, absorber_key)].
    """
    items = list(d.items())
    kept, dropped = {}, []
    for i, (ki, (Pi, ci)) in enumerate(items):
        absorber = None
        for j, (kj, (Pj, cj)) in enumerate(items):
            if i == j:
                continue
            sub_ji = _ideal_subset(Pj, Pi)
            if sub_ji and not _ideal_subset(Pi, Pj):          # Pj ⊊ Pi  -> drop Pi
                absorber = kj; break
            if sub_ji and _ideal_subset(Pi, Pj) and j < i:    # exact dup -> keep earlier
                absorber = kj; break
        if absorber is None:
            kept[ki] = (Pi, ci)
        else:
            dropped.append((ki, Pi, ci, absorber))
    return kept, dropped


# Cache reduction+decomposition by parameter stratum (cells that zero the same
# constants share an identical specialized ansatz, hence identical reduction and
# primes; only their inequations -- and thus the pruning -- differ).
strata_cache = {}
union_primes = {}            # prime_key -> (ideal, contributing cells)  -- nontrivial
trivial_primes = {}          # prime_key -> (ideal, contributing cells)  -- Psi==0 strata

# The Maple base name of the order-0 wavefunction jet (Ps -> Psi).
_PSI_MAPLE_BASE = next(k for k, v in tc.JET_NAMES.items() if v == 'Psi')

def cell_forces_psi_zero(cell):
    """True iff the cell's defining equations include  <Psi>(x,y,z) = 0 -- i.e.
    DifferentialThomas ITSELF forced the wavefunction's order-0 jet to vanish (a
    Psi==0 stratum).  This is authoritative and immune to the Ritt-multiplier
    artifact that fools the reduction-based test below: on a stratum where v==0
    the PDE's multiplier h = 64*(a0+a1*v)*r^2 vanishes, so differential_prem
    returns a false remainder==0 AND the Psi-vs-ansatz reduction keeps Psi a live
    indeterminate -- yet the PDE forces Psi=0 (it reads -(1+E*r)*Psi there).
    Cell 14 (the spurious prime TH_E) is exactly this case; Thomas records it as
    `Ps(x,y,z) = 0`, which we honor here."""
    target = f"{_PSI_MAPLE_BASE}(x,y,z)"
    for eq in tc._maple_list(cell['eqs_str']):
        if '=' not in eq:
            continue
        lhs, rhs = (s.replace(' ', '') for s in eq.split('=', 1))
        if rhs == '0' and lhs == target:
            return True
    return False

for cell in cells:
    num = cell['num']
    if WANT_CELLS is not None and num not in WANT_CELLS:
        continue
    cp = tc.cell_polys(cell, parse, const_syms)
    Z = sorted((p for p in cp['param_eqs']), key=str)        # constants set to 0
    Zkey = tuple(map(str, Z))

    if Zkey not in strata_cache:
        sub = {p: 0 for p in Z}
        spec = [e for e in (sympy.expand(a.subs(sub)) for a in ansatz0) if e != 0]
        h, rem = DiffRing.differential_prem(PDE, spec)
        # Does the specialized ansatz force the wavefunction to vanish?  If Psi
        # reduces to 0 modulo the stratum's ansatz, the only "solution" there is
        # the trivial Psi == 0 -- a degenerate stratum the paper's generic
        # reduction never reports (it keeps Psi a live indeterminate).
        _, psi_rem = DiffRing.differential_prem(Psi, spec)
        # v ≡ 0 stratum: if all of v1..v4 vanish, the ansatz forces Psi spatially
        # constant (Psi_x = DPsi*v_x = 0, ...), and the PDE on a constant Psi reads
        # -(1+E*r)*Psi, which has NO nonzero solution -> Psi ≡ 0.  This is the
        # degenerate locus where the Ritt multiplier h = 64*(a0+a1*v)*r^2 vanishes
        # through v=0, so differential_prem returns a false remainder==0 with Psi
        # kept live; the reduction-based test below cannot see it.  (Cells 14/15/16
        # are this case -> the spurious primes TH_E and its sub-strata.)
        Zideal = ideal([PolyRing(p) for p in Z]) if Z else ideal(PolyRing.zero())
        v_zero = all(PolyRing(vv) in Zideal for vv in (v1, v2, v3, v4))
        trivial = (psi_rem == 0) or v_zero
        if rem == 0:
            eqns = ()
        else:
            eqns = build_system_of_equations(PolyRing(rem), PolyRing_constants)
        # ideal = projected coefficient conditions  +  the stratum (Z = 0)
        gens = list(eqns) + [PolyRing(p) for p in Z]
        I = ideal(gens) if gens else ideal(PolyRing.zero())
        primes = I.minimal_associated_primes()
        strata_cache[Zkey] = dict(spec_len=len(spec), rem=rem, eqns=eqns,
                                  primes=primes, trivial=trivial)

    sc = strata_cache[Zkey]
    # Per-cell triviality: the cached reduction-based test (per parameter stratum)
    # OR the cell's own Ps(x,y,z)=0 equation.  The latter catches the v==0 strata
    # where the multiplier artifact makes the reduction test wrongly report the
    # stratum as a nontrivial solution (e.g. cell 14 -> the spurious TH_E).
    cell_trivial = sc['trivial'] or cell_forces_psi_zero(cell)

    # cell inequations: a prime is excluded when an inequation polynomial lies
    # in it (then the inequation vanishes identically on that component -- the
    # cell forbids it).  param_ineqs are in QQ[constants][x,y,z]; jet_ineqs
    # (e.g. Ps != 0) are "nontrivial solution" conditions, not parameter ones.
    ineq_polys = [PolyRing(f) for f in cp['param_ineqs']]
    survivors = []
    for P in sc['primes']:
        if P.is_one():                       # unit ideal => empty variety
            continue
        if any(g in P for g in ineq_polys):  # inequation vanishes identically here
            continue
        survivors.append(P)

    tag = "TRIVIAL (Psi==0 forced)" if cell_trivial else "nontrivial"
    print(f"\n--- cell {num}: zero {{{', '.join(Zkey) or '(none, generic)'}}}; "
          f"ansatz {sc['spec_len']} eqs; "
          f"{len(cp['param_ineqs'])} param-ineqs, {len(cp['jet_ineqs'])} jet-ineqs; {tag} ---")
    if VERBOSE_REM:
        print("  remainder:", sc['rem'])
    if sc['rem'] == 0 and not cell_trivial:
        print("  PDE reduces to 0: the whole stratum solves the PDE (nontrivially)")
    print("  system of equations:", *sc['eqns'], sep="\n    " if sc['eqns'] else " (none)")
    if not survivors:
        print("  surviving solution varieties: NONE (all pruned by inequations / empty)")
    bucket = trivial_primes if cell_trivial else union_primes
    for P in survivors:
        print("   V:", P, "" if not cell_trivial else "   [trivial Psi==0]")
        key = prime_key(P)
        bucket.setdefault(key, (P, []))[1].append(num)

# --- union over all cells --------------------------------------------------
def dump_union(title, d):
    # Refine to the largest varieties: drop ideals that strictly contain another
    # in the list (their variety is already covered -> redundant in the union).
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
dump_union(f"NONTRIVIAL solution varieties over {scope} cells", union_primes)
print("\n" + "-" * 72)
dump_union(f"TRIVIAL (Psi==0) strata over {scope} cells", trivial_primes)
