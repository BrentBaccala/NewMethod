# -*- mode: python -*-
#
# Sage script for the paper in the Journal of Computational Algebra.
#
# Native-Sage variant of joca-thomas.sage.  Instead of loading the differential
# Thomas cells of the hydrogen ansatz from Maple's DifferentialThomas log
# (~/thomas-experiments/hydrogen_thomas.log), it computes them *in Sage* with
# regularchains-sage's differential_thomas_decomposition (param_field=False:
# the constants are ranked variables, so the decomposition stratifies on them).
# Everything downstream -- specialize the ansatz to each cell's parameter
# stratum, Ritt-reduce the PDE, project, GTZ minimal primes, prune by
# inequations, union -- is byte-for-byte the joca-thomas.sage pipeline.
#
# This removes the last Maple dependency from the completeness oracle: the
# cells are produced by the same engine (BLAD + Layer-0 algebraic Thomas) the
# rest of the pipeline already uses.
#
# Author: Brent Baccala
# Date: June 2026
# Tested on Ubuntu 24 with Sage 10.7 and DifferentialAlgebra 1.14.
#
#   sage joca-thomas-native.sage [--max-rounds N] [--cells-out PATH]
#                                [--verbose-remainder]

import os, sys, time

import sympy

try:
    import DifferentialAlgebra
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(ex.msg + "\nInstall it with '%pip install DifferentialAlgebra'")

# regularchains-sage: the native differential Thomas engine
sys.path.insert(0, os.path.expanduser('~/regularchains-sage'))
from regularchains.diffthomas import (differential_thomas_decomposition,
                                       appearing_derivatives)

# --- command-line options -------------------------------------------------
def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default

MAX_ROUNDS = int(_argval('--max-rounds', '5000'))
CELLS_OUT = _argval('--cells-out',
                    os.path.expanduser('~/thomas-experiments/hydrogen_thomas_native.cells'))
VERBOSE_REM = '--verbose-remainder' in sys.argv

# Claude Sonnet 4's solution to printing "DPsi" as "\Psi'" in LaTeX
def patch_latex_varify():
    from sage.misc.latex import latex_varify
    import sage.misc.latex
    original = latex_varify
    def custom(a, is_fname=False):
        return r"\Psi'" if a == "DPsi" else original(a, is_fname=is_fname)
    sage.misc.latex.latex_varify = custom
patch_latex_varify()

# --- ring / PDE / ansatz (identical to joca.sage / joca-thomas.sage) -------
x, y, z = sympy.var('x,y,z')
E = sympy.var('E')
v1, v2, v3, v4 = sympy.var('v1,v2,v3,v4')
a0, a1, b0, b1, c0, c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]
# strata constants = everything but E (E lives only in the PDE, not the ansatz)
strata_constants = [v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]

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
    Identical to joca.sage / joca-rg.sage / joca-thomas.sage."""
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


# --- the differential-Thomas cells, computed natively ---------------------
print(f"\nComputing native differential Thomas decomposition of the hydrogen "
      f"ansatz ({len(ansatz0)} equations, param_field=False, "
      f"max_rounds={MAX_ROUNDS}) ...", flush=True)
_t0 = time.time()
cells_ds = differential_thomas_decomposition(DiffRing, ansatz0, char=0,
                                             param_field=False,
                                             max_rounds=MAX_ROUNDS,
                                             progress=True)
print(f"\nNative decomposition: {time.time() - _t0:.1f}s -> "
      f"{len(cells_ds)} cells\n" + "=" * 72, flush=True)

# checkpoint the raw cells so the (expensive) decomposition is never lost
try:
    with open(CELLS_OUT, 'w') as fh:
        for i, ds in enumerate(cells_ds, 1):
            fh.write("--- cell %d : %d eqs, %d ineqs ---\n"
                     % (i, len(ds.eqs), len(ds.ineqs)))
            fh.write("EQS: %s\nINEQS: %s\n\n"
                     % ([str(e) for e in ds.eqs], [str(q) for q in ds.ineqs]))
    print(f"Wrote raw native cells to {CELLS_OUT}", flush=True)
except Exception as ex:
    print(f"(could not write cells file: {ex})", flush=True)


def has_jet(p):
    """True iff the polynomial involves a differential indeterminate (jet),
    i.e. a Psi/DPsi/DDPsi/v/r derivative -- as opposed to a pure parameter
    relation in the constants."""
    return len(appearing_derivatives(DiffRing, [p])) > 0


def adapt_cell(ds):
    """A native DifferentialSystem -> the {param_eqs, param_ineqs, jet_ineqs}
    dict that the joca-thomas pipeline consumes.  The inverse of
    thomas_cells.cell_polys's classification: split by whether a polynomial
    involves a jet (then it's a differential/solution condition) or only the
    constants (then it's a parameter-stratum / parameter-inequation)."""
    param_eqs, param_ineqs, jet_ineqs = [], [], []
    for e in ds.eqs:
        if has_jet(e):
            continue                       # solved-form jet eq -- pipeline re-derives it
        param_eqs.append(sympy.expand(e))
    for q in ds.ineqs:
        if has_jet(q):
            jet_ineqs.append(q)            # e.g. Psi != 0 (nontrivial-solution)
        else:
            param_ineqs.append(sympy.expand(q))
    return dict(param_eqs=param_eqs, param_ineqs=param_ineqs, jet_ineqs=jet_ineqs)


def specialize(param_eqs):
    """The cell-specialized ansatz: solve the parameter relations for as many
    constants as possible and substitute into ansatz0.  Generalizes
    joca-thomas.sage's `sub = {p:0 for p in Z}` (which assumed each stratum is a
    subset of constants set to 0) to arbitrary parameter relations."""
    sub = {}
    if param_eqs:
        try:
            sols = sympy.solve(list(param_eqs), list(strata_constants), dict=True)
            if sols:
                sub = sols[0]
        except Exception:
            sub = {}
        if not sub:
            sub = {p: 0 for p in param_eqs if p in set(strata_constants)}
    spec = [e for e in (sympy.expand(a.subs(sub)) for a in ansatz0) if e != 0]
    return sub, spec


def prime_key(P):
    """A canonical, sortable string for a prime ideal (for dedup / union)."""
    return tuple(sorted(str(g) for g in P.gens()))


# Cache reduction+decomposition by parameter stratum (cells with the same
# parameter relations share an identical specialized ansatz, reduction, primes;
# only their inequations -- and thus the pruning -- differ).
strata_cache = {}
union_primes = {}            # prime_key -> (ideal, contributing cells)  -- nontrivial
trivial_primes = {}          # prime_key -> (ideal, contributing cells)  -- Psi==0 strata

for num, ds in enumerate(cells_ds, 1):
    cp = adapt_cell(ds)
    Z = sorted((p for p in cp['param_eqs']), key=str)
    Zkey = tuple(map(str, Z))

    if Zkey not in strata_cache:
        sub, spec = specialize(cp['param_eqs'])
        h, rem = DiffRing.differential_prem(PDE, spec)
        # Does the specialized ansatz force the wavefunction to vanish?
        _, psi_rem = DiffRing.differential_prem(Psi, spec)
        trivial = (psi_rem == 0)
        if rem == 0:
            eqns = ()
        else:
            eqns = build_system_of_equations(PolyRing(rem), PolyRing_constants)
        gens = list(eqns) + [PolyRing(p) for p in Z]
        I = ideal(gens) if gens else ideal(PolyRing.zero())
        primes = I.minimal_associated_primes()
        strata_cache[Zkey] = dict(spec_len=len(spec), rem=rem, eqns=eqns,
                                  primes=primes, trivial=trivial)

    sc = strata_cache[Zkey]

    ineq_polys = [PolyRing(f) for f in cp['param_ineqs']]
    survivors = []
    for P in sc['primes']:
        if P.is_one():
            continue
        if any(g in P for g in ineq_polys):
            continue
        survivors.append(P)

    tag = "TRIVIAL (Psi==0 forced)" if sc['trivial'] else "nontrivial"
    print(f"\n--- cell {num}: zero {{{', '.join(Zkey) or '(none, generic)'}}}; "
          f"ansatz {sc['spec_len']} eqs; "
          f"{len(cp['param_ineqs'])} param-ineqs, {len(cp['jet_ineqs'])} jet-ineqs; {tag} ---")
    if VERBOSE_REM:
        print("  remainder:", sc['rem'])
    if sc['rem'] == 0 and not sc['trivial']:
        print("  PDE reduces to 0: the whole stratum solves the PDE (nontrivially)")
    print("  system of equations:", *sc['eqns'], sep="\n    " if sc['eqns'] else " (none)")
    if not survivors:
        print("  surviving solution varieties: NONE (all pruned by inequations / empty)")
    bucket = trivial_primes if sc['trivial'] else union_primes
    for P in survivors:
        print("   V:", P, "" if not sc['trivial'] else "   [trivial Psi==0]")
        key = prime_key(P)
        bucket.setdefault(key, (P, []))[1].append(num)

# --- union over all cells --------------------------------------------------
def dump_union(title, d):
    print(f"\n{title} ({len(d)} distinct primes):\n")
    for key, (P, cells_for) in sorted(d.items(), key=lambda kv: str(kv[0])):
        print("  V:", P)
        print("       (from cells:", ", ".join(map(str, sorted(set(cells_for)))) + ")")

print("\n" + "=" * 72)
dump_union("NONTRIVIAL solution varieties over all cells", union_primes)
print("\n" + "-" * 72)
dump_union("TRIVIAL (Psi==0) strata over all cells", trivial_primes)
