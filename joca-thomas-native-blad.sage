# -*- mode: python -*-
#
# Sage script for the paper in the Journal of Computational Algebra.
#
# BLAD-NATIVE sibling of joca-thomas-native.sage.  Same pipeline (specialize the
# ansatz to each cell's parameter stratum, Ritt-reduce the PDE, project, GTZ
# minimal primes, prune by inequations, union) -- but the differential-Thomas
# layer is the standalone, native-Sage BLAD type
# (sage_differential_polynomial.DifferentialPolynomialRing +
# regularchains.diffthomas_blad) instead of
# DifferentialAlgebra.DifferentialRing + regularchains.diffthomas.
#
# The DifferentialPolynomial elements (BLAD `bap` objects) are the ring/ansatz/
# PDE carriers; the downstream prime pipeline stays sympy/Sage (build_system_of_
# equations / specialize / minimal_associated_primes operate on sympy IndexedBase
# jet expressions and the Sage PolyRing).  The two are joined by `_elt_to_sympy`,
# a term-walk (`_blad.read_terms`) that lowers a DifferentialPolynomial to the
# exact sympy jet expression the old DifferentialAlgebra layer produced
# (IndexedBase `Psi[x,x]`, params/coords as Symbols).  See the BRIDGE notes
# below.
#
#   *** RESERVED -- NOT RUN IN THE PORTING TASK ***
#
# This file is the N=29 cell-count follow-up.  The hydrogen decomposition is a
# ~18-minute run (29 components in the open-maple emulator) and is OUT OF SCOPE
# for the port that created this file; the port only verified that the ring,
# ansatz, and PDE construct + import cleanly on the BLAD-native path.  Running
# the actual decomposition is a separate, explicit top-level interactive run.
#
# Two representation differences from the DifferentialAlgebra layer, both forced
# by the package and both downstream-harmless (they preserve the cells' varieties
# and the rem==0 test):
#
#  * RANKING.  The old ring used blocks=[[DDPsi,DPsi,Psi,v,r], constants], which
#    is a SINGLE jet block (grlexA-within-block: order dominates, so DPsi[x]
#    outranks DDPsi).  The package mirror is therefore
#    ranking={'blocks': [['DDPsi','DPsi','Psi','v','r']]} -- ONE block, NOT a
#    block-per-indeterminate split.  (Block-per-indeterminate would make DDPsi
#    outrank every derivative of DPsi, giving leader(DPsi[x]-DDPsi*v[x]) = DDPsi
#    instead of the old DPsi[x].  Leader-agreement was checked against the old
#    ring -- see the port report.)
#
#  * COEFFICIENTS.  The package's BLAD `bap` body is over ZZ, not QQ: the string
#    parser truncates `/` (1/2 -> integer) and rational scalars collapse to int.
#    So the PDE's -1/2 factor is CLEARED by writing 2*PDE.  Scaling the PDE by a
#    nonzero rational does not change rem==0 nor the variety of the coefficient
#    relations build_system_of_equations extracts, so the cell decomposition and
#    the minimal primes are unaffected.
#
# Author: Brent Baccala
# Date: June 2026
# Tested on Ubuntu 24 with Sage 10.7 (conda env `sage`) and the BLAD-native
# sage_differential_polynomial / regularchains-sage packages.
#
#   sage joca-thomas-native-blad.sage [--max-rounds N] [--cells-out PATH]
#                                     [--verbose-remainder]

import os, sys, time

import sympy

# regularchains-sage: the BLAD-native differential Thomas engine, and the
# standalone native-Sage differential-polynomial ring it operates on.
sys.path.insert(0, os.path.expanduser('~/regularchains-sage'))
sys.path.insert(0, os.path.expanduser('~/sage-differential-polynomial/src'))
from sage_differential_polynomial import DifferentialPolynomialRing, _blad
from regularchains.diffthomas_blad import (differential_thomas_decomposition,
                                           appearing_derivatives)

# --- command-line options -------------------------------------------------
def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default

MAX_ROUNDS = int(_argval('--max-rounds', '5000'))
CELLS_OUT = _argval('--cells-out',
                    os.path.expanduser('~/thomas-experiments/hydrogen_thomas_native_blad.cells'))
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

# --- ring / PDE / ansatz --------------------------------------------------
# Sympy symbols for the downstream (prime) pipeline -- identical to joca.sage.
x, y, z = sympy.var('x,y,z')
E = sympy.var('E')
v1, v2, v3, v4 = sympy.var('v1,v2,v3,v4')
a0, a1, b0, b1, c0, c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]
# strata constants = everything but E (E lives only in the PDE, not the ansatz)
strata_constants = [v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]

# Sympy IndexedBase jets, mirroring DifferentialAlgebra.indexedbase, so the
# bridged-down remainder has exactly the form the old pipeline consumed.
Psi, DPsi, DDPsi = (sympy.IndexedBase('Psi'), sympy.IndexedBase('DPsi'),
                    sympy.IndexedBase('DDPsi'))
v = sympy.IndexedBase('v')
r = sympy.IndexedBase('r')

# The BLAD-native differential ring.  SINGLE jet block (see header): mirrors the
# old blocks=[[DDPsi,DPsi,Psi,v,r], constants] elimination/order ranking.
PARAM_NAMES = ['E', 'v1', 'v2', 'v3', 'v4', 'a0', 'a1', 'b0', 'b1', 'c0', 'c1']
R = DifferentialPolynomialRing(
    QQ,
    ['DDPsi', 'DPsi', 'Psi', 'v', 'r'],
    ['x', 'y', 'z'],
    parameters=PARAM_NAMES,
    ranking={'blocks': [['DDPsi', 'DPsi', 'Psi', 'v', 'r']]})

# PDE, denominators cleared (bap is over ZZ): this is 2 * the joca.sage PDE
#   joca:  -1/2*(Psi[x,x]+Psi[y,y]+Psi[z,z])*r - Psi - E*r*Psi
#   here:  -(Psi[x,x]+Psi[y,y]+Psi[z,z])*r - 2*Psi - 2*E*r*Psi
PDE = R('-(Psi[x,x] + Psi[y,y] + Psi[z,z])*r - 2*Psi - 2*E*r*Psi')
print("PDE (denominators cleared, x2):", PDE)

# Ansatz as BLAD-native package elements (jet string notation, e.g. Psi[x]).
ansatz0 = [R('Psi[x] - DPsi*v[x]'),
           R('Psi[y] - DPsi*v[y]'),
           R('Psi[z] - DPsi*v[z]'),
           R('DPsi[x] - DDPsi*v[x]'),
           R('DPsi[y] - DDPsi*v[y]'),
           R('DPsi[z] - DDPsi*v[z]'),
           R('(a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi'),
           R('v - (v1*x + v2*y + v3*z + v4*r)'),
           R('r^2 - x^2 - y^2 - z^2')]

# --- the BLAD <-> sympy bridge --------------------------------------------
# Maps for lowering a DifferentialPolynomial back to a sympy IndexedBase jet
# expression (the form build_system_of_equations / specialize / PolyRing want).
_IB = {'Psi': Psi, 'DPsi': DPsi, 'DDPsi': DDPsi, 'v': v, 'r': r}
_DERIV = {'x': x, 'y': y, 'z': z}
_PARAM = {p: sympy.Symbol(p) for p in PARAM_NAMES}


def _elt_to_sympy(e):
    """Lower a package DifferentialPolynomial to the sympy jet expression the
    downstream prime pipeline consumes.

    Walk the BLAD term iterator (`_blad.read_terms`, yielding
    ``(int_coeff, [(blad_name, deg), ...])``) and rebuild each monomial:

      * a jet ``'Psi[x,x]'`` -> sympy ``Psi[x, x]`` (IndexedBase indexed by the
        derivation symbols), mirroring DifferentialAlgebra's jet notation;
      * an order-zero head ``'Psi'`` / ``'v'`` / ``'r'`` -> the bare IndexedBase
        object (matches the old order-0 jet, which prints as the base);
      * a coordinate ``'x'`` -> the coord Symbol;
      * a parameter ``'E'`` / ``'v1'`` ... -> the constant Symbol.

    Coefficients are integers (the bap body is over ZZ); the PDE was scaled to
    clear denominators, so no rational ever needs to round-trip.
    """
    out = sympy.Integer(0)
    for coeff, term in _blad.read_terms(e._h()):
        mon = sympy.Integer(int(coeff))
        for nm, deg in term:
            if '[' in nm:
                head, rest = nm.split('[', 1)
                idx = tuple(_DERIV[d] for d in rest.rstrip(']').split(','))
                sym = _IB[head][idx]
            elif nm in _IB:
                sym = _IB[nm]            # order-0 head: bare IndexedBase
            elif nm in _DERIV:
                sym = _DERIV[nm]
            elif nm in _PARAM:
                sym = _PARAM[nm]
            else:
                raise KeyError("no sympy image for BLAD name %r" % (nm,))
            mon *= sym ** int(deg)
        out += mon
    return sympy.expand(out)


def _sympy_to_blad_str(expr):
    """Render a sympy jet expression as a BLAD ring-parser string.

    Indexed jets ``Psi[x, x]`` -> ``Psi[x,x]``; everything else prints as sympy
    does (symbols/coords/params keep their names, ``**`` -> ``^``).  Only used
    on the specialized ansatz, whose coefficients are integers after the
    denominator-cleared scaling, so no ``/`` ever appears."""
    import re
    s = sympy.sstr(expr)
    # sympy prints Indexed as 'Psi[x, x]'; BLAD wants 'Psi[x,x]'
    s = re.sub(r'\[([^\]]*)\]',
               lambda m: '[' + m.group(1).replace(' ', '') + ']', s)
    return s.replace('**', '^')


# --- Sage polynomial ring (identical to joca.sage) ------------------------
# Base generators only (derivations + heads + constants); the jets that survive
# in a fully-reduced remainder are all over these base generators.
PolyRing = PolynomialRing(QQ, names=(['x', 'y', 'z']
                                     + ['DDPsi', 'DPsi', 'Psi', 'v', 'r']
                                     + PARAM_NAMES))
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


# --- the differential-Thomas cells, computed natively (BLAD path) ---------
print(f"\nComputing BLAD-native differential Thomas decomposition of the "
      f"hydrogen ansatz ({len(ansatz0)} equations, param_field=False, "
      f"max_rounds={MAX_ROUNDS}) ...", flush=True)
_t0 = time.time()
cells_ds = differential_thomas_decomposition(R, ansatz0, char=0,
                                             param_field=False,
                                             max_rounds=MAX_ROUNDS,
                                             progress=True)
print(f"\nBLAD-native decomposition: {time.time() - _t0:.1f}s -> "
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
    return len(appearing_derivatives(R, [p])) > 0


def adapt_cell(ds):
    """A native DifferentialSystem -> the {param_eqs, param_ineqs, jet_ineqs}
    dict that the joca-thomas pipeline consumes.  Split each polynomial by
    whether it involves a jet (a differential/solution condition) or only the
    constants (a parameter-stratum / parameter-inequation).  Parameter
    polynomials are lowered to sympy via the bridge."""
    param_eqs, param_ineqs, jet_ineqs = [], [], []
    for e in ds.eqs:
        if has_jet(e):
            continue                       # solved-form jet eq -- pipeline re-derives it
        param_eqs.append(_elt_to_sympy(e))
    for q in ds.ineqs:
        if has_jet(q):
            jet_ineqs.append(q)            # e.g. Psi != 0 (nontrivial-solution)
        else:
            param_ineqs.append(_elt_to_sympy(q))
    return dict(param_eqs=param_eqs, param_ineqs=param_ineqs, jet_ineqs=jet_ineqs)


def specialize(param_eqs):
    """The cell-specialized ansatz: solve the parameter relations for as many
    constants as possible and substitute into ansatz0.

    Returns ``(sub, spec)`` where ``sub`` is the sympy substitution dict and
    ``spec`` is the specialized ansatz as a list of PACKAGE elements (the
    substitution is applied in sympy, then each substituted ansatz equation is
    re-parsed back into the ring), so it can feed R.differential_prem directly."""
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
    # ansatz0 are package elements; lower to sympy, substitute, re-lift to R.
    spec = []
    for a in ansatz0:
        a_sym = sympy.expand(_elt_to_sympy(a).subs(sub))
        if a_sym == 0:
            continue
        spec.append(R(_sympy_to_blad_str(a_sym)))
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
        # R.differential_prem returns (remainder, cofactor) -- NOTE the order is
        # flipped vs the old DifferentialAlgebra (which returned (cofactor, rem)).
        rem_elt, h = R.differential_prem(PDE, spec)
        # Does the specialized ansatz force the wavefunction to vanish?
        psi_rem_elt, _ = R.differential_prem(R('Psi'), spec)
        trivial = psi_rem_elt.is_zero()
        rem = _elt_to_sympy(rem_elt)        # bridge down for the prime pipeline
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
