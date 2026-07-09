# -*- mode: python -*-
#
# Sage script for the paper in the Journal of Computational Algebra.
#
# NATIVE-DIFFERENTIAL-THOMAS sibling of joca-thomas-native-blad.sage.  Same
# downstream pipeline (specialize the ansatz to each cell's parameter stratum,
# Ritt-reduce the PDE, project, GTZ minimal primes, prune by inequations, union)
# -- but the differential-Thomas layer is the standalone native-Sage port of
# Lange-Hegermann's DifferentialThomas (``~/DifferentialThomas-sage``:
# ``differentialthomas.differential_thomas_decomposition``), NOT the
# regularchains CTD-based ``diffthomas_blad`` engine.
#
# WHY THIS DRIVER EXISTS.  ``joca-thomas-native-blad.sage`` uses
# ``regularchains.diffthomas_blad`` which, at ``param_field=False`` (parameters
# ranked -- the true 29-cell parametric decomposition the paper needs), suffers a
# branch-count explosion (42 000+ branches, 24 GB, never completes -- see
# ``~/project/reports/`` regularchains hydrogen memory profiles).  The
# DifferentialThomas port is a faithful completion engine (selection strategy +
# eager inconsistency pruning + implied-inequation suppression) and computes the
# SAME 29 cells without exploding: peak work-queue ~14, ~4 min.  Its 29 cells
# match the open-maple reference decomposition cell-for-cell (verified:
# ``~/project/reports/differentialthomas-sage-phase6-hydrogen.md``).
#
# The differential-Thomas cells are native-Sage ``DifferentialPolynomial``
# elements over ONE ``sage_differential_polynomial.DifferentialPolynomialRing``;
# the downstream prime pipeline stays sympy/Sage.  They are joined by
# ``_elt_to_sympy`` -- the same term-walk bridge (``_blad.read_terms``) as the
# blad driver -- which lowers a DifferentialPolynomial to the sympy IndexedBase
# jet expression the old DifferentialAlgebra layer produced.
#
# Author: Brent Baccala (AI assistant: Claude)
# Date: July 2026
# Tested on Ubuntu 24 with Sage 10.7 (conda env `sage`) and the native-Sage
# DifferentialThomas-sage / sage_differential_polynomial packages.
#
#   sage joca-thomas-native-dt.sage [--cells-out PATH] [--verbose-remainder]
#                                   [--decompose-only]

import os, sys, time

import sympy

# The native-Sage DifferentialThomas port and the differential-polynomial ring
# it operates on.
sys.path.insert(0, os.path.expanduser('~/DifferentialThomas-sage'))
sys.path.insert(0, os.path.expanduser('~/sage-differential-polynomial/src'))
import differentialthomas as dt
from sage_differential_polynomial import _blad


# --- command-line options -------------------------------------------------
def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default

CELLS_OUT = _argval('--cells-out',
                    os.path.expanduser('~/thomas-experiments/hydrogen_thomas_native_dt.cells'))
VERBOSE_REM = '--verbose-remainder' in sys.argv
DECOMPOSE_ONLY = '--decompose-only' in sys.argv
MAX_CELLS = int(_argval('--max-cells', '0'))     # 0 = all cells


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

# Sympy IndexedBase jets, mirroring DifferentialAlgebra.indexedbase.
Psi, DPsi, DDPsi = (sympy.IndexedBase('Psi'), sympy.IndexedBase('DPsi'),
                    sympy.IndexedBase('DDPsi'))
v = sympy.IndexedBase('v')
r = sympy.IndexedBase('r')

# The differential ranking + ring.  FLAT DegRevLex over [jets..., E, params...]
# -- exactly the ranking the open-maple reference ex4_hydrogen.mpl uses
# (``Ranking(ivars, [op(jets), op(pars)])``: jets high, params low, all ranked =
# param_field=False), which is the ranking under which the port reproduces the
# verified 29 cells.  E is ranked at the bottom of the jets and gets NO constancy
# equation; it never appears in the ansatz decomposition (it enters only the
# downstream PDE) and stays inert, so the 29-cell result is unaffected.
JETS = ['DDPsi', 'DPsi', 'Psi', 'v', 'r']
ANSATZ_PARAMS = ['v1', 'v2', 'v3', 'v4', 'a0', 'a1', 'b0', 'b1', 'c0', 'c1']
DVAR = JETS + ['E'] + ANSATZ_PARAMS
IVAR = ['x', 'y', 'z']
rk = dt.compute_ranking(IVAR, DVAR)
R = rk.ring

# PDE, denominators cleared (bap is over ZZ): this is 2 * the joca.sage PDE.
PDE = R('-(Psi[x,x] + Psi[y,y] + Psi[z,z])*r - 2*Psi - 2*E*r*Psi')
print("PDE (denominators cleared, x2):", PDE)

# Ansatz as native ring elements (jet string notation, e.g. Psi[x]).
ansatz0 = [R('Psi[x] - DPsi*v[x]'),
           R('Psi[y] - DPsi*v[y]'),
           R('Psi[z] - DPsi*v[z]'),
           R('DPsi[x] - DDPsi*v[x]'),
           R('DPsi[y] - DDPsi*v[y]'),
           R('DPsi[z] - DDPsi*v[z]'),
           R('(a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi'),
           R('v - (v1*x + v2*y + v3*z + v4*r)'),
           R('r^2 - x^2 - y^2 - z^2')]

# Parameter-constancy equations: every first derivative of each ANSATZ parameter
# is zero (the 10 params are constants).  These mirror the reference's ``pconst``
# (10 params x 3 ivars = 30 equations).  E gets none (inert).
# BLAD derivation notation for a first derivative is ``p[x]`` (not the Maple
# exponent form ``p[1,0,0]``): the substrate ring parser consumes derivation
# names.
pconst = [R('%s[%s]' % (p, iv)) for p in ANSATZ_PARAMS for iv in IVAR]

# --- the native <-> sympy bridge ------------------------------------------
_IB = {'Psi': Psi, 'DPsi': DPsi, 'DDPsi': DDPsi, 'v': v, 'r': r}
_DERIV = {'x': x, 'y': y, 'z': z}
_PARAM = {p: sympy.Symbol(p) for p in (['E'] + ANSATZ_PARAMS)}
_JET_HEADS = set(_IB)                     # heads that make a poly "have a jet"


def _elt_to_sympy(e):
    """Lower a native DifferentialPolynomial to the sympy jet expression the
    downstream prime pipeline consumes (identical term-walk to the blad
    driver)."""
    out = sympy.Integer(0)
    for coeff, term in _blad.read_terms(e._h()):
        mon = sympy.Integer(int(coeff))
        for nm, deg in term:
            if '[' in nm:
                head, rest = nm.split('[', 1)
                idx = tuple(_DERIV[d] for d in rest.rstrip(']').split(','))
                sym = _IB[head][idx]
            elif nm in _IB:
                sym = _IB[nm]
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
    """Render a sympy jet expression as a BLAD ring-parser string."""
    import re
    s = sympy.sstr(expr)
    s = re.sub(r'\[([^\]]*)\]',
               lambda m: '[' + m.group(1).replace(' ', '') + ']', s)
    return s.replace('**', '^')


def has_jet(p):
    """True iff the polynomial involves a differential indeterminate (a jet
    Psi/DPsi/DDPsi/v/r or a derivative of one), as opposed to a pure parameter
    relation in the constants.  Replaces the blad driver's
    ``appearing_derivatives``: walk the term iterator and look for any name whose
    head is a jet head."""
    for _coeff, term in _blad.read_terms(p._h()):
        for nm, _deg in term:
            head = nm.split('[', 1)[0]
            if head in _JET_HEADS:
                return True
    return False


def is_param_constancy(p):
    """True iff the polynomial is a parameter-constancy condition -- it involves
    a DERIVATIVE of a parameter (e.g. ``v1[x]``, i.e. dv1/dx = 0).  Because the
    port ranks the 10 parameters as full differential variables (param_field=
    False), each cell carries the explicit constancy equations ``p[x]=p[y]=
    p[z]=0``; the reference's ``parameters=`` declaration made these implicit.
    They carry NO parameter-stratum information (the downstream pipeline already
    treats each parameter as a single constant symbol) and have no sympy image,
    so ``adapt_cell`` drops them."""
    for _coeff, term in _blad.read_terms(p._h()):
        for nm, _deg in term:
            if '[' in nm and nm.split('[', 1)[0] in _PARAM:
                return True
    return False


# --- Sage polynomial ring (identical to joca.sage) ------------------------
PolyRing = PolynomialRing(QQ, names=(['x', 'y', 'z']
                                     + ['DDPsi', 'DPsi', 'Psi', 'v', 'r']
                                     + ['E'] + ANSATZ_PARAMS))
PolyRing_constants = list(map(PolyRing, constants))


def build_system_of_equations(eqn, constants):
    """Factor each term into constant * non-constant parts, gather like
    non-constant monomials; return the (constant-only) coefficient polynomials.
    Identical to joca.sage."""
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


# --- the differential-Thomas cells, computed by the native port -----------
print(f"\nComputing native DifferentialThomas decomposition of the hydrogen "
      f"ansatz ({len(ansatz0)} ansatz + {len(pconst)} constancy = "
      f"{len(ansatz0) + len(pconst)} equations, param_field=False) ...",
      flush=True)
print("  (set DT_BRANCH_TRACE=1 for the branch/RSS trajectory on stderr)",
      flush=True)
_t0 = time.time()
cells_ds = dt.differential_thomas_decomposition(ansatz0 + pconst, [], rk)
_wall = time.time() - _t0
print(f"\nNative DifferentialThomas decomposition: {_wall:.1f}s -> "
      f"{len(cells_ds)} cells\n" + "=" * 72, flush=True)

if len(cells_ds) != 29:
    print(f"WARNING: expected 29 cells (the verified reference count), got "
          f"{len(cells_ds)}.", flush=True)


def cell_eqs(ds):
    return list(dt.differential_system_equations(ds))


def cell_ineqs(ds):
    return list(dt.differential_system_inequations(ds))


# --- per-cell soundness cross-check (triangularity) -----------------------
# Engine-independent Thomas-system soundness: within a cell the equations have
# pairwise-distinct leaders (a triangular set).  Disjointness of the whole
# decomposition is inherited from the cell-for-cell match to the open-maple
# reference (verified independently disjoint/complete); see the phase-6 report.
_sound = True
for i, ds in enumerate(cells_ds, 1):
    leaders = [e.leader() for e in cell_eqs(ds) if not e.is_zero()]
    jet_leaders = [L for L in leaders if L is not None]
    if len(jet_leaders) != len(set(jet_leaders)):
        _sound = False
        print(f"  SOUNDNESS FAIL cell {i}: repeated leader among {jet_leaders}",
              flush=True)
print(f"per-cell triangularity (distinct leaders): "
      f"{'OK' if _sound else 'FAIL'} over {len(cells_ds)} cells", flush=True)

# checkpoint the raw cells so the decomposition is never lost
try:
    with open(CELLS_OUT, 'w') as fh:
        for i, ds in enumerate(cells_ds, 1):
            eqs, ineqs = cell_eqs(ds), cell_ineqs(ds)
            fh.write("--- cell %d : %d eqs, %d ineqs ---\n"
                     % (i, len(eqs), len(ineqs)))
            fh.write("EQS: %s\nINEQS: %s\n\n"
                     % ([str(e) for e in eqs], [str(q) for q in ineqs]))
    print(f"Wrote raw native cells to {CELLS_OUT}", flush=True)
except Exception as ex:
    print(f"(could not write cells file: {ex})", flush=True)

if DECOMPOSE_ONLY:
    print("\n--decompose-only: stopping before the prime pipeline.", flush=True)
    sys.exit(0)


def adapt_cell(ds):
    """A native DifferentialSystem -> the {param_eqs, param_ineqs, jet_ineqs}
    dict the joca-thomas pipeline consumes."""
    param_eqs, param_ineqs, jet_ineqs = [], [], []
    for e in cell_eqs(ds):
        if e.is_zero():
            continue
        if has_jet(e):
            continue                       # solved-form jet eq -- re-derived
        if is_param_constancy(e):
            continue                       # p[x]=0 constancy -- no stratum info
        param_eqs.append(_elt_to_sympy(e))
    for q in cell_ineqs(ds):
        if has_jet(q):
            jet_ineqs.append(q)            # e.g. Psi != 0 (nontrivial-solution)
        elif is_param_constancy(q):
            continue
        else:
            param_ineqs.append(_elt_to_sympy(q))
    return dict(param_eqs=param_eqs, param_ineqs=param_ineqs,
                jet_ineqs=jet_ineqs)


def specialize(param_eqs):
    """The cell-specialized ansatz: solve the parameter relations and substitute
    into ansatz0.  Returns ``(sub, spec)`` (spec = specialized ansatz as ring
    elements)."""
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
    spec = []
    for a in ansatz0:
        a_sym = sympy.expand(_elt_to_sympy(a).subs(sub))
        if a_sym == 0:
            continue
        spec.append(R(_sympy_to_blad_str(a_sym)))
    return sub, spec


def reductors_for(spec):
    """The differential-prem reductor set: the specialized ansatz PLUS the
    parameter-constancy equations ``p[x]=p[y]=p[z]=0``.  The constancy reductors
    are essential here (unlike the blad driver, where ``parameters=`` made the
    params order-zero): without them the reduction of the PDE by the ansatz
    introduces parameter derivatives (``v1[x]`` from differentiating
    ``v - (v1*x+...)``), which have no downstream image.  With them, every
    parameter derivative reduces to 0, exactly reproducing the order-zero-
    parameter behaviour of the reference."""
    return list(spec) + list(pconst)


def prime_key(P):
    return tuple(sorted(str(g) for g in P.gens()))


# Cache reduction+decomposition by parameter stratum.
strata_cache = {}
union_primes = {}
trivial_primes = {}

_downstream_cells = cells_ds if MAX_CELLS <= 0 else cells_ds[:MAX_CELLS]
if MAX_CELLS > 0:
    print(f"\n(--max-cells {MAX_CELLS}: running the prime pipeline on the first "
          f"{len(_downstream_cells)} of {len(cells_ds)} cells to demonstrate "
          f"the hand-off)", flush=True)

for num, ds in enumerate(_downstream_cells, 1):
    cp = adapt_cell(ds)
    Z = sorted((p for p in cp['param_eqs']), key=str)
    Zkey = tuple(map(str, Z))

    if Zkey not in strata_cache:
        sub, spec = specialize(cp['param_eqs'])
        reductors = reductors_for(spec)
        rem_elt, h = R.differential_prem(PDE, reductors)
        psi_rem_elt, _ = R.differential_prem(R('Psi'), reductors)
        trivial = psi_rem_elt.is_zero()
        rem = _elt_to_sympy(rem_elt)
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
          f"{len(cp['param_ineqs'])} param-ineqs, {len(cp['jet_ineqs'])} "
          f"jet-ineqs; {tag} ---")
    if VERBOSE_REM:
        print("  remainder:", sc['rem'])
    if sc['rem'] == 0 and not sc['trivial']:
        print("  PDE reduces to 0: the whole stratum solves the PDE (nontrivially)")
    print("  system of equations:", *sc['eqns'],
          sep="\n    " if sc['eqns'] else " (none)")
    if not survivors:
        print("  surviving solution varieties: NONE (all pruned / empty)")
    bucket = trivial_primes if sc['trivial'] else union_primes
    for P in survivors:
        print("   V:", P, "" if not sc['trivial'] else "   [trivial Psi==0]")
        key = prime_key(P)
        bucket.setdefault(key, (P, []))[1].append(num)


def dump_union(title, d):
    print(f"\n{title} ({len(d)} distinct primes):\n")
    for key, (P, cells_for) in sorted(d.items(), key=lambda kv: str(kv[0])):
        print("  V:", P)
        print("       (from cells:",
              ", ".join(map(str, sorted(set(cells_for)))) + ")")

print("\n" + "=" * 72)
dump_union("NONTRIVIAL solution varieties over all cells", union_primes)
print("\n" + "-" * 72)
dump_union("TRIVIAL (Psi==0) strata over all cells", trivial_primes)
