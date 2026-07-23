# -*- mode: python -*-
#
# thomas-ansatz-solve.sage
# ------------------------------------------------------------------------
# Solve an (ansatz, PDE) pair by the NewMethod staged pipeline of
# prolongation-projection-algorithm.tex (Phase III-forall / membership, run
# per Thomas cell -- the "staged route" of section 7):
#
#   1. pick one ansatz + one PDE                    (ansatz-library.sage)
#   2. differential-Thomas-decompose the ANSATZ ALONE  -> disjoint cells
#   3. reduce the PDE modulo each cell              (differential pseudo-remainder)
#   4. forall-project each remainder onto the constants: collect like terms in
#      the independents + parametric jets, zero the constant coefficients, take
#      minimal associated primes, prune by the cell's inequations
#   5. union the surviving varieties over all cells
#   6. print the union
#
# The ansatz and PDE are supplied by ansatz-library.sage in the differential-
# algebra formulation (differential-polynomial equations), NOT by helium.sage.
# The downstream (steps 2-6) is joca-thomas-native-dt.sage's proven pipeline,
# generalised to the problem's jet/param names.
#
# Validation: `--pde hydrogen --ansatz 5` must reproduce the verified 29-cell
# hydrogen decomposition and its known union of solution varieties.
#
#   sage thomas-ansatz-solve.sage --pde hydrogen --ansatz 5 [--decompose-only]
#   sage thomas-ansatz-solve.sage --pde helium   --ansatz 5
#
# Author: Brent Baccala (AI assistant: Claude).  July 2026.

import os, re, sys, time
import sympy

sys.path.insert(0, os.path.expanduser('~/DifferentialThomas-sage'))
sys.path.insert(0, os.path.expanduser('~/sage-differential-polynomial/src'))
import differentialthomas as dt
from sage_differential_polynomial import _blad

_HERE = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv[0] else '.'
load(os.path.join(_HERE, 'ansatz-library.sage'))


# --- command-line options -------------------------------------------------
def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default

PDE_NAME = _argval('--pde', 'hydrogen')
ANSATZ = _argval('--ansatz', '5')
ANSATZ = float(ANSATZ) if '.' in str(ANSATZ) else int(ANSATZ)
DECOMPOSE_ONLY = '--decompose-only' in sys.argv
VERBOSE_REM = '--verbose-remainder' in sys.argv
MAX_CELLS = int(_argval('--max-cells', '0'))
# Differential ranking: 'orderly' (degrevlex, coordinate-order dominates) or
# 'elimination' (block ranking, each jet its own block so a high jet like DDPsi
# outranks all lower-jet derivatives).  The ranking changes the decomposition,
# so the cells file name records which one was used.
RANKING = _argval('--ranking', 'orderly')
# By default, prune GENUINE varieties whose zero-set is contained in another
# genuine variety's, printing only the maximal (enclosing) ones.  Two "genuine"
# primes coming out of different cells are often nested (e.g. the c1=0 wall of a
# larger E=-1/2 component reappears as its own prime), which double-counts one
# solution family.  --keep-enclosed restores the raw per-cell union.
PRUNE_ENCLOSED = '--keep-enclosed' not in sys.argv
CELLS_OUT = _argval('--cells-out',
                    os.path.expanduser('~/thomas-experiments/%s_ansatz%s_%s.cells'
                                       % (PDE_NAME, ANSATZ, RANKING)))
os.makedirs(os.path.dirname(CELLS_OUT), exist_ok=True)


def patch_latex_varify():
    from sage.misc.latex import latex_varify
    import sage.misc.latex
    original = latex_varify
    def custom(a, is_fname=False):
        return r"\Psi'" if a == "DPsi" else original(a, is_fname=is_fname)
    sage.misc.latex.latex_varify = custom
patch_latex_varify()


# ==========================================================================
# the problem (from the differential-algebra library)
# ==========================================================================
print("Building %s / ansatz %s from the differential-algebra library ..."
      % (PDE_NAME, ANSATZ), flush=True)
prob = build_problem(PDE_NAME, ANSATZ, ranking=RANKING)

R = prob['R']
COORDS = prob['coords']
ROOT_NAMES = [rn for rn, _ in prob['roots']]
PARAMS = prob['params']
JETS = prob['jets']                          # all differential indeterminates
ansatz0 = prob['ansatz_eqs']
pconst = prob['pconst']
PDE = prob['pde']

print("PDE:", PDE)
print("ansatz (%d eqs):" % len(ansatz0))
for s in prob['ansatz_eqs_str']:
    print("   ", s)
print("params:", PARAMS)


# ==========================================================================
# native <-> sympy bridge (dynamic names), then the DT pipeline
# ==========================================================================
constants = [sympy.Symbol('E')] + [sympy.Symbol(p) for p in PARAMS]
strata_constants = [sympy.Symbol(p) for p in PARAMS]

_IB = {nm: sympy.IndexedBase(nm) for nm in JETS}
_DERIV = {c: sympy.Symbol(c) for c in COORDS}
_PARAM = {p: sympy.Symbol(p) for p in (['E'] + PARAMS)}
_JET_HEADS = set(_IB)


def _elt_to_sympy(e):
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
    s = sympy.sstr(expr)
    s = re.sub(r'\[([^\]]*)\]',
               lambda m: '[' + m.group(1).replace(' ', '') + ']', s)
    return s.replace('**', '^')


def has_jet(p):
    for _coeff, term in _blad.read_terms(p._h()):
        for nm, _deg in term:
            if nm.split('[', 1)[0] in _JET_HEADS:
                return True
    return False


def is_param_constancy(p):
    for _coeff, term in _blad.read_terms(p._h()):
        for nm, _deg in term:
            if '[' in nm and nm.split('[', 1)[0] in _PARAM:
                return True
    return False


# Under the ELIMINATION ranking, full_prem's normal form may legitimately retain
# parametric DERIVATIVE jets: the block ranking makes the order-0 jets the
# leaders of the chain rules (DPsi outranks Psi[R1], so `Psi[R1] - DPsi*v[R1]`
# rewrites DPsi), leaving derivative jets like Psi[R12] under the staircase.
# Those surviving jets are the cell's parametric derivatives -- free on the cell
# exactly like the order-0 jets are under the orderly ranking -- so they get
# PolyRing generators of their own (bracket-free mangled names; the BLAD-name ->
# generator map is patched into _GEN_IDX below).  Order <= 3 covers everything
# the order-2 PDEs can leave; a higher-order survivor still raises the explicit
# TypeError in _elt_to_polyring.  Gated on the ranking so the orderly path is
# byte-identical to before.
EXTRA_JET_PAIRS = []                 # [(blad_name, mangled_generator_name)]
if RANKING in ('elimination', 'block', 'elim'):
    from itertools import combinations_with_replacement
    for _head in prob['jets']:
        for _k in (1, 2, 3):
            for _idx in combinations_with_replacement(COORDS, _k):
                EXTRA_JET_PAIRS.append(('%s[%s]' % (_head, ','.join(_idx)),
                                        '%s_%s' % (_head, '_'.join(_idx))))

PolyRing = PolynomialRing(QQ, names=(COORDS + prob['jets']
                                     + [m for _b, m in EXTRA_JET_PAIRS]
                                     + ['E'] + PARAMS))
PolyRing_constants = list(map(PolyRing, [str(c) for c in constants]))
V_PARAM_GENS = [PolyRing(p) for p in prob['v_params']]
AMP_PARAM_GENS = [PolyRing(p) for p in prob.get('amp_params', [])]

# sympy images of PolyRing's generators, in generator order.  _elt_to_sympy maps
# coords -> Symbol, bare jets -> IndexedBase, E/params -> Symbol, so this list is
# the exact mirror of PolyRing.gens().
def _blad_jet_to_sympy(blad_name):
    head, rest = blad_name.split('[', 1)
    return _IB[head][tuple(_DERIV[d] for d in rest.rstrip(']').split(','))]


_SYMPY_GENS = ([_DERIV[c] for c in COORDS]
               + [_IB[j] for j in prob['jets']]
               + [_blad_jet_to_sympy(b) for b, _m in EXTRA_JET_PAIRS]
               + [_PARAM['E']]
               + [_PARAM[p] for p in PARAMS])
assert len(_SYMPY_GENS) == PolyRing.ngens()


def _sympy_to_polyring(expr):
    """Convert a sympy expression to PolyRing WITHOUT a string round-trip.

    Sage has no sympy -> libsingular conversion, so `PolyRing(expr)` falls
    through to `self(str(expr))` and hands the result to `eval`.  CPython parses
    a sum as a left-nested tree of binary ops and its *compiler* recurses once
    per term, so that path dies with

        RecursionError: maximum recursion depth exceeded during compilation

    at ~2995 additive terms (python 3.11, recursionlimit 1000) -- regardless of
    how simple the polynomial is.  helium/ansatz 20.1 crossed that ceiling.
    Building from an exponent->coefficient dict never invokes the parser, so it
    has no such limit (measured: 4000 monomials in 0.04s).
    """
    expr = sympy.expand(expr)
    if expr == 0:
        return PolyRing.zero()

    # A derivative jet (Psi[R1,R1]) has no PolyRing generator -- the reduction is
    # meant to eliminate every one of them.  The old string path masked a leftover
    # behind the RecursionError above (compile fails before name resolution), so
    # check explicitly and name the offenders.
    # A bare jet's image is IndexedBase('Psi'), whose .label is Symbol('Psi'), and
    # that label shows up in .atoms(Symbol).  Admit the labels too, or every bare
    # jet reads as unknown; a derivative jet is an Indexed and is still caught.
    allowed = set(_SYMPY_GENS)
    allowed |= {g.label for g in _SYMPY_GENS if isinstance(g, sympy.IndexedBase)}
    unknown = {a for a in expr.atoms(sympy.Symbol, sympy.Indexed, sympy.IndexedBase)
               if a not in allowed}
    if unknown:
        raise TypeError("no PolyRing generator for: %s"
                        % ", ".join(sorted(map(str, unknown))))

    p = sympy.Poly(expr, *_SYMPY_GENS)
    d = {}
    for mon, c in zip(p.monoms(), p.coeffs()):
        r = sympy.Rational(c)
        d[tuple(int(e) for e in mon)] = QQ(int(r.p)) / QQ(int(r.q))
    return PolyRing(d)


# BLAD differential polynomial -> PolyRing, DIRECTLY (no sympy round-trip).
# `_blad.read_terms` already yields (coeff, [(name, degree), ...]) -- exactly an
# exponent/coefficient dict -- so we build the PolyRing element in one O(terms)
# pass, grouping like terms as we go.  This replaces the _elt_to_sympy ->
# _sympy_to_polyring detour, whose incremental `out += mon` was O(terms^2)
# (sympy re-sorts the whole Add on every addition) and blew up on the large
# remainders the cell_eqs reducer produces.  A bracketed derivative jet has no
# generator (it should have been eliminated by the reduction) and raises, the
# same guard _sympy_to_polyring gave.
_GEN_IDX = {str(g): i for i, g in enumerate(PolyRing.gens())}
# derivative-jet generators are keyed by their BLAD name (`Psi[R1,R2]`), which
# is what read_terms yields; the mangled python-safe name never appears there.
for _b, _m in EXTRA_JET_PAIRS:
    _GEN_IDX[_b] = _GEN_IDX.pop(_m)


def _elt_to_polyring(e):
    n = PolyRing.ngens()
    d = {}
    for coeff, term in _blad.read_terms(e._h()):
        exps = [0] * n
        for nm, deg in term:
            i = _GEN_IDX.get(nm)
            if i is None:
                raise TypeError("no PolyRing generator for %r "
                                "(derivative jet survived reduction?)" % (nm,))
            exps[i] += int(deg)
        key = tuple(exps)
        d[key] = d.get(key, QQ(0)) + QQ(int(coeff))
    return PolyRing(d)


def forces_v_zero(P):
    """True iff the prime forces the inner variable v == 0 identically (every
    inner-variable coefficient lies in P) -- i.e. the ansatz has collapsed to a
    constant and the 'solution' is degenerate."""
    return all(g in P for g in V_PARAM_GENS)


def forces_psi_zero(P):
    """True iff the prime forces Psi == 0 via the amplitude collapsing (product
    ansatz Psi = A*F with every A-coefficient in P).  This is a trivial solution
    the per-cell Psi-reduction check misses, because A==0 is imposed by the
    variety, not by the cell."""
    return bool(AMP_PARAM_GENS) and all(g in P for g in AMP_PARAM_GENS)


def build_system_of_equations(eqn, constants):
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


# --- differential-Thomas decomposition of the ANSATZ ALONE ----------------
print("\nComputing native DifferentialThomas decomposition of the ansatz "
      "(%d ansatz + %d constancy eqs) ..." % (len(ansatz0), len(pconst)), flush=True)
_t0 = time.time()
cells_ds = dt.differential_thomas_decomposition(ansatz0 + pconst, [], prob['rk'])
_wall = time.time() - _t0
print("-> %d cells in %.1fs\n" % (len(cells_ds), _wall) + "=" * 72, flush=True)


def cell_eqs(ds):
    return list(dt.differential_system_equations(ds))


def cell_ineqs(ds):
    return list(dt.differential_system_inequations(ds))


for i, ds in enumerate(cells_ds, 1):
    leaders = [e.leader() for e in cell_eqs(ds) if not e.is_zero()]
    jl = [L for L in leaders if L is not None]
    if len(jl) != len(set(jl)):
        print("  SOUNDNESS FAIL cell %d: repeated leader" % i, flush=True)

try:
    with open(CELLS_OUT, 'w') as fh:
        for i, ds in enumerate(cells_ds, 1):
            fh.write("--- cell %d ---\nEQS: %s\nINEQS: %s\n\n"
                     % (i, [str(e) for e in cell_eqs(ds)],
                        [str(q) for q in cell_ineqs(ds)]))
    print("Wrote raw cells to", CELLS_OUT, flush=True)
except Exception as ex:
    print("(could not write cells file: %s)" % ex, flush=True)

if DECOMPOSE_ONLY:
    print("\n--decompose-only: stopping before the prime pipeline.", flush=True)
    sys.exit(0)


def adapt_cell(ds):
    param_eqs, param_ineqs, jet_ineqs = [], [], []
    for e in cell_eqs(ds):
        if e.is_zero() or has_jet(e) or is_param_constancy(e):
            continue
        param_eqs.append(_elt_to_polyring(e))
    for q in cell_ineqs(ds):
        if has_jet(q):
            jet_ineqs.append(q)
        elif is_param_constancy(q):
            continue
        else:
            param_ineqs.append(_elt_to_polyring(q))
    return dict(param_eqs=param_eqs, param_ineqs=param_ineqs, jet_ineqs=jet_ineqs)


def specialize(param_eqs):
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
    return list(spec) + list(pconst)


def full_prem(p, reductors, max_passes=64):
    """Reduce to a FIXPOINT.  R.differential_prem makes a single pass over the
    reductor list; reducing a high derivative by one reductor can re-expose a
    lower derivative reducible by an EARLIER reductor (e.g. Psi[R1,R1] -> ... ->
    DPsi[R1] -> n0*Psi[R1], and Psi[R1] is the leader of an earlier chain rule
    already passed).  A single pass therefore leaves first-order jets un-reduced
    whenever the PDE has first-derivative terms -- invisible for hydrogen's pure
    Laplacian, wrong for helium's 2/Ri d/dRi terms.  Looping to a fixpoint gives
    the true normal form (and is a no-op once the remainder is fully reduced)."""
    r = p if isinstance(p, type(R.one())) else R(p)
    h = R.one()
    for _ in range(max_passes):
        r2, h2 = R.differential_prem(r, reductors)
        h = h * h2
        if r2 == r:
            return r2, h
        r = r2
    return r, h


def prime_key(P):
    return tuple(sorted(str(g) for g in P.gens()))


strata_cache = {}
union_primes = {}          # GENUINE nontrivial (v != 0) solution varieties
degenerate_primes = {}     # nontrivial but v == 0 (ansatz collapsed to a constant)
trivial_primes = {}        # Psi == 0 forced

_cells = cells_ds if MAX_CELLS <= 0 else cells_ds[:MAX_CELLS]

for num, ds in enumerate(_cells, 1):
    cp = adapt_cell(ds)
    Z = sorted((p for p in cp['param_eqs']), key=str)
    Zkey = tuple(map(str, Z))

    # Reduce the PDE against the cell's OWN differential-triangular equations
    # (`cell_eqs`, polynomial form, initials carried as cofactors/inequations by
    # the Thomas decomposition) instead of a sympy-solved re-specialization of
    # the ansatz.  Dropping `specialize`/`sympy.solve` avoids its radicals,
    # RootOf objects, injected denominators, arbitrary branch choice, and
    # zero-substitution fallback.  Keyed on the cell's equations (not just the
    # parametric stratum Zkey), since the reduction now depends on the full cell.
    ce = cell_eqs(ds)
    cache_key = tuple(sorted(str(e) for e in ce))
    if cache_key not in strata_cache:
        # pconst dropped: cell_eqs already carries the (triangularized) constancy
        # relations, so `+ pconst` was redundant reductors (extra per-pass cost).
        reductors = list(ce)
        # Flushed phase markers with timings, so a stall is diagnosable from the
        # LAST line: stuck after "entering full_prem" => in the pseudo-reduction;
        # stuck after "entering GTZ" => in minimal_associated_primes (primdec).
        print("  [cell %d] entering full_prem: %d reductors ..." % (num, len(reductors)),
              flush=True)
        _t = time.time()
        rem_elt, h = full_prem(PDE, reductors)
        psi_rem_elt, _ = full_prem(R('Psi'), reductors)
        t_prem = time.time() - _t
        trivial = psi_rem_elt.is_zero()
        rem = _elt_to_polyring(rem_elt)
        if rem.is_zero():
            eqns = ()
        else:
            eqns = build_system_of_equations(rem, PolyRing_constants)
        gens = list(eqns) + list(Z)
        I = ideal(gens) if gens else ideal(PolyRing.zero())
        print("  [cell %d] full_prem %.1fs (%d eqns); entering GTZ minimal_associated_primes"
              " (%d gens) ..." % (num, t_prem, len(eqns), len(gens)), flush=True)
        _t = time.time()
        primes = I.minimal_associated_primes()
        t_gtz = time.time() - _t
        print("  [cell %d] GTZ %.1fs -> %d primes" % (num, t_gtz, len(primes)), flush=True)
        strata_cache[cache_key] = dict(spec_len=len(reductors), rem=rem, eqns=eqns,
                                       primes=primes, trivial=trivial)

    sc = strata_cache[cache_key]
    ineq_polys = list(cp['param_ineqs'])
    survivors = []
    for P in sc['primes']:
        if P.is_one():
            continue
        if any(g in P for g in ineq_polys):
            continue
        survivors.append(P)

    tag = "TRIVIAL (Psi==0 forced)" if sc['trivial'] else "nontrivial"
    print("\n--- cell %d: zero {%s}; ansatz %d eqs; %d param-ineqs, %d jet-ineqs; %s ---"
          % (num, ', '.join(Zkey) or '(none, generic)', sc['spec_len'],
             len(cp['param_ineqs']), len(cp['jet_ineqs']), tag), flush=True)
    if VERBOSE_REM:
        print("  remainder:", sc['rem'], flush=True)
    if sc['rem'].is_zero() and not sc['trivial']:
        print("  PDE reduces to 0: the whole stratum solves the PDE (nontrivially)", flush=True)
    if not survivors:
        print("  surviving solution varieties: NONE (all pruned / empty)", flush=True)
    for P in survivors:
        triv = sc['trivial'] or forces_psi_zero(P)
        deg = forces_v_zero(P)
        label = ("  [TRIVIAL: Psi=0]" if triv else
                 "  [DEGENERATE: v=0]" if deg else "  [GENUINE: v!=0]")
        print("   V:", P, label, flush=True)
        bucket = (trivial_primes if triv
                  else degenerate_primes if deg else union_primes)
        bucket.setdefault(prime_key(P), (P, []))[1].append(num)


def prune_enclosed(d):
    r"""
    Keep only the maximal (enclosing) varieties in a bucket of solution primes.

    A prime `P_i` is dropped when some other prime `P_j` defines a larger
    variety, `V(P_i) \subseteq V(P_j)`.  For the radical (prime) ideals here
    that containment is exactly the reverse ideal containment
    `P_j \subseteq P_i`, tested with Sage's native ``P_j <= P_i``.  Equal ideals
    -- mutual containment, e.g. the same variety surfacing from two different
    cells -- are a tie, broken by keeping the lowest-indexed entry so duplicates
    collapse to a single one.  Each dropped variety's source-cell list is folded
    into the entry that encloses it, so no cell provenance is lost.

    INPUT:

    - ``d`` -- a bucket ``{key: (P, cells)}`` mapping a prime's
      :func:`prime_key` to the prime ideal ``P`` and the list of cell numbers
      that produced it (the shape built for ``union_primes`` and its siblings)

    OUTPUT:

    a pair ``(kept, dropped)``, where ``kept`` is a bucket of the same shape
    holding only the enclosing varieties (with absorbed cells merged in) and
    ``dropped`` is the list of ``(dropped_key, encloser_key)`` pairs

    EXAMPLES:

    The `d_0 = d_1 = 0` wall of an `E = -1/2` component is contained in the full
    component, so it is pruned and its cell folded into the encloser::

        sage: R.<v4,a0,a1,b0,b1,c0,c1,E> = PolynomialRing(QQ)
        sage: wall = ideal(v4*c0 - b0, c1, b1, a1, a0, E + 1/2)
        sage: full = ideal(-b1*c0 + b0*c1, v4*c1 - b1, v4*c0 - b0, a1, a0, E + 1/2)
        sage: d = {prime_key(wall): (wall, [12]), prime_key(full): (full, [7])}
        sage: kept, dropped = prune_enclosed(d)
        sage: len(kept)
        1
        sage: P, cells = next(iter(kept.values()))
        sage: P == full
        True
        sage: cells
        [7, 12]
    """
    entries = [(k, P, list(cells)) for k, (P, cells) in d.items()]
    n = len(entries)
    encloser = [None] * n
    for i in range(n):
        Pi = entries[i][1]
        for j in range(n):
            if i == j:
                continue
            Pj = entries[j][1]
            if Pj <= Pi:                        # ideal Pj ⊆ Pi ⟺ V(Pi) ⊆ V(Pj)
                if not (Pi <= Pj) or j < i:      # strict, else tie -> low index
                    encloser[i] = j
                    break
    kept = {k: (P, list(cells)) for i, (k, P, cells) in enumerate(entries)
            if encloser[i] is None}
    dropped = []
    for i, (k, P, cells) in enumerate(entries):
        if encloser[i] is None:
            continue
        j = encloser[i]
        while encloser[j] is not None:                   # walk to a kept ancestor
            j = encloser[j]
        ek = entries[j][0]
        Pk, ck = kept[ek]
        kept[ek] = (Pk, sorted(set(ck) | set(cells)))
        dropped.append((k, ek))
    return kept, dropped


def dump_union(title, d, prune=False):
    r"""
    Print a bucket of solution varieties under ``title``, one prime per line.

    INPUT:

    - ``title`` -- a heading string printed above the list

    - ``d`` -- a bucket ``{key: (P, cells)}`` (see :func:`prune_enclosed`)

    - ``prune`` -- boolean (default: ``False``); when ``True`` and the global
      ``PRUNE_ENCLOSED`` is set, first drop the varieties contained in a larger
      one (via :func:`prune_enclosed`) so that only the enclosing varieties are
      printed, with a count of how many sub-varieties were pruned

    OUTPUT:

    the bucket actually printed -- the pruned one when pruning applies, otherwise
    ``d`` unchanged -- so the caller can report its size (e.g. in the verdict)

    """
    dropped = []
    if prune and PRUNE_ENCLOSED:
        d, dropped = prune_enclosed(d)
    note = ("" if not dropped
            else ", %d enclosed sub-variety(ies) pruned" % len(dropped))
    print("\n%s (%d distinct primes%s):\n" % (title, len(d), note))
    for key, (P, cells_for) in sorted(d.items(), key=lambda kv: str(kv[0])):
        print("  V:", P, "  (from cells:",
              ", ".join(map(str, sorted(set(cells_for)))) + ")")
    return d


print("\n" + "=" * 72)
genuine_primes = dump_union(
          "GENUINE solution varieties over all cells (v != 0 -- the real union)",
          union_primes, prune=True)
print("\n" + "-" * 72)
dump_union("DEGENERATE strata (nontrivial but v == 0, ansatz collapsed)",
          degenerate_primes)
print("\n" + "-" * 72)
dump_union("TRIVIAL (Psi==0) strata over all cells", trivial_primes)

print("\n" + "=" * 72)
if genuine_primes:
    print("VERDICT: %d GENUINE solution variety(ies) found for %s / ansatz %s."
          % (len(genuine_primes), PDE_NAME, ANSATZ))
else:
    print("VERDICT: NO genuine solution for %s / ansatz %s "
          "(every stratum is degenerate v=0 or trivial Psi=0)."
          % (PDE_NAME, ANSATZ))
