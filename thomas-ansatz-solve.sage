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
CELLS_OUT = _argval('--cells-out',
                    os.path.expanduser('~/thomas-experiments/%s_ansatz%s.cells'
                                       % (PDE_NAME, ANSATZ)))
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
prob = build_problem(PDE_NAME, ANSATZ)

R = prob['R']
COORDS = prob['coords']
ROOT_NAMES = [rn for rn, _ in prob['roots']]
PARAMS = prob['params']
JET_TOWER = prob['tower']                    # [Psi, DPsi(, DDPsi)]
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

_IB = {nm: sympy.IndexedBase(nm) for nm in (JET_TOWER + ['v'] + ROOT_NAMES)}
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


PolyRing = PolynomialRing(QQ, names=(COORDS + prob['jets'] + ['E'] + PARAMS))
PolyRing_constants = list(map(PolyRing, [str(c) for c in constants]))


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
        param_eqs.append(_elt_to_sympy(e))
    for q in cell_ineqs(ds):
        if has_jet(q):
            jet_ineqs.append(q)
        elif is_param_constancy(q):
            continue
        else:
            param_ineqs.append(_elt_to_sympy(q))
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
union_primes = {}
trivial_primes = {}

_cells = cells_ds if MAX_CELLS <= 0 else cells_ds[:MAX_CELLS]

for num, ds in enumerate(_cells, 1):
    cp = adapt_cell(ds)
    Z = sorted((p for p in cp['param_eqs']), key=str)
    Zkey = tuple(map(str, Z))

    if Zkey not in strata_cache:
        sub, spec = specialize(cp['param_eqs'])
        reductors = reductors_for(spec)
        rem_elt, h = full_prem(PDE, reductors)
        psi_rem_elt, _ = full_prem(R('Psi'), reductors)
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
    print("\n--- cell %d: zero {%s}; ansatz %d eqs; %d param-ineqs, %d jet-ineqs; %s ---"
          % (num, ', '.join(Zkey) or '(none, generic)', sc['spec_len'],
             len(cp['param_ineqs']), len(cp['jet_ineqs']), tag))
    if VERBOSE_REM:
        print("  remainder:", sc['rem'])
    if sc['rem'] == 0 and not sc['trivial']:
        print("  PDE reduces to 0: the whole stratum solves the PDE (nontrivially)")
    if not survivors:
        print("  surviving solution varieties: NONE (all pruned / empty)")
    bucket = trivial_primes if sc['trivial'] else union_primes
    for P in survivors:
        print("   V:", P)
        bucket.setdefault(prime_key(P), (P, []))[1].append(num)


def dump_union(title, d):
    print("\n%s (%d distinct primes):\n" % (title, len(d)))
    for key, (P, cells_for) in sorted(d.items(), key=lambda kv: str(kv[0])):
        print("  V:", P, "  (from cells:",
              ", ".join(map(str, sorted(set(cells_for)))) + ")")


print("\n" + "=" * 72)
dump_union("NONTRIVIAL solution varieties over all cells (the union)", union_primes)
print("\n" + "-" * 72)
dump_union("TRIVIAL (Psi==0) strata over all cells", trivial_primes)
