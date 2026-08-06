# -*- mode: python -*-
#
# thomas-ansatz-solve.sage
# ------------------------------------------------------------------------
# Solve an (ansatz, PDE) pair by the NewMethod staged pipeline of
# prolongation-projection-algorithm.tex (Phase III-forall / membership, run
# per Thomas cell -- the "staged route" of section 7):
#
#   1. pick one ansatz + one PDE system             (ansatz-library.sage)
#   2. differential-Thomas-decompose the ANSATZ ALONE  -> disjoint cells
#   3. reduce EACH PDE of the system modulo each cell  (differential
#      pseudo-remainder against the same cell)
#   4. forall-project each remainder onto the constants: collect like terms in
#      the independents + parametric jets, zero the constant coefficients,
#      combine the equations of all the remainders into ONE system, take
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
#   sage thomas-ansatz-solve.sage --pde hydrogen --ansatz 5 --latex
#   sage thomas-ansatz-solve.sage --pde navier-stokes --ansatz 25.3 --generic-cell
#
# --generic-cell skips step 2 entirely and runs steps 3-6 on a single cell built
# by hand from the ansatz -- see GenericCell below for what that cell is and
# what it is not.
#
# Author: Brent Baccala (AI assistant: Claude).  July 2026.

import hashlib, os, re, shutil, signal, subprocess, sys, time

USAGE = r"""usage: sage thomas-ansatz-solve.sage [--pde NAME] [--ansatz N] [options]

Solve an (ansatz, PDE) pair by the staged NewMethod pipeline: differential-
Thomas-decompose the ansatz alone into disjoint cells, reduce each PDE of the
system modulo each cell, project the remainders onto the constants, take
minimal associated primes, prune by the cell's inequations, and print the union
of the surviving solution varieties, classified GENUINE / DEGENERATE / TRIVIAL.

Choosing the problem
  --pde NAME         hydrogen | helium | navier-stokes | navier-stokes-nd
                     (default: hydrogen).  navier-stokes-nd is the
                     nondimensional rho = mu = 1 variant of navier-stokes,
                     which carries both as free constants the locus may
                     constrain.  coordinate_system() in ansatz-library.sage is
                     the authority on this list.
  --ansatz N         which ansatz (default: 5).  Integers and decimal variants
                     both work: 5, 20.1, 25.3, 25.34.  The catalogue, with a
                     comment block per family, is ansatz-library.sage.

Choosing what to compute
  --generic-cell     skip the Thomas decomposition and run the pipeline on the
                     ansatz's own generic cell: the ansatz equations plus the
                     constancy relations, with the initials, separants and
                     discriminants carried as inequations.  This is the way to
                     get a membership locus out of an ansatz whose
                     decomposition does not terminate -- but it covers only the
                     stratum where no initial vanishes, and the other cells can
                     hold solutions of their own.  See the GenericCell
                     docstring for what it does and does not entitle you to
                     conclude.
  --decompose-only   stop after the decomposition, before the prime pipeline.
                     Cheap way to find out whether a decomposition terminates
                     at all, and in what memory.
  --ranking NAME     orderly (default) or elimination.  elimination is a block
                     ranking -- each jet its own block, so a high jet like
                     DDPsi outranks every lower-jet derivative.  It changes the
                     decomposition, and is much the more expensive of the two.
  --max-cells N      process only the first N cells (default 0 = all).

Running the prime step
  --gtz-subprocess   run minimal_associated_primes in a standalone Singular
                     subprocess with option(prot) instead of in-process
                     libsingular.  The in-process call is silent -- primdec.lib's
                     minAssGTZ carries no dbprint instrumentation, and
                     option(prot) does not reach stdout through Sage's
                     libsingular wrapper -- so a GTZ call that runs for hours
                     prints nothing at all.  Under this flag you get: a protocol
                     log per ideal (degrees, pair counts, sub-algorithm names,
                     memory) you can tail while it runs; a killable child
                     process, optionally bounded by --gtz-timeout; and the
                     result cached on disk, so a re-run after an interruption
                     skips every GTZ call that already finished.
  --gtz-dir PATH     where the Singular input, protocol log and cached result
                     go, one triple per ideal (default:
                     ~/thomas-experiments/gtz/<pde>_ansatz<N>[_generic]/).
  --gtz-timeout SECS abandon a GTZ call after SECS and abort with the log path
                     (default 0 = no limit).  Only meaningful with
                     --gtz-subprocess.
  --singular-bin PATH  the Singular to run (default: the first on PATH, else
                     <sys.prefix>/bin/Singular).

Output
  --verbose-remainder  print each PDE's remainder for every cell.
  --keep-enclosed    keep GENUINE varieties contained in another genuine one.
                     The default prints only the maximal (enclosing) ones,
                     since a smaller prime surfacing from a second cell is
                     usually the same solution family seen again.
  --latex            re-print the genuine union as a LaTeX subequations block
                     in the form the paper uses, denominators cleared.
  --cells-out PATH   where to write the raw cells (default:
                     ~/thomas-experiments/<pde>_ansatz<N>_<ranking>.cells, with
                     _generic appended under --generic-cell).
  --help, -h         this message.

Examples
  sage thomas-ansatz-solve.sage --pde hydrogen --ansatz 5
  sage thomas-ansatz-solve.sage --pde hydrogen --ansatz 5 --latex
  sage thomas-ansatz-solve.sage --pde navier-stokes-nd --ansatz 25.34
  sage thomas-ansatz-solve.sage --pde navier-stokes --ansatz 25.3 --generic-cell
  sage thomas-ansatz-solve.sage --pde helium --ansatz 15 --decompose-only
"""

if '--help' in sys.argv or '-h' in sys.argv:
    sys.stdout.write(USAGE)
    sys.stdout.flush()
    # os._exit, not sys.exit: the sage runner swallows SystemExit and returns 1
    # regardless of its code (the same quirk --decompose-only exits 1 under), so
    # a plain `sys.exit(0)` here would make `--help` look like a failure.
    os._exit(0)

import sympy

_T_START = time.time()

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
# Skip the differential Thomas decomposition and run the downstream pipeline on
# the ansatz's own generic cell (see GenericCell).  For an ansatz whose
# decomposition does not terminate, this is the only way to get a membership
# locus at all; for one whose does, it computes the single cell the
# decomposition would have called generic.
GENERIC_CELL = '--generic-cell' in sys.argv
# By default, prune GENUINE varieties whose zero-set is contained in another
# genuine variety's, printing only the maximal (enclosing) ones.  Two "genuine"
# primes coming out of different cells are often nested (e.g. the c1=0 wall of a
# larger E=-1/2 component reappears as its own prime), which double-counts one
# solution family.  --keep-enclosed restores the raw per-cell union.
PRUNE_ENCLOSED = '--keep-enclosed' not in sys.argv
# Re-print the final GENUINE union as a LaTeX `subequations`/`align` block, in
# the form the paper uses (one \left(...\right) per prime, each \label'd
# `ideal:N`).  Coefficients are cleared of denominators first, so `a1 - 1/2*b0`
# prints as `2 a_{1} - b_{0}`.
LATEX_OUT = '--latex' in sys.argv
# Run the minimal-associated-primes step in a standalone Singular subprocess
# (option(prot), killable, cached) rather than in-process through libsingular.
# See minimal_associated_primes_gtz for why the in-process call cannot be made
# to talk.
GTZ_SUBPROCESS = '--gtz-subprocess' in sys.argv
GTZ_TIMEOUT = int(_argval('--gtz-timeout', '0'))
GTZ_DIR = _argval('--gtz-dir',
                  os.path.expanduser('~/thomas-experiments/gtz/%s_ansatz%s%s'
                                     % (PDE_NAME, ANSATZ,
                                        '_generic' if GENERIC_CELL else '')))
SINGULAR_BIN = _argval('--singular-bin',
                       shutil.which('Singular')
                       or os.path.join(sys.prefix, 'bin', 'Singular'))
CELLS_OUT = _argval('--cells-out',
                    os.path.expanduser('~/thomas-experiments/%s_ansatz%s_%s%s.cells'
                                       % (PDE_NAME, ANSATZ, RANKING,
                                          '_generic' if GENERIC_CELL else '')))
os.makedirs(os.path.dirname(CELLS_OUT), exist_ok=True)
if GTZ_SUBPROCESS:
    os.makedirs(GTZ_DIR, exist_ok=True)


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
PDES = prob['pdes']                          # the PDE system (list; often 1)
PDE_PARAMS = prob['pde_params']              # the PDE's own constants (E; rho,mu)


# --- consistent variable rendering ----------------------------------------
# Three notations reach the output otherwise: `str()` on a BLAD differential
# polynomial writes a derivative jet with UNDERSCORES (`Psi_x_x`), the ansatz
# strings that built the ring use DifferentialAlgebra's BRACKET input notation
# (`Psi[x,x]`), and a bare python list of names prints QUOTED.  Print bracket
# notation everywhere -- it is what the ring parser accepts, so the output can
# be fed back in.  (The `elimination` ranking's mangled PolyRing generator names
# are underscore-form too, so they normalise through the same function.)
_COORD_SET = set(COORDS)


def to_bracket(s):
    r"""
    Rewrite underscore-form derivative jets in ``s`` as bracket-form ones.

    A token is rewritten only when every underscore-separated piece after the
    head is a coordinate name, so ordinary names pass through untouched and a
    string already in bracket form is a fixed point.

    INPUT:

    - ``s`` -- anything with a ``str()`` (a BLAD element, a polynomial, a string)

    OUTPUT: a string

    EXAMPLES::

        sage: to_bracket('-Psi_x_x*r - 2*Psi*r*E')      # not tested (needs COORDS)
        '-Psi[x,x]*r - 2*Psi*r*E'
    """
    def repl(m):
        head, *idx = m.group(0).split('_')
        if idx and head and all(p in _COORD_SET for p in idx):
            return '%s[%s]' % (head, ','.join(idx))
        return m.group(0)
    return re.sub(r'[A-Za-z]\w*', repl, str(s))


for _P in PDES:
    print("PDE:", to_bracket(_P))
print("ansatz (%d eqs):" % len(ansatz0))
for s in prob['ansatz_eqs_str']:
    print("   ", to_bracket(s))
print("params:", ", ".join(PARAMS))
# Excluded-locus header: ONLY for ansatz variants that record one (the
# 25.3x normalization rungs).  Emitting it unconditionally would perturb
# the byte-identity regression gate on hydrogen/5 and helium/9.
if prob.get('excludes'):
    print("excluded locus (results lift to the parent ansatz off this locus):")
    for _x in prob['excludes']:
        print("   ", _x)


# ==========================================================================
# native <-> sympy bridge (dynamic names), then the DT pipeline
# ==========================================================================
constants = [sympy.Symbol(p) for p in (PDE_PARAMS + PARAMS)]
strata_constants = [sympy.Symbol(p) for p in PARAMS]

_IB = {nm: sympy.IndexedBase(nm) for nm in JETS}
_DERIV = {c: sympy.Symbol(c) for c in COORDS}
_PARAM = {p: sympy.Symbol(p) for p in (PDE_PARAMS + PARAMS)}
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
                                     + PDE_PARAMS + PARAMS))
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
               + [_PARAM[p] for p in PDE_PARAMS]
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
    constant and the 'solution' is degenerate.  When the ansatz declares NO
    v_params (e.g. the normalization rungs 25.32+ pin v5 = 1, so the amplitude
    vector cannot vanish), the test does not apply and nothing is degenerate --
    without the bool() guard, all() on the empty list is vacuously True and
    every nontrivial prime would be mislabeled DEGENERATE."""
    return bool(V_PARAM_GENS) and all(g in P for g in V_PARAM_GENS)


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


# --- cells ----------------------------------------------------------------
# Everything downstream touches a cell only through these two functions, so a
# cell can be anything that answers them -- a differential system from the
# decomposition, or the hand-built GenericCell below.
class GenericCell(object):
    r"""
    The ansatz's own generic cell, built without decomposing anything.

    Its equations are the ansatz equations together with the parameter-constancy
    relations, and its inequations are the nondegeneracy conditions that make
    that system simple: the initial and the separant of every equation, plus the
    discriminant of every equation of degree `\geq 2` in its leader.  That is the
    stratum on which no initial, separant or discriminant vanishes -- the one
    branch the Thomas decomposition never has to split, and the one it would
    label generic.

    The inequations are the point of building it this way.  Reducing the PDE
    against the bare ansatz and taking minimal associated primes returns
    varieties that are only there because some initial vanishes on them -- the
    ODE degenerating to lower order, say -- and those are artifacts of ignoring
    the case distinction, not solutions.  Carrying the inequations lets the
    downstream pruning drop them, exactly as it does for a real cell.

    Two things it is NOT:

    - **Not the whole answer.**  The decomposition's other cells are where the
      initials do vanish, and they can carry solutions of their own (that is how
      ansatz 25 reaches its viscous solutions, on `a_0 = a_1 = 0`).  A
      `--generic-cell` run says nothing about them.

    - **Not checked for coherence, and possibly not passive.**  The
      decomposition completes each cell to a passive system, adjoining the
      integrability conditions between equations whose leaders are derivatives
      of the same dependent variable.  This cell takes the ansatz as given and
      adjoins nothing.  **Nothing here tests passivity**: the repeated-leader
      check below establishes only that the system is triangular, which is a
      weaker property entirely -- a triangular system can still have an
      integrability condition that does not reduce to zero.

      When the cell is not passive, `full_prem` against it is not a normal
      form: a PDE that lies in the completed cell's differential ideal can
      still leave a nonzero remainder here.  The projection then reads that
      remainder's coefficients as constraints, so the membership locus comes
      out too SMALL -- `--generic-cell` under-reports, and silently.

      Verified passive so far: `navier-stokes` / 25.34 (66 integrability
      conditions, all reducing to zero).  Not yet checked: 25, 25.2, 25.3 --
      the check is itself expensive on those.  Until it is run, treat their
      loci as lower bounds.
    """
    def __init__(self, eqs, ineqs):
        self.eqs = list(eqs)
        self.ineqs = list(ineqs)


def cell_eqs(ds):
    if isinstance(ds, GenericCell):
        return list(ds.eqs)
    return list(dt.differential_system_equations(ds))


def cell_ineqs(ds):
    if isinstance(ds, GenericCell):
        return list(ds.ineqs)
    return list(dt.differential_system_inequations(ds))


def _discriminant_in_leader(e):
    r"""
    The discriminant of ``e`` with respect to its own leader, or ``None``.

    ``None`` is returned when the discriminant is not a condition: a constant
    equation, degree `1` in the leader (the discriminant is then a unit, and the
    separant already says what the initial does not), or an equation carrying a
    derivative jet, which has no image in ``PolyRing``.  The one place degree
    `\geq 2` in the leader occurs in this library is the nonlinear ODE of ansatz
    25, whose jets are all of order `0` -- so the restriction costs nothing
    there.
    """
    L = e.leader()
    if L is None or '[' in L:
        return None
    for _coeff, term in _blad.read_terms(e._h()):
        if any('[' in nm for nm, _deg in term):
            return None                  # a derivative jet: no PolyRing image
    p = _elt_to_polyring(e)
    g = PolyRing(L)
    if p.degree(g) < 2:
        return None
    d = p.discriminant(g)
    return None if d.is_zero() or d.is_unit() else R(str(d))


def generic_cell(eqs):
    """Build the :class:`GenericCell` of the equations ``eqs``."""
    keep = [e for e in eqs if not e.is_zero()]
    ineqs, seen = [], set()

    def add(q):
        # leader() is None for a rational constant: `1 != 0` says nothing.
        if q is None or q.is_zero() or q.leader() is None:
            return
        if str(q) not in seen:
            seen.add(str(q))
            ineqs.append(q)

    for e in keep:
        add(e.initial())
        add(e.separant())
        add(_discriminant_in_leader(e))
    return GenericCell(keep, ineqs)


# --- differential-Thomas decomposition of the ANSATZ ALONE ----------------
if GENERIC_CELL:
    print("\n--generic-cell: skipping the differential Thomas decomposition; "
          "building the ansatz's generic cell (%d ansatz + %d constancy eqs) ..."
          % (len(ansatz0), len(pconst)), flush=True)
    cells_ds = [generic_cell(list(ansatz0) + list(pconst))]
    print("-> 1 generic cell, %d inequations:" % len(cell_ineqs(cells_ds[0])),
          flush=True)
    for _q in cell_ineqs(cells_ds[0]):
        print("     ", to_bracket(_q), "!= 0", flush=True)
    print("   (the other cells of the decomposition -- where these vanish -- "
          "are NOT computed)\n" + "=" * 72, flush=True)
else:
    print("\nComputing native DifferentialThomas decomposition of the ansatz "
          "(%d ansatz + %d constancy eqs) ..." % (len(ansatz0), len(pconst)),
          flush=True)
    _t0 = time.time()
    cells_ds = dt.differential_thomas_decomposition(ansatz0 + pconst, [], prob['rk'])
    _wall = time.time() - _t0
    print("-> %d cells in %.1fs\n" % (len(cells_ds), _wall) + "=" * 72, flush=True)


for i, ds in enumerate(cells_ds, 1):
    leaders = [e.leader() for e in cell_eqs(ds) if not e.is_zero()]
    jl = [L for L in leaders if L is not None]
    if len(jl) != len(set(jl)):
        print("  SOUNDNESS FAIL cell %d: repeated leader" % i, flush=True)

try:
    with open(CELLS_OUT, 'w') as fh:
        for i, ds in enumerate(cells_ds, 1):
            fh.write("--- cell %d ---\nEQS: %s\nINEQS: %s\n\n"
                     % (i, [to_bracket(e) for e in cell_eqs(ds)],
                        [to_bracket(q) for q in cell_ineqs(ds)]))
    print("Wrote raw cells to", CELLS_OUT, flush=True)
except Exception as ex:
    print("(could not write cells file: %s)" % ex, flush=True)

def print_total_time():
    t = time.time() - _T_START
    print("Total time: %.1fs (%d:%02d:%05.2f)"
          % (t, int(t) // 3600, (int(t) % 3600) // 60, t % 60), flush=True)


if DECOMPOSE_ONLY:
    print("\n--decompose-only: stopping before the prime pipeline.", flush=True)
    print_total_time()
    sys.exit(0)


def ineq_coeff_set(q):
    r"""
    The constant coefficients of an inequation ``q``, as a tuple, or ``None``.

    An inequation is a condition on *functions*: `q \neq 0` fails exactly where
    `q` vanishes identically, which is where every coefficient of `q` on the jet
    monomials vanishes.  Collecting those coefficients turns the inequation into
    the condition on the constants that :func:`prune` can test -- so a jet-
    carrying inequation such as `a_0 + a_1 s \neq 0` prunes the variety
    `a_0 = a_1 = 0` instead of being ignored for having a jet in it.

    An inequation in the constants alone collects to ``(q,)``, so the test
    ``all(g in P for g in ineq_coeff_set(q))`` reduces on those to the plain
    ``q in P`` it generalises.

    ``None`` is returned when ``q`` has no image in ``PolyRing`` -- a derivative
    jet survived in it -- and such an inequation prunes nothing, as before.
    """
    try:
        qp = _elt_to_polyring(q)
    except TypeError:
        return None
    if qp.is_zero():
        return None
    return tuple(build_system_of_equations(qp, PolyRing_constants))


def adapt_cell(ds):
    param_eqs, param_ineqs, jet_ineqs, ineq_coeffs = [], [], [], []
    for e in cell_eqs(ds):
        if e.is_zero() or has_jet(e) or is_param_constancy(e):
            continue
        param_eqs.append(_elt_to_polyring(e))
    for q in cell_ineqs(ds):
        if is_param_constancy(q):
            continue
        cs = ineq_coeff_set(q)
        if cs is not None:
            ineq_coeffs.append(cs)
        if has_jet(q):
            jet_ineqs.append(q)
        else:
            param_ineqs.append(_elt_to_polyring(q))
    return dict(param_eqs=param_eqs, param_ineqs=param_ineqs,
                jet_ineqs=jet_ineqs, ineq_coeffs=ineq_coeffs)


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


# --- the minimal-associated-primes step -----------------------------------
#
# I.minimal_associated_primes() calls primdec.lib's minAssGTZ through
# libsingular, in this process.  Cheap to call, and completely opaque: the
# minAssGTZ -> minAssGTZ_i chain carries no dbprint instrumentation (raising
# printlevel adds nothing), and option(prot) does not reach stdout through
# Sage's libsingular wrapper -- verified both ways, opt['prot'] = True and
# opt_ctx(prot=True), on minAssGTZ and on a plain cyclic-6 std: no output
# either way, while the same computation in a standalone Singular prints the
# full protocol stream.  So an in-process GTZ call that runs for hours prints
# nothing whatsoever.  That is how the 4 Aug navier-stokes / ansatz 25
# --generic-cell run died: 32 hours of CPU, RSS flat at 10.34 GB, killed by a
# power cut, and the only record in the log was "entering GTZ".
#
# --gtz-subprocess runs the same primdec.lib minAssGTZ in a standalone Singular
# with option(prot) and option(mem), and buys three things the in-process call
# cannot give: a live protocol log per ideal, a child process that is killable
# and can be bounded by --gtz-timeout, and a result cached on disk so an
# interrupted run resumes instead of recomputing.

_GTZ_BEGIN, _GTZ_PRIME, _GTZ_END = '===PRIMES_BEGIN===', '===PRIME===', '===PRIMES_END==='


def _gtz_key(I):
    """Content hash of an ideal: its generators and its ring's variables."""
    h = hashlib.sha1()
    h.update(','.join(sorted(str(g) for g in I.gens())).encode())
    h.update(('|' + ','.join(str(v) for v in I.ring().gens())).encode())
    return h.hexdigest()[:16]


def _gtz_script(I, primes_path):
    """A standalone Singular program computing minAssGTZ(I) under option(prot).

    The primes go to their own file rather than to stdout, so the protocol
    stream and the result never have to be untangled from one another.
    """
    # Every name here is gtz-prefixed: primdec.lib and the libraries it pulls
    # in occupy a lot of the top-level namespace (a plain `res` collides), and
    # a collision is only reported as "identifier in use" after the expensive
    # part has already run.
    gens = [str(g) for g in I.gens() if not g.is_zero()] or ['0']
    return '\n'.join([
        'option(prot);',
        'option(mem);',
        'LIB "primdec.lib";',
        'ring gtzring = 0,(%s),dp;' % ','.join(str(v) for v in I.ring().gens()),
        'ideal gtzI =\n  %s;' % ',\n  '.join(gens),
        'list gtzP = minAssGTZ(gtzI);',
        'string gtzout = "%s";' % _GTZ_BEGIN,
        'for (int gtzi = 1; gtzi <= size(gtzP); gtzi++)',
        '{ gtzout = gtzout + newline + "%s" + newline + string(gtzP[gtzi]); }' % _GTZ_PRIME,
        'gtzout = gtzout + newline + "%s";' % _GTZ_END,
        'write(":w %s", gtzout);' % primes_path,
        'quit;',
        ''])


def _gtz_errors(log_path):
    """Singular's error lines, if any.  It exits 0 even after a hard error, so
    the log is the only place a failed run announces itself."""
    try:
        with open(log_path, errors='replace') as fh:
            return [ln.strip() for ln in fh if ln.lstrip().startswith('?')]
    except OSError:
        return []


def _gtz_parse(path, R):
    """Read back the primes file written by _gtz_script.

    Singular may fold a long `string(ideal)` across several lines, so a prime's
    body is rejoined before splitting on the generator commas (a polynomial
    never contains one).  Absence of the closing marker means the file is a
    torn write from an interrupted run -- treated as no result at all.
    """
    lines = open(path).read().splitlines()
    if not lines or lines[0].strip() != _GTZ_BEGIN or lines[-1].strip() != _GTZ_END:
        raise RuntimeError('%s: no %s marker (truncated?)' % (path, _GTZ_END))
    blocks, cur = [], None
    for line in lines[1:-1]:
        if line.strip() == _GTZ_PRIME:
            cur = []
            blocks.append(cur)
        elif cur is not None:
            cur.append(line.strip())
    return [R.ideal([R(g) for g in ''.join(b).split(',') if g.strip()]) for b in blocks]


def minimal_associated_primes_gtz(I, tag):
    """minAssGTZ(I), in-process by default or in a Singular subprocess.

    `tag` labels this ideal's files in --gtz-dir; it only has to be unique
    within a run, the cache key proper is the hash of the ideal itself.
    """
    if not GTZ_SUBPROCESS:
        return I.minimal_associated_primes()

    R = I.ring()
    if all(g.is_zero() for g in I.gens()):
        # The zero ideal of a domain is prime; minAssGTZ says so too, but
        # spawning a Singular to be told that is silly.
        return [R.ideal(R.zero())]

    base = os.path.join(GTZ_DIR, '%s_%s' % (tag, _gtz_key(I)))
    sing_path, log_path, primes_path = base + '.sing', base + '.log', base + '.primes'

    if os.path.exists(primes_path):
        try:
            primes = _gtz_parse(primes_path, R)
            print('  [%s] GTZ cached -> %d primes (%s)' % (tag, len(primes), primes_path),
                  flush=True)
            return primes
        except Exception as exc:
            print('  [%s] cached GTZ result unusable (%s); recomputing' % (tag, exc),
                  flush=True)

    with open(sing_path, 'w') as fh:
        fh.write(_gtz_script(I, primes_path))
    print('  [%s] GTZ in Singular subprocess; tail the protocol at %s' % (tag, log_path),
          flush=True)

    with open(log_path, 'wb') as log:
        # start_new_session: the child leads its own process group, so a
        # timeout can kill the whole group (Singular forks for factorization).
        proc = subprocess.Popen([SINGULAR_BIN, '-q', sing_path],
                                stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            rc = proc.wait(timeout=GTZ_TIMEOUT or None)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            raise RuntimeError('[%s] GTZ exceeded --gtz-timeout %ds; protocol log: %s'
                               % (tag, GTZ_TIMEOUT, log_path))
    errors = _gtz_errors(log_path)
    if rc != 0 or errors:
        raise RuntimeError('[%s] Singular failed (exit %d)%s; protocol log: %s'
                           % (tag, rc,
                              ''.join('\n    ' + e for e in errors[:5]), log_path))
    return _gtz_parse(primes_path, R)


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
        # Leading blank line so each cell's timing block is separated from the
        # previous cell's variety list (the "--- cell N ---" header supplies the
        # blank line before the result block).
        print("\n  [cell %d] entering full_prem: %d reductors ..." % (num, len(reductors)),
              flush=True)
        _t = time.time()
        # Reduce EACH PDE of the system against the same cell.  For the
        # membership locus every PDE must lie in the cell's differential ideal,
        # so the constant-coefficient equations of all the remainders are
        # combined into ONE system, solved once per cell.
        rem_elts = [full_prem(P_, reductors)[0] for P_ in PDES]
        psi_rem_elt, _ = full_prem(R('Psi'), reductors)
        t_prem = time.time() - _t
        trivial = psi_rem_elt.is_zero()
        rems = [_elt_to_polyring(re_) for re_ in rem_elts]
        nonzero = [r_ for r_ in rems if not r_.is_zero()]
        if not nonzero:
            eqns = ()
        elif len(nonzero) == 1:
            eqns = build_system_of_equations(nonzero[0], PolyRing_constants)
        else:
            # deduplicate the combined list the same way
            # build_system_of_equations already does internally.
            _all = []
            for r_ in nonzero:
                _all.extend(build_system_of_equations(r_, PolyRing_constants))
            eqns = tuple(set(_all))
        gens = list(eqns) + list(Z)
        I = ideal(gens) if gens else ideal(PolyRing.zero())
        print("  [cell %d] full_prem %.1fs (%d eqns); entering GTZ minimal_associated_primes"
              " (%d gens) ..." % (num, t_prem, len(eqns), len(gens)), flush=True)
        _t = time.time()
        primes = minimal_associated_primes_gtz(I, 'cell%d' % num)
        t_gtz = time.time() - _t
        print("  [cell %d] GTZ %.1fs -> %d primes" % (num, t_gtz, len(primes)), flush=True)
        strata_cache[cache_key] = dict(spec_len=len(reductors), rems=rems, eqns=eqns,
                                       primes=primes, trivial=trivial)

    sc = strata_cache[cache_key]
    # A prime is dropped when some inequation of the cell vanishes identically
    # on the whole of it -- every coefficient of the inequation lying in the
    # prime.  See ineq_coeff_set: on inequations in the constants alone this is
    # the `q in P` test it replaces.
    survivors = []
    for P in sc['primes']:
        if P.is_one():
            continue
        if any(all(g in P for g in cs) for cs in cp['ineq_coeffs']):
            continue
        survivors.append(P)

    tag = "TRIVIAL (Psi==0 forced)" if sc['trivial'] else "nontrivial"
    print("\n--- cell %d: zero {%s}; ansatz %d eqs; %d param-ineqs, %d jet-ineqs; %s ---"
          % (num, ', '.join(Zkey) or '(none, generic)', sc['spec_len'],
             len(cp['param_ineqs']), len(cp['jet_ineqs']), tag), flush=True)
    if VERBOSE_REM:
        for _i, _r in enumerate(sc['rems'], 1):
            print("  remainder[pde %d]:" % _i, to_bracket(_r), flush=True)
    if all(r_.is_zero() for r_ in sc['rems']) and not sc['trivial']:
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


def clear_denominators(g):
    r"""
    Scale a polynomial over `\QQ` by the lcm of its coefficient denominators.

    The paper writes generators with integer coefficients, so ``a1 - 1/2*b0``
    must appear as `2 a_{1} - b_{0}`.  Signs are left alone (a leading minus is
    part of the generator as GTZ returned it).

    EXAMPLES::

        sage: R.<a0,a1,b0> = PolynomialRing(QQ)
        sage: clear_denominators(a1 - 1/2*b0)
        2*a1 - b0
        sage: clear_denominators(R.zero())
        0
    """
    dens = [c.denominator() for c in g.coefficients()]
    return g * lcm(dens) if dens else g


def latex_union(d, label='ideal'):
    r"""
    Print a bucket of solution varieties as a LaTeX ``subequations`` block.

    Each prime becomes one ``align`` row, ``& \left(g_1, \ldots, g_k\right)``
    followed by ``\label{<label>:N}``, numbered ``1..n`` in the same order
    :func:`dump_union` prints them.  Coefficients are cleared of denominators
    by :func:`clear_denominators`.

    INPUT:

    - ``d`` -- a bucket ``{key: (P, cells)}`` (see :func:`prune_enclosed`) --
      normally the pruned genuine union returned by :func:`dump_union`

    - ``label`` -- string (default: ``'ideal'``); the ``\label`` prefix, also
      used as the block's own ``\label``

    OUTPUT: ``None`` (the block is printed)
    """
    items = sorted(d.items(), key=lambda kv: str(kv[0]))
    rows = []
    for i, (_key, (P, _cells)) in enumerate(items, 1):
        gens = ", ".join(latex(clear_denominators(g)) for g in P.gens())
        rows.append(r"& \left(%s\right)\label{%s:%d}" % (gens, label, i))
    print(r"\begin{subequations}")
    print(r"\label{%s}" % label)
    print(r"\begin{align}")
    print(" \\\\\n".join(rows))
    print(r"\end{align}")
    print(r"\end{subequations}")


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

if LATEX_OUT and genuine_primes:
    print("\n" + "-" * 72)
    print("GENUINE solution varieties, LaTeX (paper) form:\n")
    latex_union(genuine_primes)

print("\n" + "=" * 72)
print_total_time()
