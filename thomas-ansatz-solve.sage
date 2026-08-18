# -*- mode: python -*-
#
# thomas-ansatz-solve.sage
# ------------------------------------------------------------------------
# Solve an (ansatz, PDE) pair by NewMethod.tex's Algorithm MembershipLocus
# (paper section "Computing the Membership Locus"), run per Thomas component:
#
#   1. pick one ansatz + one PDE system.  BOTH may declare inequations -- the
#      ansatz its non-degeneracy forms, the system its P^!= -- and Q is the
#      product of all of them                        (ansatz-library.sage)
#   2. differential-Thomas-decompose the ansatz TOGETHER WITH Q as the
#      system's inequation, (A(c), {Q}) -- so Q participates in the Thomas
#      decomposition's own case-splitting, and every resulting component's
#      inequations already carry whatever of Q's non-degeneracy survives on
#      it -- into disjoint components
#   3. reduce EACH PDE of the system, unmultiplied, modulo each component's
#      own equations (differential pseudo-remainder).  There is no separate
#      guard step: Q's non-degeneracy was already imposed when the component
#      was built, in step 2, so a plain reduction is all the membership test
#      needs.
#   4. forall-project each remainder onto the constants: collect like terms in
#      the independents + parametric jets, zero the constant coefficients,
#      combine the equations of all the remainders into ONE system, take
#      minimal associated primes, prune by the component's own inequations
#      (which, since Q was fed into the decomposition, already excludes the
#      loci where Q degenerates -- no separate Q-tracking is needed)
#   5. union the surviving varieties over all components
#   6. print the union
#
# The ansatz and PDE are supplied by ansatz-library.sage in the differential-
# algebra formulation (differential-polynomial equations), NOT by helium.sage.
# The downstream (steps 3-6) is joca-thomas-native-dt.sage's proven pipeline,
# generalised to the problem's jet/param names.
#
# Validation: `--pde hydrogen --ansatz 5` must reproduce the hydrogen
# decomposition and known union of solution varieties under the CURRENT
# algorithm -- feeding Q into the decomposition (step 2) changes both the
# component count and, in general, which components exist, relative to the
# pre-2026-08 script that decomposed the ansatz alone and multiplied Q into
# each PDE reduction instead.  See the commit that introduced this change for
# the new baseline count; do not compare it against an old "29 cells" figure
# recorded before the paper's algorithm changed.
#
#   sage thomas-ansatz-solve.sage --pde hydrogen --ansatz 5 [--decompose-only]
#   sage thomas-ansatz-solve.sage --pde helium   --ansatz 5
#   sage thomas-ansatz-solve.sage --pde hydrogen --ansatz 5 --latex
#   sage thomas-ansatz-solve.sage --pde navier-stokes --ansatz 25.3 --generic-cell
#
# --generic-cell skips step 2's decomposition entirely and runs steps 3-6 on a
# single cell built by hand from the ansatz -- see GenericCell below for what
# that cell is and what it is not.  Q is folded into that cell's inequations
# too, for consistency with the decomposed path.
#
# Author: Brent Baccala (AI assistant: Claude).  July-August 2026.

import hashlib, itertools, os, re, shutil, signal, subprocess, sys, time

USAGE = r"""usage: sage thomas-ansatz-solve.sage [--pde NAME] [--ansatz N] [options]

Solve an (ansatz, PDE) pair by NewMethod.tex's Algorithm MembershipLocus:
differential-Thomas-decompose the ansatz together with Q (the product of all
declared inequations) into disjoint components, reduce each PDE of the system
modulo each component's own equations, project the remainders onto the
constants, take minimal associated primes, prune by the component's own
inequations (which already carry Q's contribution, since Q was fed into the
decomposition), and print the union of the surviving GENUINE solution
varieties.

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
  --cells-out PATH   write the raw cells to PATH.  Omitted, no cells file is
                     written; nothing reads one back, so it is a debugging
                     artifact rather than an output.
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
    r"""
    The command-line argument following ``flag``, or ``default``.

    A deliberately minimal stand-in for :mod:`argparse`: every option of this
    script is either a bare flag, tested directly with ``in sys.argv``, or a
    flag followed by exactly one value.  An unrecognized argument is ignored
    rather than rejected, so the script can be given extra flags by a wrapper.

    INPUT:

    - ``flag`` -- string; the option to look for, leading dashes included

    - ``default`` -- value returned when ``flag`` does not appear (default:
      ``None``)

    OUTPUT: the argument following ``flag`` as a string, or ``default``

    EXAMPLES::

        sage: _argval('--pde', 'hydrogen')       # not tested (reads sys.argv)
        'hydrogen'
    """
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
# The raw cells are a debugging artifact, not an output of the algorithm, and
# nothing reads the file back.  Written only when asked for: no --cells-out, no
# file.  (It used to default to a path built from the pde/ansatz/ranking, which
# also meant the --decompose-with-inequations and plain runs of the same ansatz
# collided on one filename.)
CELLS_OUT = _argval('--cells-out')
if CELLS_OUT:
    os.makedirs(os.path.dirname(os.path.abspath(CELLS_OUT)), exist_ok=True)
if GTZ_SUBPROCESS:
    os.makedirs(GTZ_DIR, exist_ok=True)


def patch_latex_varify():
    r"""
    Make Sage's LaTeX printer render the jet ``DPsi`` as `\Psi'`.

    Sage would otherwise typeset it as `\mathit{DPsi}`, which is how a variable
    named ``DPsi`` should print but not how the paper writes the first
    derivative of `\Psi`.  The patch is monkeyed onto
    :func:`sage.misc.latex.latex_varify` because that is the single point every
    ``latex()`` call routes variable names through; there is no per-ring hook.

    OUTPUT: ``None`` (the module is patched in place)

    EXAMPLES::

        sage: patch_latex_varify()               # not tested (global side effect)
    """
    from sage.misc.latex import latex_varify
    import sage.misc.latex
    original = latex_varify
    def custom(a, is_fname=False):
        r"""The patched :func:`latex_varify`: ``DPsi`` specially, else as before."""
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
PDE_INEQS = prob['pde_ineqs']                # its inequations (P^!=; ['Psi'] or [])
PDE_PARAMS = prob['pde_params']              # the PDE's own constants (E; rho,mu)

# Q = the product of ALL declared inequations (the ansatz's own non-degeneracy
# forms together with the target system's P^!=), the paper's Q.  Algorithm
# MembershipLocus feeds Q into the differential Thomas decomposition itself,
# as the system's inequation (A(c), {Q}) -- so every component's own
# inequations already carry whatever of Q's non-degeneracy survives on it, and
# no PDE reduction needs Q multiplied in separately.  Q = 1 when the system
# declares no inequations, in which case feeding it into the decomposition
# would add a vacuous "1 != 0" and is skipped instead.
ANSATZ_INEQS = prob['ansatz_ineqs']          # the ansatz's own non-degeneracy forms
ALL_INEQS = list(ANSATZ_INEQS) + list(PDE_INEQS)
Q_INEQ = R(1)
for _q in ALL_INEQS:
    Q_INEQ = Q_INEQ * _q


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
        r"""Rewrite one matched name to bracket form, or leave it alone."""
        head, *idx = m.group(0).split('_')
        if idx and head and all(p in _COORD_SET for p in idx):
            return '%s[%s]' % (head, ','.join(idx))
        return m.group(0)
    return re.sub(r'[A-Za-z]\w*', repl, str(s))


for _P in PDES:
    print("PDE:", to_bracket(_P))
for _q in ANSATZ_INEQS:
    print("ansatz inequation:", to_bracket(_q), "!= 0")
for _q in PDE_INEQS:
    print("PDE inequation:", to_bracket(_q), "!= 0")
if not ALL_INEQS:
    print("inequations: none (Q = 1, membership locus unguarded)")
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
    r"""
    The sympy image of a BLAD differential polynomial.

    Walks the term list rather than the printed form, so no parsing is
    involved.  Each BLAD name is mapped by its kind: a bracketed derivative
    ``Psi[x,y]`` to the ``IndexedBase`` element ``_IB['Psi'][x, y]``, a bare jet
    to its ``IndexedBase``, a coordinate to its ``Symbol`` in ``_DERIV``, and a
    parameter to its ``Symbol`` in ``_PARAM``.

    INPUT:

    - ``e`` -- a BLAD differential polynomial (an element of ``R``)

    OUTPUT: the expanded sympy expression equal to ``e``

    A :class:`KeyError` is raised for a BLAD name of none of those kinds, which
    would mean the problem was built with an indeterminate this module never
    learned about.

    EXAMPLES::

        sage: _elt_to_sympy(R('Psi[x]'))         # not tested (needs R, _IB)
        Psi[x]
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
    r"""
    A sympy expression rendered in BLAD's input syntax.

    Two differences from :func:`sympy.sstr` have to be repaired: sympy prints an
    ``Indexed`` with spaces after the index commas (``Psi[x, y]``) where BLAD
    accepts none, and it writes powers as ``**`` where BLAD wants ``^``.

    INPUT:

    - ``expr`` -- a sympy expression

    OUTPUT: a string parseable by ``R(...)``

    EXAMPLES::

        sage: _sympy_to_blad_str(sympy.Symbol('a')**2)
        'a^2'
        sage: _sympy_to_blad_str(sympy.IndexedBase('Psi')[sympy.Symbol('x'),
        ....:                                             sympy.Symbol('y')])
        'Psi[x,y]'
    """
    s = sympy.sstr(expr)
    s = re.sub(r'\[([^\]]*)\]',
               lambda m: '[' + m.group(1).replace(' ', '') + ']', s)
    return s.replace('**', '^')


def has_jet(p):
    r"""
    Whether ``p`` involves any jet indeterminate.

    A name counts as a jet when its head -- the part before any ``[`` -- is one
    of the problem's jets, so both the bare ``Psi`` and the derivative
    ``Psi[x,y]`` answer ``True``.  Used to separate a cell's conditions on the
    constants alone from those that still carry an unknown function.

    INPUT:

    - ``p`` -- a BLAD differential polynomial

    OUTPUT: boolean

    EXAMPLES::

        sage: has_jet(R('Psi[x] - DPsi*v[x]'))    # not tested (needs R)
        True
        sage: has_jet(R('a0 + a1'))               # not tested (needs R)
        False
    """
    for _coeff, term in _blad.read_terms(p._h()):
        for nm, _deg in term:
            if nm.split('[', 1)[0] in _JET_HEADS:
                return True
    return False


def is_param_constancy(p):
    r"""
    Whether ``p`` is a constancy equation for a parameter.

    The decomposition is handed `\delta_i c = 0` for every constant `c`, and
    those equations come back inside the cells as terms in a *derivative of a
    parameter* -- a bracketed name whose head is in ``_PARAM``.  They carry no
    information about the cell (they hold identically) and so are filtered out
    wherever a cell's real content is read.

    INPUT:

    - ``p`` -- a BLAD differential polynomial

    OUTPUT: boolean

    EXAMPLES::

        sage: is_param_constancy(R('a0[x]'))      # not tested (needs R)
        True
    """
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

# sympy images of PolyRing's generators, in generator order.  _elt_to_sympy maps
# coords -> Symbol, bare jets -> IndexedBase, E/params -> Symbol, so this list is
# the exact mirror of PolyRing.gens().
def _blad_jet_to_sympy(blad_name):
    r"""
    The sympy image of a BLAD derivative-jet *name*.

    The string form of :func:`_elt_to_sympy`'s jet case, needed while building
    ``_SYMPY_GENS`` -- at that point there is a name but no BLAD element to
    read terms from.

    INPUT:

    - ``blad_name`` -- string; a bracketed derivative jet such as ``'Psi[x,y]'``

    OUTPUT: the corresponding ``sympy.Indexed``

    EXAMPLES::

        sage: _blad_jet_to_sympy('Psi[x,y]')     # not tested (needs _IB, _DERIV)
        Psi[x, y]
    """
    head, rest = blad_name.split('[', 1)
    return _IB[head][tuple(_DERIV[d] for d in rest.rstrip(']').split(','))]


_SYMPY_GENS = ([_DERIV[c] for c in COORDS]
               + [_IB[j] for j in prob['jets']]
               + [_blad_jet_to_sympy(b) for b, _m in EXTRA_JET_PAIRS]
               + [_PARAM[p] for p in PDE_PARAMS]
               + [_PARAM[p] for p in PARAMS])
assert len(_SYMPY_GENS) == PolyRing.ngens()


def _sympy_to_polyring(expr):
    r"""
    Convert a sympy expression to ``PolyRing`` without a string round-trip.

    .. WARNING::

        Unused.  The pipeline reaches ``PolyRing`` from BLAD directly, through
        :func:`_elt_to_polyring`, which never builds a sympy expression at all
        (see the comment above ``_GEN_IDX``); this is the older sympy detour,
        kept because it documents the failure it was written to avoid and
        because it is the only converter that accepts a sympy input.

    Sage has no sympy `\to` libsingular conversion, so ``PolyRing(expr)`` falls
    through to ``self(str(expr))`` and hands the result to :func:`eval`.
    CPython parses a sum as a left-nested tree of binary operations and its
    *compiler* recurses once per term, so that path dies with

    .. CODE-BLOCK:: text

        RecursionError: maximum recursion depth exceeded during compilation

    at about 2995 additive terms (python 3.11, recursion limit 1000) --
    regardless of how simple the polynomial is.  ``helium`` / ansatz 20.1
    crossed that ceiling.  Building from an exponent-to-coefficient dictionary
    never invokes the parser, so it has no such limit (measured: 4000 monomials
    in 0.04s).

    INPUT:

    - ``expr`` -- a sympy expression in the images of ``PolyRing``'s generators

    OUTPUT: the corresponding element of ``PolyRing``

    A :class:`TypeError` naming the offenders is raised when ``expr`` involves
    a symbol with no ``PolyRing`` generator -- typically a derivative jet that
    survived a reduction meant to eliminate every one of them.  The old string
    path masked such a leftover behind the ``RecursionError`` above, since
    compilation fails before name resolution.

    EXAMPLES::

        sage: _sympy_to_polyring(sympy.Integer(0))   # not tested (needs PolyRing)
        0
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
    r"""
    The ``PolyRing`` image of a BLAD differential polynomial.

    The direct BLAD `\to` libsingular path: each term's exponent vector is
    accumulated into a dictionary keyed by exponent tuple, which is handed to
    ``PolyRing`` in one call.  No sympy expression and no string are built on
    the way, so neither the quadratic ``out += mon`` of the sympy detour nor
    the parser ceiling described in :func:`_sympy_to_polyring` applies.

    INPUT:

    - ``e`` -- a BLAD differential polynomial

    OUTPUT: the corresponding element of ``PolyRing``

    A :class:`TypeError` is raised when a BLAD name has no ``PolyRing``
    generator.  Under the orderly ranking that means a derivative jet survived
    a reduction that should have eliminated it; under the elimination ranking
    the parametric derivative jets legitimately survive and are given
    generators of their own, so only an unexpectedly high order reaches this.

    EXAMPLES::

        sage: _elt_to_polyring(R('a0'))          # not tested (needs R, PolyRing)
        a0
    """
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


def build_system_of_equations(eqn, constants):
    r"""
    The coefficients of ``eqn`` regarded as a polynomial in the non-constants.

    Requiring a differential polynomial to vanish *identically* in the unknown
    functions is requiring each of its coefficients on the jet monomials to
    vanish, and those coefficients are polynomials in the constants alone.
    This performs that split: each term is divided into a power product in the
    indeterminates outside ``constants`` and a coefficient inside them, and the
    coefficients of equal power products are summed.

    INPUT:

    - ``eqn`` -- a multivariate polynomial

    - ``constants`` -- the generators to treat as constants; anything not in
      this collection is a non-constant the coefficients are taken over

    OUTPUT: a tuple of the distinct coefficients, deduplicated

    EXAMPLES::

        sage: R4.<x,a,b> = PolynomialRing(QQ)
        sage: sorted(build_system_of_equations(a*x + b*x + a, [a, b]), key=str)
        [a, a + b]

    A polynomial already free of non-constants is its own only coefficient::

        sage: build_system_of_equations(a + b, [a, b])
        (a + b,)
    """
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
        r"""
        Initialize the cell from its equations and inequations.

        INPUT:

        - ``eqs`` -- the cell's equations

        - ``ineqs`` -- the cell's inequations, each read as `\neq 0`

        EXAMPLES::

            sage: c = GenericCell([1, 2], [3])
            sage: c.eqs, c.ineqs
            ([1, 2], [3])
        """
        self.eqs = list(eqs)
        self.ineqs = list(ineqs)


def cell_eqs(ds):
    r"""
    The equations of a cell.

    One of the two accessors everything downstream uses, so that a cell can be
    either a differential system from the decomposition or a hand-built
    :class:`GenericCell`.

    INPUT:

    - ``ds`` -- a cell: a ``differentialthomas`` system or a
      :class:`GenericCell`

    OUTPUT: a list of the cell's equations

    EXAMPLES::

        sage: cell_eqs(GenericCell([1, 2], [3]))
        [1, 2]
    """
    if isinstance(ds, GenericCell):
        return list(ds.eqs)
    return list(dt.differential_system_equations(ds))


def cell_ineqs(ds):
    r"""
    The inequations of a cell, each read as `\neq 0`.

    The companion of :func:`cell_eqs`; see there.

    INPUT:

    - ``ds`` -- a cell: a ``differentialthomas`` system or a
      :class:`GenericCell`

    OUTPUT: a list of the cell's inequations

    EXAMPLES::

        sage: cell_ineqs(GenericCell([1, 2], [3]))
        [3]
    """
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

    INPUT:

    - ``e`` -- a BLAD differential polynomial

    OUTPUT:

    the discriminant as a BLAD element, or ``None`` when it is not a condition

    EXAMPLES::

        sage: _discriminant_in_leader(R('DDPsi^2 - a0'))   # not tested (needs R)
        4*a0
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
    r"""
    The :class:`GenericCell` of the equations ``eqs``.

    The cell a Thomas decomposition would call generic, built directly: the
    equations unchanged, and as inequations the conditions under which they are
    a regular triangular system -- each equation's initial and separant, plus
    its discriminant in its own leader where that is a condition at all (see
    :func:`_discriminant_in_leader`).  Duplicates are dropped by printed form,
    and a rational constant is skipped since `1 \neq 0` says nothing.

    INPUT:

    - ``eqs`` -- the ansatz equations, as BLAD differential polynomials

    OUTPUT: a :class:`GenericCell`

    EXAMPLES::

        sage: generic_cell([R('Psi[x] - DPsi*v[x]')])   # not tested (needs R)
        <GenericCell ...>
    """
    keep = [e for e in eqs if not e.is_zero()]
    ineqs, seen = [], set()

    def add(q):
        r"""Record ``q`` as an inequation, skipping zeros, constants, repeats."""
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


# --- differential-Thomas decomposition of the ansatz, WITH Q as its
# --- inequation (Algorithm MembershipLocus, line 1: (A(c), {Q})) ----------
if GENERIC_CELL:
    print("\n--generic-cell: skipping the differential Thomas decomposition; "
          "building the ansatz's generic cell (%d ansatz + %d constancy eqs) ..."
          % (len(ansatz0), len(pconst)), flush=True)
    cells_ds = [generic_cell(list(ansatz0) + list(pconst))]
    if ALL_INEQS:
        # Fold Q into this cell's inequations too, for consistency with the
        # decomposed path (both now prune on a component whose inequations
        # already carry Q).
        cells_ds[0].ineqs.append(Q_INEQ)
    print("-> 1 generic cell, %d inequations:" % len(cell_ineqs(cells_ds[0])),
          flush=True)
    for _q in cell_ineqs(cells_ds[0]):
        print("     ", to_bracket(_q), "!= 0", flush=True)
    print("   (the other cells of the decomposition -- where these vanish -- "
          "are NOT computed)\n" + "=" * 72, flush=True)
else:
    # Feed Q in as the system's inequation whenever there is one: {Q_INEQ}, a
    # single-element set containing the combined product, matching the
    # paper's own abbreviation (section "The Constant Loci") rather than the
    # raw list of individual inequations.  Q = 1 (ALL_INEQS empty) adds
    # nothing, so it is skipped rather than handed in as a vacuous "1 != 0".
    decompose_ineqs = [Q_INEQ] if ALL_INEQS else []
    print("\nComputing native DifferentialThomas decomposition of the ansatz "
          "together with Q (%d ansatz + %d constancy eqs%s) ..."
          % (len(ansatz0), len(pconst),
             ", Q from %d inequations" % len(ALL_INEQS) if ALL_INEQS
             else ""), flush=True)
    _t0 = time.time()
    cells_ds = dt.differential_thomas_decomposition(
        ansatz0 + pconst, decompose_ineqs, prob['rk'])
    _wall = time.time() - _t0
    print("-> %d cells in %.1fs\n" % (len(cells_ds), _wall) + "=" * 72, flush=True)


for i, ds in enumerate(cells_ds, 1):
    leaders = [e.leader() for e in cell_eqs(ds) if not e.is_zero()]
    jl = [L for L in leaders if L is not None]
    if len(jl) != len(set(jl)):
        print("  SOUNDNESS FAIL cell %d: repeated leader" % i, flush=True)

if CELLS_OUT:
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
    r"""
    Print the elapsed wall time since the script started.

    Called from each of the exit paths -- the ``--decompose-only`` early
    return, the normal end, and the error path -- so that a run always reports
    its cost even when it stops short.

    OUTPUT: ``None`` (the line is printed)

    EXAMPLES::

        sage: print_total_time()                 # not tested (needs _T_START)
        Total time: 4457.3s (1:14:17.29)
    """
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

    INPUT:

    - ``q`` -- a BLAD differential polynomial, read as the inequation
      `q \neq 0`

    OUTPUT:

    a tuple of elements of ``PolyRing``, read as "these do not all vanish", or
    ``None``

    EXAMPLES::

        sage: ineq_coeff_set(R('a0 + a1*v'))     # not tested (needs R, PolyRing)
        (a0, a1)
    """
    try:
        qp = _elt_to_polyring(q)
    except TypeError:
        return None
    if qp.is_zero():
        return None
    return tuple(build_system_of_equations(qp, PolyRing_constants))


def adapt_cell(ds):
    r"""
    A cell's conditions on the constants, sorted into the four kinds used later.

    The cell's equations and inequations arrive as differential polynomials
    mixing jets, coordinates and constants; downstream only their content in
    the constants matters, in three different ways.  Parameter-constancy
    equations are dropped throughout, holding identically (see
    :func:`is_param_constancy`).

    INPUT:

    - ``ds`` -- a cell: a ``differentialthomas`` system or a
      :class:`GenericCell`

    OUTPUT:

    a dictionary with four keys:

    - ``'param_eqs'`` -- the equations free of jets, in ``PolyRing``; the
      stratum of constant space the cell sits over

    - ``'param_ineqs'`` -- the inequations free of jets, in ``PolyRing``

    - ``'jet_ineqs'`` -- the inequations that do carry a jet, left as BLAD
      elements

    - ``'ineq_coeffs'`` -- every inequation's coefficient set (see
      :func:`ineq_coeff_set`), jet-carrying or not, which is the form the
      pruning and the piece conditions consume.  An inequation with no
      ``PolyRing`` image contributes nothing.

    EXAMPLES::

        sage: sorted(adapt_cell(GenericCell([], [])))
        ['ineq_coeffs', 'jet_ineqs', 'param_eqs', 'param_ineqs']
    """
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
    r"""
    The ansatz with a cell's parametric equations solved and substituted in.

    .. WARNING::

        Unused.  The pipeline reduces against the cell's own differential
        triangular equations instead, which avoids everything this function has
        to cope with: :func:`sympy.solve` introduces radicals, ``RootOf``
        objects and injected denominators, picks an arbitrary branch when the
        system has several solutions, and falls back to substituting zero when
        it cannot solve at all.  Kept as the record of the earlier approach.

    INPUT:

    - ``param_eqs`` -- the cell's equations in the constants

    OUTPUT:

    a pair ``(sub, spec)`` -- the substitution as a sympy dictionary, and the
    ansatz equations under it as BLAD elements, identically-zero ones dropped

    EXAMPLES::

        sage: specialize([])                     # not tested (needs the ansatz)
        ({}, [...])
    """
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
    r"""
    The reductor list for a specialized ansatz: its equations plus ``pconst``.

    .. WARNING::

        Unused, along with :func:`specialize` whose output it was written to
        consume.  The pipeline builds its reductors from the cell's own
        equations; appending ``pconst`` there proved to be redundant reductors
        at an extra cost per pass, the cell already carrying the constancy
        relations.

    INPUT:

    - ``spec`` -- the specialized ansatz equations

    OUTPUT: a list of BLAD differential polynomials to reduce against

    EXAMPLES::

        sage: reductors_for([])                  # not tested (needs pconst)
        [...]
    """
    return list(spec) + list(pconst)


def full_prem(p, reductors, max_passes=64):
    r"""
    Ritt's full reduction of ``p`` against ``reductors``, iterated to a fixpoint.

    ``R.differential_prem`` makes a single pass over the reductor list, and one
    pass is not enough: reducing a high derivative by one reductor can re-expose
    a lower derivative that an *earlier* reductor handles -- for instance
    ``Psi[R1,R1]`` reduces down to ``DPsi[R1]`` and then to ``n0*Psi[R1]``,
    whose leader belongs to a chain rule the pass has already gone by.  A single
    pass therefore leaves first-order jets unreduced whenever the system has
    first-derivative terms: invisible for ``hydrogen``'s pure Laplacian, wrong
    for ``helium``'s `2/R_i \, \partial/\partial R_i` terms.  Looping to a
    fixpoint gives the true normal form, and costs one extra no-op pass once the
    remainder is fully reduced.

    INPUT:

    - ``p`` -- the differential polynomial to reduce; coerced into ``R`` if it
      is not already an element

    - ``reductors`` -- the differential polynomials to reduce against

    - ``max_passes`` -- integer (default: 64); the iteration cap, a backstop
      against a non-terminating reduction rather than a limit ever reached in
      practice

    OUTPUT:

    a pair ``(r, h)`` -- the remainder, and the product of the initials the
    reduction multiplied through by, so that `h \cdot p \equiv r`

    EXAMPLES::

        sage: full_prem(R('Psi[x]'), [])         # not tested (needs R)
        (Psi[x], 1)
    """
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
    r"""
    A hashable identity for an ideal: its generators as a sorted string tuple.

    The buckets are keyed on this, so a prime surfacing from several cells is
    merged into one entry rather than reported once per cell.  Sorting makes the
    key independent of the order minAss returned the generators in.

    .. NOTE::

        This identifies ideals by *presentation*, not mathematically: two
        different generating sets of the same ideal get different keys.  Nothing
        here depends on that, since every key comes from minAss, which returns a
        canonical basis.

    INPUT:

    - ``P`` -- an ideal

    OUTPUT: a tuple of strings

    EXAMPLES::

        sage: R4.<x,y> = PolynomialRing(QQ)
        sage: prime_key(R4.ideal(y, x))
        ('x', 'y')
        sage: prime_key(R4.ideal(x, y)) == prime_key(R4.ideal(y, x))
        True
    """
    return tuple(sorted(str(g) for g in P.gens()))


def fmt_ideal(P):
    r"""
    An ideal as its generators alone -- ``Ideal (g1, ..., gk)``.

    Sage's own repr appends ``of Multivariate Polynomial Ring in x, y, z, ...
    over Rational Field``, the same 150-character tail on every ideal printed.
    Every ideal here lives in :data:`PolyRing`, so the tail carries no
    information and buries the generators, which are the part that differs.

    INPUT:

    - ``P`` -- an ideal

    OUTPUT: a string

    EXAMPLES::

        sage: R4.<x,y> = PolynomialRing(QQ)
        sage: fmt_ideal(R4.ideal(x, y))
        'Ideal (x, y)'
        sage: fmt_ideal(R4.ideal(R4.zero()))
        'Ideal (0)'
    """
    return "Ideal (%s)" % ", ".join(str(g) for g in P.gens())


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
    r"""
    A content hash of an ideal: its generators and its ring's variables.

    The cache key for the subprocess GTZ results.  The ring's variables are
    hashed alongside the generators because the same generator strings in a
    different variable order are a different computation; the tag the caller
    supplies only has to be unique within a run, and it is this hash that makes
    a cached result safe to reuse across runs.

    INPUT:

    - ``I`` -- an ideal

    OUTPUT: a 16-character hexadecimal string

    EXAMPLES::

        sage: R4.<x,y> = PolynomialRing(QQ)
        sage: len(_gtz_key(R4.ideal(x, y)))
        16
        sage: _gtz_key(R4.ideal(x, y)) == _gtz_key(R4.ideal(y, x))
        True
        sage: _gtz_key(R4.ideal(x, y)) == _gtz_key(R4.ideal(x))
        False
    """
    h = hashlib.sha1()
    h.update(','.join(sorted(str(g) for g in I.gens())).encode())
    h.update(('|' + ','.join(str(v) for v in I.ring().gens())).encode())
    return h.hexdigest()[:16]


def _gtz_script(I, primes_path):
    r"""
    A standalone Singular program computing ``minAssGTZ(I)`` under ``option(prot)``.

    The primes are written to their own file rather than to stdout, so the
    protocol stream and the result never have to be untangled from one another:
    the log stays tailable while the computation runs, and the result parses
    without having to skip past megabytes of protocol.

    INPUT:

    - ``I`` -- the ideal to decompose

    - ``primes_path`` -- string; where the Singular program should write its
      result, in the format :func:`_gtz_parse` reads

    OUTPUT: the Singular program, as a string

    EXAMPLES::

        sage: R4.<x,y> = PolynomialRing(QQ)
        sage: s = _gtz_script(R4.ideal(x*y), '/tmp/p')
        sage: 'ring gtzring = 0,(x,y),dp;' in s
        True
        sage: 'list gtzP = minAssGTZ(gtzI);' in s
        True

    The zero ideal still has to generate a syntactically valid program::

        sage: 'ideal gtzI =\n  0;' in _gtz_script(R4.ideal(R4.zero()), '/tmp/p')
        True
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
    r"""
    Singular's error lines from a protocol log, if any.

    Singular exits 0 even after a hard error, so the exit status says nothing
    and the log is the only place a failed run announces itself.  Its errors
    are the lines beginning ``?``.

    INPUT:

    - ``log_path`` -- string; path to the protocol log

    OUTPUT:

    a list of the stripped error lines; empty when there are none, and also
    when the log cannot be read at all, an unreadable log being the caller's
    problem to report rather than this function's

    EXAMPLES::

        sage: p = tmp_filename()
        sage: _ = open(p, 'w').write('std in ...\n? ideal not zero-dimensional\n')
        sage: _gtz_errors(p)
        ['? ideal not zero-dimensional']
        sage: _gtz_errors('/nonexistent/path')
        []
    """
    try:
        with open(log_path, errors='replace') as fh:
            return [ln.strip() for ln in fh if ln.lstrip().startswith('?')]
    except OSError:
        return []


def _gtz_parse(path, R):
    r"""
    Read back the primes file written by :func:`_gtz_script`.

    Singular may fold a long ``string(ideal)`` across several lines, so each
    prime's body is rejoined before being split on the generator commas -- safe
    because a polynomial never contains one.

    INPUT:

    - ``path`` -- string; the primes file to read

    - ``R`` -- the ring to build the ideals in

    OUTPUT: a list of ideals of ``R``, one per prime

    A :class:`RuntimeError` is raised when the closing marker is missing, which
    means a torn write from an interrupted run.  That is treated as no result at
    all rather than as a short one, so a half-written cache file causes a
    recomputation instead of silently dropping primes.

    EXAMPLES::

        sage: R4.<x,y> = PolynomialRing(QQ)
        sage: p = tmp_filename()
        sage: _ = open(p, 'w').write('\n'.join(
        ....:     ['===PRIMES_BEGIN===', '===PRIME===', 'x,y', '===PRIMES_END===']))
        sage: _gtz_parse(p, R4)
        [Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field]

    A truncated file is rejected rather than parsed::

        sage: _ = open(p, 'w').write('===PRIMES_BEGIN===\n===PRIME===\nx,y\n')
        sage: _gtz_parse(p, R4)
        Traceback (most recent call last):
        ...
        RuntimeError: ...no ===PRIMES_END=== marker (truncated?)
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
    r"""
    The minimal associated primes of ``I``, in-process or in a subprocess.

    By default this is Sage's ``I.minimal_associated_primes()``, which reaches
    the same ``primdec.lib`` ``minAssGTZ`` through libsingular.  Under
    ``--gtz-subprocess`` it instead writes a standalone Singular program and
    runs it, which buys three things the in-process call cannot give: a
    tailable ``option(prot)`` protocol stream (libsingular swallows it -- see
    the comment above), a wall-clock bound via ``--gtz-timeout``, and a result
    cached on disk so an interrupted run resumes instead of recomputing.

    INPUT:

    - ``I`` -- the ideal to decompose

    - ``tag`` -- string labelling this ideal's files in ``--gtz-dir``; it need
      only be unique within a run, the cache key proper being
      :func:`_gtz_key` of the ideal

    OUTPUT: a list of the minimal associated primes of ``I``

    A :class:`RuntimeError` is raised when Singular fails, when it exceeds
    ``--gtz-timeout``, or when its log carries error lines; an unusable cached
    result is reported and recomputed rather than raising.

    EXAMPLES::

        sage: minimal_associated_primes_gtz(I, 'cell1')   # not tested (needs the flags)
        [...]
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

# The cell inequations restricting each reported piece, keyed as the buckets are.
#
# A piece is not the variety V(P): the cell it came from carries inequations, and
# saturating by them at the GTZ step deleted only the primes lying WHOLLY inside
# their zero locus -- each surviving prime still carries the inequations as a
# residual condition cutting a proper closed subset out of it.  Writing
# U_c = union_j V(cs_j) for the locus cell c's inequations exclude (cs_j being
# the coefficient set of the j-th inequation, see ineq_coeff_set), the piece
# reported for prime P out of cells c1, c2, ... is
#
#     V(P) \ ( U_{c1} INTERSECT U_{c2} INTERSECT ... ) ,
#
# an intersection because the same prime surfacing from several cells contributes
# the UNION of their pieces.  piece_excl[key] records those U_c, one per
# contributing cell, as a pair (cell number, tuple of coefficient sets); an empty
# tuple of coefficient sets means that cell excluded nothing and the piece is all
# of V(P).  The cell number is carried so :func:`piece_conditions` can label each
# clause with the cell it came from -- the bucket's own cell list cannot supply
# it, since :func:`prune_enclosed` merges an absorbed piece's cells into the
# encloser's list without merging piece_excl (correctly: absorption happens only
# when the absorbed piece is contained, so the encloser's own conditions already
# describe the union).
piece_excl = {}


def _piece_unrestricted(key):
    r"""
    True when the piece reported for ``key`` is the whole variety `V(P)`.

    That happens when some contributing cell excluded nothing, since the removed
    locus is the intersection over contributing cells.  Keys absent from
    ``piece_excl`` -- as in the doctest of :func:`prune_enclosed`, which builds a
    bucket by hand -- count as unrestricted; every path that fills a bucket fills
    this alongside it.

    INPUT:

    - ``key`` -- a :func:`prime_key`, indexing :data:`piece_excl`

    OUTPUT: boolean

    EXAMPLES::

        sage: piece_excl.clear()
        sage: _piece_unrestricted(('x',))            # absent: counts as whole
        True
        sage: piece_excl[('x',)] = [(7, ())]         # cell 7 excluded nothing
        sage: _piece_unrestricted(('x',))
        True
        sage: piece_excl[('x',)] = [(7, ((1,),))]    # cell 7 excluded something
        sage: _piece_unrestricted(('x',))
        False
        sage: piece_excl.clear()
    """
    excl = piece_excl.get(key)
    if not excl:
        return True
    return any(len(u) == 0 for _num, u in excl)


def _excluded_components(key, R):
    r"""
    The locus excluded from the piece reported for ``key``, as a list of ideals.

    The piece is `V(P) \setminus D` with
    `D = \bigcap_c U_c` and `U_c = \bigcup_j V(I_{c,j})`, `I_{c,j}` being the
    ideal generated by the `j`-th inequation's coefficient set in cell `c` (see
    :data:`piece_excl`).  Distributing the intersection over the unions,

    .. MATH::

        D = \bigcup_{f} V\Bigl( \sum_c I_{c, f(c)} \Bigr),

    one selection `f` per choice of an inequation from each contributing cell.
    The returned list is those `\sum_c I_{c,f(c)}`, so `D` is the union of their
    varieties; an empty list means `D = \emptyset`.

    This is the ideal-valued form of the paper's single `h`: where the paper
    writes `D = V(h)` for `h` the product of a cell's constant inequations, a
    jet-carrying inequation collects to a coefficient *set* rather than one
    polynomial (see :func:`ineq_coeff_set`), so each `U_c` is a union of
    varieties of ideals instead of one hypersurface.

    INPUT:

    - ``key`` -- a :func:`prime_key`, indexing :data:`piece_excl`

    - ``R`` -- the ring to build the ideals in

    OUTPUT:

    a list of ideals whose varieties union to `D`; empty when `D = \emptyset`

    EXAMPLES::

        sage: R4.<x,y> = PolynomialRing(QQ)
        sage: piece_excl.clear()
        sage: _excluded_components(('nothing',), R4)         # no entry: D empty
        []
        sage: piece_excl[('k',)] = [(7, ((x,), (y,)))]        # one cell, two ineqs
        sage: _excluded_components(('k',), R4)
        [Ideal (x) of ..., Ideal (y) of ...]

    Two contributing cells give the pairwise sums, one per selection::

        sage: piece_excl[('k',)] = [(7, ((x,),)), (8, ((y,),))]
        sage: _excluded_components(('k',), R4)
        [Ideal (x, y) of ...]
        sage: piece_excl.clear()
    """
    excl = piece_excl.get(key)
    if not excl:
        return []
    per_cell = []
    for _num, u in excl:
        if len(u) == 0:
            return []                    # this cell excluded nothing -> D empty
        per_cell.append(list(u))
    comps = []
    for choice in itertools.product(*per_cell):
        gens = [g for cs in choice for g in cs]
        comps.append(R.ideal(gens) if gens else R.ideal(R.zero()))
    return comps


def _tidy_cond(g):
    r"""
    Normalize a condition polynomial for printing: squarefree, integral, positive.

    Only `V(g)` matters for a condition `g \neq 0`, so `g` may be replaced by its
    squarefree part -- the product of its distinct irreducible factors -- which
    turns the `-v_4^2 b_1^2` that a reduction naturally produces into `v_4 b_1`.
    Denominators are then cleared and the sign normalized, neither of which
    changes the locus either.

    INPUT:

    - ``g`` -- a polynomial, read as the condition `g \neq 0`

    OUTPUT: a polynomial with the same vanishing locus

    EXAMPLES::

        sage: R4.<u,v> = PolynomialRing(QQ)
        sage: _tidy_cond(-u^2*v^2)
        u*v
        sage: _tidy_cond(-1/2*u)
        u
        sage: _tidy_cond(R4(3))                  # constants pass through
        3
    """
    if g.is_constant():
        return g
    try:
        g = prod(f for f, _ in g.factor())
    except (ArithmeticError, NotImplementedError, TypeError):
        pass                             # factorization is an optimization only
    g = clear_denominators(g)
    try:
        if g.lc() < 0:
            g = -g
    except (AttributeError, TypeError):
        pass
    return g


def piece_conditions(P, key):
    r"""
    The residual inequations cutting the reported piece out of `V(P)`.

    The algorithm's output is a pair `(\mathfrak{p}, h)` -- an irreducible
    locally closed set, not the variety `V(\mathfrak{p})`.  This returns the
    `h` half, in the form the coefficient-set representation calls for.

    A point lies in the piece when it lies in `V(P)` and, **for at least one**
    contributing cell, **every** inequation of that cell is non-vanishing:
    `V(P) \setminus \bigcap_c U_c = V(P) \cap \bigcup_c U_c^{\,c}`.  So the
    return value is a disjunction over cells of a conjunction over that cell's
    inequations -- ``[(cell number, [coefficient set, ...]), ...]``, read as
    "some cell's clauses all hold", each coefficient set read as "these do not
    all vanish".

    Each coefficient is reduced modulo `P` first, which is what makes the
    conditions readable: on `V(P)` an inequation vanishes exactly when its
    normal forms do, so the reduction loses nothing and typically collapses a
    dozen coefficients to one or two.  Two cases fall out of the reduction and
    are handled rather than printed:

    - a coefficient set reducing to a non-zero constant is satisfied everywhere
      on `V(P)`, so its clause is dropped;
    - a coefficient set reducing to all zeros vanishes identically on `V(P)`,
      so that cell contributes nothing to the piece and its whole disjunct is
      dropped.  (:func:`prune` diverts such primes to the off-`N_Q` bucket, so
      this should not arise for a genuine prime; it is handled for safety.)

    ``None`` is returned when the piece is all of `V(P)` -- no contributing
    cell restricts it -- which is the case in which the variety and the piece
    coincide and nothing need be printed.

    INPUT:

    - ``P`` -- the prime whose piece is being described

    - ``key`` -- its :func:`prime_key`, indexing :data:`piece_excl`

    OUTPUT:

    a list of pairs ``(cell number, [coefficient set, ...])``, or ``None``

    EXAMPLES::

        sage: R4.<u,v> = PolynomialRing(QQ)
        sage: P = R4.ideal(u)
        sage: piece_excl.clear()
        sage: piece_excl[prime_key(P)] = [(7, ((v,), (-v^2,)))]
        sage: piece_conditions(P, prime_key(P))     # the two clauses coincide
        [(7, [(v,)])]

    A condition vanishing on all of `V(P)` empties that cell's disjunct, and a
    condition holding everywhere on `V(P)` is dropped as vacuous::

        sage: piece_excl[prime_key(P)] = [(7, ((u,),))]
        sage: piece_conditions(P, prime_key(P)) is None
        True
        sage: piece_excl[prime_key(P)] = [(7, ((R4(2),),))]
        sage: piece_conditions(P, prime_key(P)) is None
        True
        sage: piece_excl.clear()
    """
    excl = piece_excl.get(key)
    if not excl:
        return None
    clauses = []
    for num, u in excl:
        if len(u) == 0:
            return None                  # this cell excluded nothing
        cell_clause = []
        for cs in u:
            red = [_tidy_cond(r) for r in (P.reduce(g) for g in cs)
                   if not r.is_zero()]
            if not red:
                cell_clause = None       # vanishes identically on V(P)
                break
            if any(r.is_constant() for r in red):
                continue                 # non-vanishing everywhere on V(P)
            clause = tuple(sorted(set(red), key=str))
            if clause not in cell_clause:
                cell_clause.append(clause)
        if cell_clause is None:
            continue
        if not cell_clause:
            return None                  # every clause vacuous: piece is V(P)
        clauses.append((num, cell_clause))
    return clauses or None


def _fmt_conditions(clauses):
    r"""
    Render :func:`piece_conditions` output as one line of plain text.

    A single-polynomial coefficient set prints as ``g != 0``; a larger one as
    ``(g1, g2) not all 0``, that being what a jet-carrying inequation's
    coefficient set means.  Disjuncts are separated by ``OR``, and the cell
    number is shown when more than one cell contributes.

    INPUT:

    - ``clauses`` -- the output of :func:`piece_conditions`

    OUTPUT: a string

    EXAMPLES::

        sage: R4.<u,v> = PolynomialRing(QQ)
        sage: _fmt_conditions([(7, [(u,), (u, v)])])
        'u != 0 and (u, v) not all 0'

    With more than one disjunct the cells are named::

        sage: _fmt_conditions([(7, [(u,)]), (8, [(v,)])])
        '[cell 7] u != 0   OR   [cell 8] v != 0'
    """
    def one(cs):
        r"""Render one coefficient set as a non-vanishing condition."""
        return ("%s != 0" % cs[0] if len(cs) == 1
                else "(%s) not all 0" % ", ".join(map(str, cs)))
    parts = []
    for num, cell_clause in clauses:
        body = " and ".join(one(cs) for cs in cell_clause)
        parts.append(body if len(clauses) == 1 else "[cell %d] %s" % (num, body))
    return "   OR   ".join(parts)


def _piece_contained(Pi, ki, Pj, kj):
    r"""
    Test `X_i \subseteq X_j` for the pieces `X = V(P) \setminus D` of two primes.

    This is the paper's containment test, equation (containment-test): a piece
    `(\mathfrak{p}_1, h_1)` is contained in `(\mathfrak{p}_2, h_2)` exactly when

    .. MATH::

        \mathfrak{p}_2 \subseteq \mathfrak{p}_1
        \quad\text{and}\quad
        h_1 \in \sqrt{\mathfrak{p}_1 + \langle h_2 \rangle},

    the first cheap and necessary (it says `V(P_i) \subseteq V(P_j)`), the second
    saying the smaller piece avoids the locus `D_j` the larger one excludes.  In
    the ideal-valued form :func:`_excluded_components` returns, the second
    condition is `V(P_i) \cap D_j \subseteq D_i`, which we decide exactly: for
    each component `V(J_f)` of `D_j`, split `V(P_i) \cap V(J_f)` into its
    irreducible components and require each to sit inside some component of
    `D_i`.  An irreducible set lies in a finite union of closed sets only if it
    lies in one of them, so the componentwise test is not merely sufficient.

    `D_j = \emptyset` short-circuits to ``True`` -- the sound special case the
    older :func:`_piece_unrestricted` test was restricted to.

    INPUT:

    - ``Pi`` -- the prime of the candidate contained piece

    - ``ki`` -- its :func:`prime_key`

    - ``Pj`` -- the prime of the candidate containing piece

    - ``kj`` -- its :func:`prime_key`

    OUTPUT: boolean

    EXAMPLES::

        sage: R4.<u,v> = PolynomialRing(QQ)
        sage: Pi, Pj = R4.ideal(u, v), R4.ideal(u)
        sage: piece_excl.clear()
        sage: _piece_contained(Pi, ('i',), Pj, ('j',))    # nothing excluded
        True
        sage: _piece_contained(Pj, ('j',), Pi, ('i',))    # the other way round
        False

    When the encloser excludes a locus the smaller piece lies in, containment
    fails -- the case the old test could not see::

        sage: piece_excl[('j',)] = [(8, ((v,),))]         # V(Pi) is inside V(v)
        sage: _piece_contained(Pi, ('i',), Pj, ('j',))
        False
        sage: piece_excl.clear()
    """
    if not (Pj <= Pi):                   # V(Pi) subset V(Pj); necessary
        return False
    R = Pi.ring()
    Dj = _excluded_components(kj, R)
    if not Dj:
        return True                      # nothing excluded from the encloser
    Di = _excluded_components(ki, R)
    for Jf in Dj:
        I = Pi + Jf
        if I.is_one():
            continue                     # V(Pi) misses this component of D_j
        if not Di:
            return False                 # nowhere for the intersection to hide
        for comp in I.minimal_associated_primes():
            if not any(Kg <= comp for Kg in Di):
                return False
    return True

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
        # Reduce EACH PDE of the system against the same cell, unmultiplied.
        # Q's non-degeneracy was already imposed when this cell was built --
        # Q was fed into the differential Thomas decomposition as the system's
        # inequation (see the decomposition call above) -- so the cell's own
        # equations already reflect it, and a plain reduction is the
        # membership test Algorithm MembershipLocus asks for (line 5).  The
        # constant-coefficient equations of all the remainders are combined
        # into ONE system, solved once per cell.
        rem_elts = [full_prem(P_, reductors)[0] for P_ in PDES]
        t_prem = time.time() - _t
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
                                       primes=primes)

    sc = strata_cache[cache_key]
    # A prime is dropped when some inequation of the cell vanishes identically
    # on the whole of it -- every coefficient of the inequation lying in the
    # prime.  See ineq_coeff_set: on inequations in the constants alone this is
    # the `q in P` test it replaces.  This is Algorithm MembershipLocus's g_i
    # saturation (paper, lines 11-12): cp['ineq_coeffs'] already covers every
    # inequation of the cell, Q-derived or not, since Q was fed into the
    # decomposition that produced this cell -- there is no separate Q test.
    survivors = []
    for P in sc['primes']:
        if P.is_one():
            continue
        if any(all(g in P for g in cs) for cs in cp['ineq_coeffs']):
            continue
        survivors.append(P)

    print("\n--- cell %d: zero {%s}; ansatz %d eqs; %d param-ineqs, %d jet-ineqs ---"
          % (num, ', '.join(Zkey) or '(none, generic)', sc['spec_len'],
             len(cp['param_ineqs']), len(cp['jet_ineqs'])), flush=True)
    if VERBOSE_REM:
        for _i, _r in enumerate(sc['rems'], 1):
            print("  remainder[pde %d]:" % _i, to_bracket(_r), flush=True)
    if all(r_.is_zero() for r_ in sc['rems']):
        print("  PDE reduces to 0: the whole stratum solves the PDE (nontrivially)", flush=True)
    if not survivors:
        print("  surviving solution varieties: NONE (all pruned / empty)", flush=True)
    for P in survivors:
        print("   V:", fmt_ideal(P), "  [GENUINE]", flush=True)
        union_primes.setdefault(prime_key(P), (P, []))[1].append(num)
        piece_excl.setdefault(prime_key(P), []).append(
            (num, tuple(cp['ineq_coeffs'])))


def prune_enclosed(d):
    r"""
    Keep only the maximal (enclosing) varieties in a bucket of solution primes.

    A piece is dropped when another piece of the bucket contains it.  The
    pieces are NOT the varieties `V(P)`: each is `V(P)` less the locus its cell's
    inequations exclude (see :data:`piece_excl`), so it is locally closed, and
    reverse ideal containment `P_j \subseteq P_i` is NOT by itself the test.  It
    gives `V(P_i) \subseteq V(P_j)`, but a point of the smaller piece may lie in
    the locus the LARGER one excludes, and then it is not in the larger piece at
    all.  Pruning on the prime alone therefore discards real solution families.

    The exact condition for `X_i \subseteq X_j` is `P_j \subseteq P_i` together
    with `X_i \cap D_j = \emptyset`, where `D_j` is the locus excluded from
    `X_j`; this is the paper's equation (containment-test).  It is decided in
    full by :func:`_piece_contained`.  Earlier versions tested only the sound
    special case `D_j = \emptyset` (:func:`_piece_unrestricted`, kept for the
    short-circuit it still provides), which pruned nothing whenever the
    enclosing piece's cell carried a constant-space inequation -- the common
    case, so nested pieces routinely survived into the output.

    Pruning is *pairwise*, as in the paper: a piece contained in the union of
    two others is not detected.  Making the presentation irredundant in that
    stronger sense is what SMPD does, and this branch declines to do it.

    Equal primes are already merged by :func:`prime_key` before we get here, so
    the mutual-containment tie of the older version cannot arise; the guard
    against it is kept for the hand-built buckets of the doctest.  Each dropped
    piece's source-cell list is folded into the entry that encloses it, so no
    cell provenance is lost.

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
            # The full test of equation (containment-test): ideal containment
            # `Pj <= Pi` plus the smaller piece avoiding the locus the larger
            # one excludes.  See :func:`_piece_contained`.
            if _piece_contained(Pi, entries[i][0], Pj, entries[j][0]):
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
        print("  V:", fmt_ideal(P), "  (from cells:",
              ", ".join(map(str, sorted(set(cells_for)))) + ")")
        # The piece is V(P) less what its cells exclude; printing the variety
        # alone overstates the answer on a proper closed subset of it, and
        # invites the reader to compare varieties where the objects are locally
        # closed (two pieces can be disjoint while one variety sits inside the
        # other).  So print the h half of the (p, h) pair alongside it.
        conds = piece_conditions(P, key)
        if conds:
            print("       where", _fmt_conditions(conds))
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

    EXAMPLES::

        sage: R4.<u,v> = PolynomialRing(QQ)
        sage: piece_excl.clear()
        sage: latex_union({('u',): (R4.ideal(u), [1])}, label='ex')
        \begin{subequations}
        \label{ex}
        \begin{align}
        & \left(u\right)\label{ex:1}
        \end{align}
        \end{subequations}

    A restricted piece carries its conditions alongside the generators::

        sage: P = R4.ideal(u)
        sage: piece_excl[prime_key(P)] = [(7, ((v,),))]
        sage: latex_union({prime_key(P): (P, [7])}, label='ex')  # abbreviated
        \begin{subequations}
        \label{ex}
        \begin{align}
        & \left(u\right) \quad\text{with}\quad v \neq 0\label{ex:1}
        \end{align}
        \end{subequations}
        sage: piece_excl.clear()
    """
    items = sorted(d.items(), key=lambda kv: str(kv[0]))
    rows = []
    for i, (key, (P, _cells)) in enumerate(items, 1):
        gens = ", ".join(latex(clear_denominators(g)) for g in P.gens())
        # The (p, h) pair, not p alone -- see the comment in dump_union.
        conds = piece_conditions(P, key)
        tail = ""
        if conds:
            def _one(cs):
                r"""Render one coefficient set as a LaTeX non-vanishing condition."""
                if len(cs) == 1:
                    return r"%s \neq 0" % latex(clear_denominators(cs[0]))
                return (r"\left(%s\right) \neq 0"
                        % ", ".join(latex(clear_denominators(g)) for g in cs))
            disj = r" \;\text{ or }\; ".join(
                r" ,\; ".join(_one(cs) for cs in cell_clause)
                for _num, cell_clause in conds)
            tail = r" \quad\text{with}\quad %s" % disj
        rows.append(r"& \left(%s\right)%s\label{%s:%d}" % (gens, tail, label, i))
    print(r"\begin{subequations}")
    print(r"\label{%s}" % label)
    print(r"\begin{align}")
    print(" \\\\\n".join(rows))
    print(r"\end{align}")
    print(r"\end{subequations}")


print("\n" + "=" * 72)
genuine_primes = dump_union(
          "GENUINE solution varieties over all cells (the real union)",
          union_primes, prune=True)

print("\n" + "=" * 72)
if genuine_primes:
    print("VERDICT: %d GENUINE solution variety(ies) found for %s / ansatz %s."
          % (len(genuine_primes), PDE_NAME, ANSATZ))
else:
    print("VERDICT: NO genuine solution for %s / ansatz %s "
          "(every stratum was pruned by some inequation of its cell, "
          "Q included)."
          % (PDE_NAME, ANSATZ))

if LATEX_OUT and genuine_primes:
    print("\n" + "-" * 72)
    print("GENUINE solution varieties, LaTeX (paper) form:\n")
    latex_union(genuine_primes)

print("\n" + "=" * 72)
print_total_time()
