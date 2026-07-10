# -*- mode: python -*-
#
# joca-rg-combined.sage --- stock Rosenfeld-Gröbner on the COMBINED system
# A u {P}: the hydrogen ansatz with the Schrödinger PDE adjoined.
#
# The other RG scripts in this directory run RG on the ANSATZ ALONE:
#
#   joca.sage         reduce P mod the raw ansatz, project, decompose  -> 5 primes
#   joca-rg.sage      RG(ansatz), then reduce P mod the returned chain -> 4 primes
#                     (one spurious; see rg-saturation-and-the-bad-locus.tex)
#   rg_basefield.py   RG(ansatz) with the constants in Q(c); 1 component
#
# This script instead hands RG the combined system.  Per
# membership-vs-variety-partial-strata.tex the two computations answer
# different questions:
#
#   membership locus   V = { c : P in [A(c)] }         (joca.sage)
#   joint variety      pi_c(Sol(A) n Sol(P)) = V u {partial-solution strata}
#
# so the combined decomposition is expected to return, besides the trivial cell
# Psi = 0, cells whose constant-relations are V together with extra strata on
# which only a proper sub-family of the ansatz solves the PDE.  This script
# computes both and classifies every nontrivial cell as membership or partial.
#
# ---------------------------------------------------------------------------
# WHERE DO THE CONSTANTS LIVE?
#
# The eleven constants E,v1..v4,a0,a1,b0,b1,c0,c1 can sit in two places, and the
# choice decides whether RG can see the strata at all.
#
#   default (no flag)           constants are RING VARIABLES in the lowest
#       block.  Initials that involve them are zero-divisors, so RG splits on
#       their vanishing and the parameter strata appear as constant-only
#       equations of the returned components.  This is the configuration that
#       can answer the question -- and the one that does NOT terminate on the
#       hydrogen system (see the README's performance findings: RG spins in
#       BLAD's bad_remainder_irreducible_factorwise / baz_Yun factorization,
#       under both elimination and orderly rankings).  Use --timeout / --memout.
#
#   --basefield                 constants are moved into the coefficient field
#       Q(E,v1,...,c1) via BaseFieldExtension.  Every nonzero constant
#       polynomial is a UNIT, so RG never splits on a parameter locus, and by
#       construction it can only return GENERIC cells: the parameter strata are
#       invisible.  On the ANSATZ this terminates (~27 s, one component) and
#       reproduces joca-rg.sage / rg_basefield.py.
#
# Both are worth running, and --toy exhibits the contrast in a tenth of a second
# on the note's own example (A: Psi'' - c*Psi, P: Psi' - Psi):
#
#   --toy               2 components: Psi = 0, and {c - 1, Psi[x] - Psi}
#                       -- the c = 1 partial-solution stratum, as the note says
#   --toy --basefield   1 component:  Psi = 0
#                       -- p_alpha(c) are units, p_alpha(c) = 0 is inconsistent,
#                          the stratum is GONE.  Not a bug: it is the precise
#                          sense in which the base field cannot resolve V.
#
# MEASURED, so you are not surprised (samsung, 2026-07-09):
#
#   --toy                              0.1 s   2 components
#   --toy --basefield                  0.1 s   1 component
#   --ansatz-only --basefield         27   s   1 component
#   --timeout 20                      abort    (ring constants; the known cliff)
#   --basefield  (combined)           >420 s   did NOT finish -- see below
#
# So on hydrogen the base field does NOT rescue the combined system the way it
# rescues the ansatz: RG still fails to terminate in seven minutes.  The
# collapse-to-Psi=0 argument above is proved on the toy and is what one expects
# on hydrogen; it has NOT been observed there.  Do not quote it as a hydrogen
# result without a finishing run.
#
# ---------------------------------------------------------------------------
# USAGE
#
#   sage joca-rg-combined.sage [options]
#
#     --basefield         put the constants in Q(constants) (default: ring vars)
#     --toy               run the D^2 - c example instead of hydrogen
#     --ansatz-only       decompose A alone, not A u {P} (regression vs joca-rg)
#     --timeout SECS      RosenfeldGroebner wall limit (0 = none; default 0)
#     --memout MB         RosenfeldGroebner memory limit (0 = none; default 0)
#     --no-membership     skip the joca.sage reference computation
#     --print-inequations print each component's initials and separants
#     --rg-verbose        trace RG's splitting (patched DifferentialAlgebra)
#     --rg-dot            dump RG's splitting tree as graphviz dot
#
# Recommended first runs:
#
#   sage joca-rg-combined.sage --toy                       # strata appear
#   sage joca-rg-combined.sage --toy --basefield           # strata vanish
#   sage joca-rg-combined.sage --ansatz-only --basefield   # 27 s regression
#   sage joca-rg-combined.sage --timeout 3600              # hydrogen, the real
#                                                          # question; expect the
#                                                          # cliff, not an answer
#
# Whatever RG does, the membership locus V is printed first, from the projection
# route, which needs no RG at all: on hydrogen it reproduces joca.sage's five
# minimal primes in about a second.  A timeout therefore never costs you V.
#
# Author: Claude (Fable 5), for Brent Baccala
# Date: July 9, 2026
#
# Tested with the conda `sage` env (SageMath 10.7) + DifferentialAlgebra 5.3.

import sys
import time

import sympy

try:
    import DifferentialAlgebra
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(ex.msg + "\nInstall it with '%pip install DifferentialAlgebra'")


def _flag(name):
    return name in sys.argv


def _val(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


use_basefield     = _flag('--basefield')
toy               = _flag('--toy')
ansatz_only       = _flag('--ansatz-only')
skip_membership   = _flag('--no-membership')
print_inequations = _flag('--print-inequations')
rg_verbose        = _flag('--rg-verbose')
rg_dot            = _flag('--rg-dot')
rg_timeout        = _val('--timeout', 0)
rg_memout         = _val('--memout', 0)

T0 = time.time()


def stamp(msg):
    print(f"[{time.time() - T0:7.2f}s] {msg}", flush=True)


# ---------------------------------------------------------------------------
# The system: hydrogen (default) or the note's two-line example (--toy).
#
# Coefficients must be sympy Rationals, never Python floats: BLAD's parser
# rejects `0.5`.  (joca.sage writes -int(1)/int(2), which is a float and works
# only because that script never hands the PDE to RosenfeldGroebner -- it goes
# to differential_prem instead.)
# ---------------------------------------------------------------------------

if toy:
    x = sympy.var('x')
    derivations = [x]
    c = sympy.var('c')
    constants = [c]

    Psi = DifferentialAlgebra.indexedbase('Psi')
    dependents = [Psi]

    def build_system():
        A = [Psi[x, x] - c * Psi]                     # ansatz: two-dimensional
        P = Psi[x] - Psi                              # PDE: first order
        return A, P

    system_name = "toy (membership-vs-variety-partial-strata.tex, Example 1)"

else:
    x, y, z = sympy.var('x,y,z')
    derivations = [x, y, z]

    E = sympy.var('E')
    v1, v2, v3, v4 = sympy.var('v1,v2,v3,v4')
    a0, a1, b0, b1, c0, c1 = sympy.var('a0,a1,b0,b1,c0,c1')
    constants = [E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]

    Psi, DPsi, DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
    v = DifferentialAlgebra.indexedbase('v')
    r = DifferentialAlgebra.indexedbase('r')
    dependents = [DDPsi, DPsi, Psi, v, r]

    def build_system():
        A = [Psi[x] - DPsi * v[x],
             Psi[y] - DPsi * v[y],
             Psi[z] - DPsi * v[z],
             DPsi[x] - DDPsi * v[x],
             DPsi[y] - DDPsi * v[y],
             DPsi[z] - DDPsi * v[z],
             (a0 + a1 * v) * DDPsi + (b0 + b1 * v) * DPsi + (c0 + c1 * v) * Psi,
             v - (v1 * x + v2 * y + v3 * z + v4 * r),
             r**2 - x**2 - y**2 - z**2]
        P = (-sympy.Rational(1, 2) * (Psi[x, x] + Psi[y, y] + Psi[z, z]) * r
             - Psi - E * r * Psi)
        return A, P

    system_name = "hydrogen (JOCA ansatz + Schrödinger PDE, 11 constants)"


DiffRing = DifferentialAlgebra.DifferentialRing(derivations=derivations,
                                                blocks=[dependents, constants],
                                                parameters=constants,
                                                notation='jet')

ansatz, PDE = build_system()
ansatz = list(map(sympy.expand, ansatz))
PDE = sympy.expand(PDE)

print("=" * 78)
print("System:   ", system_name)
print("Constants:", "base field Q(%s)" % ", ".join(map(str, constants))
      if use_basefield else "ring variables, lowest block")
print("Decompose:", "ansatz alone" if ansatz_only else "ansatz u {PDE}  (COMBINED)")
print("=" * 78)
print("\nPDE:", PDE)
print("\nAnsatz:", *ansatz, sep='\n')

# ---------------------------------------------------------------------------
# Sympy -> Sage conversion.
#
# CAREFUL.  joca.sage builds its polynomial ring from
# `DiffRing.indets(selection='all')`, which returns only the BASE
# indeterminates (x,y,z,Psi,DPsi,DDPsi,v,r,constants) -- no jet variables.
# That is sound there only because Ritt reduction against the hydrogen ansatz
# eliminates every derivative, so the remainder is free of jets.  It is NOT
# sound in general: whenever the PDE's leader ranks BELOW the ansatz's, the
# remainder keeps a jet variable (the --toy case, where the remainder is
# Psi[x] - Psi), and Sage's conversion then silently drops it -- `Psi[x] - Psi`
# converts to `-Psi`, corrupting the projection.
#
# So we rename jet variables to legal identifiers (Psi[x,x] -> Psi_x_x) and
# build the ring from the atoms actually present.
# ---------------------------------------------------------------------------

# The ring of the constants.  A one-constant system (--toy) would give a
# UNIVARIATE ring, whose ideals do not implement minimal_associated_primes; a
# dummy generator keeps us in the multivariate implementation.  It never occurs
# in any polynomial we build, so it changes no ideal and no prime.
_const_names = [str(k) for k in constants]
ConstRing = PolynomialRing(QQ, names=_const_names + (['dummy0'] if len(_const_names) == 1 else []))
_const_gens = dict(zip(_const_names, ConstRing.gens()))


def to_const(expr):
    """Convert a sympy expression or a Sage polynomial in the constants into
    ConstRing.  Goes through the string form because sympy objects do not
    coerce into a Sage polynomial ring directly."""
    return ConstRing(sage_eval(str(expr), locals=_const_gens))


def _dejet(expr):
    """Replace every jet atom Psi[x,x] by a plain symbol Psi_x_x."""
    expr = sympy.expand(expr)
    return expr.xreplace({a: sympy.Symbol("%s_%s" % (a.base, "_".join(map(str, a.indices))))
                          for a in expr.atoms(sympy.Indexed)})


def _numer(expr):
    """Clear denominators; base-field mode yields rational coefficients."""
    return sympy.fraction(sympy.together(sympy.expand(expr)))[0]


def poly_ring_for(expr):
    """A Sage polynomial ring with a generator for every symbol of `expr`,
    plus every constant (which may not occur in `expr`)."""
    names = {str(s) for s in _dejet(_numer(expr)).free_symbols}
    names |= {str(k) for k in constants}
    return PolynomialRing(QQ, names=sorted(names))


def is_constant_only(expr):
    """True iff `expr` involves none of the differential indeterminates.

    A jet variable Psi[x] contributes {Psi, x} to free_symbols and a bare
    dependent Psi contributes {Psi}; the derivations contribute {x,y,z}.  So
    testing free_symbols against the constants is exact.
    """
    return bool(expr.free_symbols) and expr.free_symbols <= set(constants)


def build_system_of_equations(eqn, consts):
    """joca.sage's projection: group the remainder's terms by their
    non-constant power product; the coefficients are the projection
    equations."""
    ring = eqn.parent()
    system_of_like_terms = dict()
    non_constant_sub = tuple(1 if ring.gen(n) in consts else ring.gen(n)
                             for n in range(ring.ngens()))
    for coeff, monomial in eqn:
        non_constant_part = monomial(non_constant_sub)
        constant_part = coeff * monomial // non_constant_part
        system_of_like_terms[non_constant_part] = \
            system_of_like_terms.get(non_constant_part, ring.zero()) + constant_part
    return tuple(set(system_of_like_terms.values()))


def minimal_primes(gens):
    """Minimal associated primes of the ideal generated by `gens` in ConstRing,
    handling the two degenerate ideals Singular's minAssGTZ refuses."""
    gens = [to_const(g) for g in gens]
    gens = [g for g in gens if g != 0]
    if not gens:                                  # zero ideal: all of C^n
        return [ConstRing.ideal(ConstRing.zero())]
    I = ConstRing.ideal(gens)
    if I.is_one():                                # unit ideal: the empty locus
        return []
    primes = I.minimal_associated_primes()
    primes.sort(key=lambda p: str(p))
    return primes


# ---------------------------------------------------------------------------
# Reference: the membership locus V, by joca.sage's route (Ritt-reduce the PDE
# modulo the RAW ansatz, project onto the constants, take minimal primes).
# Cheap (~5 s on hydrogen) and independent of everything RG does below, so it
# runs first: if the RG call later times out we still have V in the log.
# ---------------------------------------------------------------------------

membership_primes = None
membership_eqns = None

if not skip_membership and not ansatz_only:
    stamp("membership locus V: differential_prem(PDE, ansatz)")
    _, remainder = DiffRing.differential_prem(PDE, ansatz)
    print("\nReduced PDE (remainder of Ritt reduction mod the ansatz):", remainder)

    PolyRing = poly_ring_for(remainder)
    PolyRing_constants = list(map(PolyRing, constants))
    membership_eqns = build_system_of_equations(PolyRing(_dejet(_numer(remainder))),
                                                PolyRing_constants)
    print("\nProjection equations (coefficients of the reduced PDE):",
          *membership_eqns, sep='\n')

    membership_primes = minimal_primes(membership_eqns)
    stamp(f"V has {len(membership_primes)} minimal prime(s)")
    if membership_primes:
        print("\nMembership locus V = V(p_1,...,p_N), minimal associated primes:",
              *membership_primes, sep='\n')
    else:
        print("\nMembership locus V = EMPTY: some projection coefficient is a")
        print("nonzero constant, so no c makes the whole ansatz family solve the PDE.")

# ---------------------------------------------------------------------------
# The Rosenfeld-Gröbner call.
# ---------------------------------------------------------------------------

system = ansatz if ansatz_only else ansatz + [PDE]

if use_basefield:
    F = DifferentialAlgebra.BaseFieldExtension(generators=constants, ring=DiffRing)
    stamp("base field extension built over the constants")
else:
    F = None
    print("\nNOTE: the constants are ring variables.  This is the configuration"
          "\n      in which RG can split on parameter loci -- and the one the"
          "\n      README reports does not terminate on the hydrogen ansatz."
          "\n      Consider --timeout / --memout.")

stamp(f"RosenfeldGroebner({'ansatz' if ansatz_only else 'ansatz + [PDE]'}"
      f"{', basefield=F' if use_basefield else ''}) START")

try:
    components = DiffRing.RosenfeldGroebner(system, basefield=F,
                                            timeout=rg_timeout, memout=rg_memout,
                                            verbose=rg_verbose, dot=rg_dot)
except RuntimeError as ex:
    # BLAD reports both limits as RuntimeError ("out of time error" /
    # "out of memory error").  Report and exit nonzero rather than traceback:
    # a timeout here is an expected, documented outcome of the ring-constants
    # configuration, not a crash.
    stamp(f"RosenfeldGroebner ABORTED: {ex}")
    print("\nThe combined decomposition did not finish under the given limits.")
    if not use_basefield:
        print("This is the documented ring-constants cliff (BLAD factorwise")
        print("reduction).  The membership locus V above is unaffected -- it")
        print("comes from the projection route, which needs no RG.")
    sys.exit(1)

stamp(f"RosenfeldGroebner DONE: {len(components)} component(s)")

# ---------------------------------------------------------------------------
# Report each component, and classify it.
#
#   trivial     the component forces Psi = 0 -- no solution, carries no stratum
#   membership  its constant-relations imply every projection equation, so the
#               whole ansatz family solves the PDE there:  stratum <= V
#   partial     a nontrivial cell that is NOT inside V: a proper sub-family of
#               the ansatz solves the PDE (Definition 1 of the note)
#
# The membership test is ideal-theoretic and exact: the stratum V(J) of a
# component lies in V = V(p_1,...,p_N) iff every p_alpha vanishes on V(J), i.e.
# iff p_alpha lies in the radical of J.  (Saturation caveat: J is read off the
# component's constant-only equations, so this classifies the cell's Zariski
# closure; the component's inequations carve out the open part.)
# ---------------------------------------------------------------------------

def psi_is_zero(chain):
    """The component forces the trivial solution Psi = 0."""
    return any(sympy.expand(e) == Psi for e in chain.equations())


summary = []

for i, C in enumerate(components):
    eqns = C.equations()
    print("\n" + "-" * 78)
    print(f"Component {i}:  {len(eqns)} equations")
    print("-" * 78)
    for e in eqns:
        print("   ", e)

    if print_inequations:
        seen, ineqs = set(), []
        for h in list(C.initial()) + list(C.separant()):
            h = sympy.expand(C.normal_form(h))
            if h.is_number:
                continue
            h = h.as_content_primitive()[1]
            if h not in seen:
                seen.add(h)
                ineqs.append(h)
        print("\n  Inequations (initials and separants), must be nonzero:")
        for h in ineqs:
            print(f"    {sympy.factor(h)} != 0")

    const_eqns = [e for e in eqns if is_constant_only(e)]
    trivial = psi_is_zero(C)

    print("\n  Constant-only equations (the projection onto the constants):")
    if const_eqns:
        for e in const_eqns:
            print("    ", e)
    else:
        print("     (none -- the component is generic in the constants)")

    if trivial:
        verdict = "TRIVIAL (Psi = 0): carries no solution, contributes no stratum"
    elif ansatz_only or membership_eqns is None:
        verdict = "not classified (--ansatz-only / --no-membership)"
    else:
        J = ConstRing.ideal([to_const(_numer(e)) for e in const_eqns]) if const_eqns \
            else ConstRing.ideal(ConstRing.zero())
        radJ = J.radical()
        inside_V = all(to_const(p) in radJ for p in membership_eqns)
        if inside_V:
            verdict = ("MEMBERSHIP: every projection equation vanishes on this "
                       "cell,\n              so the whole ansatz family solves "
                       "the PDE here (stratum <= V)")
        else:
            verdict = ("PARTIAL-SOLUTION STRATUM: nontrivial solutions, but some "
                       "projection\n              equation survives -- only a "
                       "proper sub-family solves the PDE")
        summary.append((i, const_eqns, inside_V))

    print(f"\n  Verdict: {verdict}")

# ---------------------------------------------------------------------------

print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)

if ansatz_only:
    print(f"RG on the ansatz alone returned {len(components)} component(s).")
elif membership_eqns is None:
    print(f"RG on the combined system returned {len(components)} component(s); "
          "no classification requested.")
else:
    n_mem = sum(1 for _, _, ok in summary if ok)
    n_par = sum(1 for _, _, ok in summary if not ok)
    n_triv = len(components) - len(summary)
    print(f"Combined system A u {{P}}: {len(components)} component(s)")
    print(f"  {n_triv:2d} trivial (Psi = 0)")
    print(f"  {n_mem:2d} inside V        -- the membership locus")
    print(f"  {n_par:2d} partial strata  -- proper sub-families solving the PDE")
    print(f"\nV itself (projection route) has {len(membership_primes)} minimal prime(s).")
    print("\npi_c(Sol(A) n Sol(P))  =  V  u  {partial strata}   [eq. (3) of the note]")
    if use_basefield and n_par == 0 and n_mem == 0:
        print("\nWith the constants in the base field every nonzero constant")
        print("polynomial is a unit, so RG could not split on any parameter locus:")
        print("the projection equations p_alpha(c) = 0 are inconsistent and the")
        print("only surviving component is the trivial one.  Re-run without")
        print("--basefield to see the strata (and read the README on the cliff).")
    elif use_basefield and n_par:
        print("\nNOTE: --basefield returned a nontrivial stratum.  That is a")
        print("surprise worth chasing: with the constants inverted RG should not")
        print("be able to split on a parameter locus at all.")

stamp("done")
