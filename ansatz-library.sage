# -*- mode: python -*-
#
# ansatz-library.sage
# ------------------------------------------------------------------------
# The ansaetze of the NewMethod / helium project, written directly in the
# DIFFERENTIAL-ALGEBRA formulation: each ansatz is a set of differential-
# polynomial equations over a differential-Thomas ranking, rather than the
# symbolic-ring substitution rules of ~/helium/helium.sage.
#
# WHY A SEPARATE FILE.  helium.sage encodes an ansatz as (a) a symbolic trial
# solution Psi and (b) hand-written `subs` rules that Ritt-reduce higher ODE
# derivatives -- e.g. ansatz 8 literally comments "a limitation of the program
# is that I have to manually calculate DD[0,0](Zeta)(B)".  That representation
# is fine for helium.sage's own reduce-then-project pipeline but it is NOT the
# object a differential-Thomas decomposition consumes.  Here every ansatz is a
# list of differential polynomials in an explicit jet ring, which is exactly
# what dt.differential_thomas_decomposition wants, and the formulation is
# uniform: an implicit (differentially triangular) family, leaders and all.
#
# WHAT AN ANSATZ IS, in this formulation.  A trial solution Psi built from one
# or more unknown "ODE functions" of one or more "inner variables", plus the
# defining relations:
#
#   * chain rules      Psi[c] - DPsi*v[c]           (per coordinate c)
#                      DPsi[c] - DDPsi*v[c]          (for a 2nd-order ODE)
#   * the ODE           (a0+a1*v)*DDPsi + ... = 0    (a diff. poly in the jets)
#   * inner-variable    v - (v1*x + v2*y + ...)      (defines the leader v)
#   * algebraic         r^2 - (x^2+y^2+z^2)          (root / extension relations)
#
# Under the ranking (jets high, params low) the leaders are DDPsi, DPsi's
# spatial derivatives, v and r; the free (parametric) jets are Psi and DPsi.
# Homogenization (helium.sage's device for forcing polynomials nonzero) is
# DROPPED: differential Thomas splits on those initials/separants itself, so
# ansatz 5 and 5.01 -- distinguished in helium.sage only by homogenization --
# are the same ansatz here.
#
# build_problem(pde_name, ansatz) returns everything the solver needs:
#   dict(R, rk, coords, roots, jets, params, ansatz_eqs, pconst, pde)
#
# Author: Brent Baccala (AI assistant: Claude).  July 2026.

import os, re, sys
from itertools import combinations_with_replacement
from collections import Counter

# native-Sage DifferentialThomas port (same dependency as joca-thomas-native-dt)
sys.path.insert(0, os.path.expanduser('~/DifferentialThomas-sage'))
sys.path.insert(0, os.path.expanduser('~/sage-differential-polynomial/src'))
import differentialthomas as dt


# ==========================================================================
# small helpers
# ==========================================================================
def _mono_str(combo):
    """A monomial (tuple of generator-name strings) as a ring-parser string."""
    if not combo:
        return '1'
    c = Counter(combo)
    return '*'.join(g if e == 1 else '%s^%d' % (g, e) for g, e in c.items())


def trial(base, gens, degree, constant=True, start=0, roots=()):
    """A native trial polynomial: sum of param*monomial over `gens` up to
    `degree`, roots capped at exponent 1 (they are square roots).  Returns
    (param_names, poly_string).  This is the differential-algebra analogue of
    helium.sage's trial_polynomial, but it emits a ring-parser string in the
    jet ring instead of a Symbolic-Ring expression."""
    mindeg = 0 if constant else 1
    terms = []
    for d in range(mindeg, degree + 1):
        for combo in combinations_with_replacement(gens, d):
            if all(combo.count(rt) < 2 for rt in roots):
                terms.append(combo)
    params = ['%s%d' % (base, i) for i in range(start, start + len(terms))]
    pieces = []
    for p, t in zip(params, terms):
        pieces.append(p if not t else '%s*%s' % (p, _mono_str(t)))
    return params, ' + '.join(pieces)


def _chain_rules(order, coords, inner='v'):
    """Psi[c]-DPsi*v[c], DPsi[c]-DDPsi*v[c], ... one per coordinate."""
    tower = ['D' * i + 'Psi' for i in range(order + 1)]
    eqs = []
    for lvl in range(order):
        lo, hi = tower[lvl], tower[lvl + 1]
        for c in coords:
            eqs.append('%s[%s] - %s*%s[%s]' % (lo, c, hi, inner, c))
    return eqs, tower


# ==========================================================================
# coordinate systems (one per PDE)
# ==========================================================================
def coordinate_system(pde_name):
    """(coords, roots) where roots = [(name, radicand_string), ...]."""
    if pde_name == 'hydrogen':
        return ['x', 'y', 'z'], [('r', 'x^2 + y^2 + z^2')]
    if pde_name == 'helium':
        # spherically-symmetric Nakatsuji S-state coordinates: no roots.
        return ['R1', 'R2', 'R12'], []
    raise ValueError("unknown pde %r" % pde_name)


# ==========================================================================
# the ansatz listing
# ==========================================================================
# Each entry returns dict(order, V, ODE, params, extra) where
#   V     : inner-variable polynomial string (defines jet `v`)
#   ODE   : the ODE as a differential polynomial in Psi/DPsi/DDPsi and v
#   params: ALL constant parameters introduced (v-coeffs then ODE-coeffs)
#   extra : any additional algebraic relations (e.g. algebraic extensions)
#
# NOTE.  These are the single-ODE-function, polynomial-inner-variable ansaetze
# (helium.sage's 4,5,5.1,5.2,5.3,8,9,10) -- the family that includes the
# validated hydrogen ansatz 5.  The product / rational-argument / algebraic-
# extension / nested-ODE ansaetze (1,2,3,6,7,11,12,13) each need their own
# differential-algebra template and are added next; see the stubs below.

def ansatz_spec(ansatz, coords, roots):
    gens = coords + [rn for rn, _ in roots]
    rset = tuple(rn for rn, _ in roots)

    if ansatz == 5:
        # 2nd-order ODE, linear coeffs, linear inner variable.
        # (helium.sage 5 / 5.01 -- homogenization dropped; Thomas splits it.)
        vp, V = trial('v', gens, 1, constant=False, start=1, roots=rset)
        ap, A = trial('a', ['v'], 1)
        bp, B = trial('b', ['v'], 1)
        cp, C = trial('c', ['v'], 1)
        ODE = '(%s)*DDPsi + (%s)*DPsi + (%s)*Psi' % (A, B, C)
        return dict(order=2, V=V, ODE=ODE, params=vp + ap + bp + cp, extra=[])

    if ansatz == 5.1:
        # 2nd-order ODE, linear coeffs, QUADRATIC inner variable.
        vp, V = trial('v', gens, 2, constant=False, roots=rset)
        dp, D = trial('d', ['v'], 1)
        mp, M = trial('m', ['v'], 1)
        np_, N = trial('n', ['v'], 1)
        ODE = '(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)
        return dict(order=2, V=V, ODE=ODE, params=vp + dp + mp + np_, extra=[])

    if ansatz == 5.2:
        # 2nd-order ODE, QUADRATIC coeffs, linear inner variable.
        vp, V = trial('v', gens, 1, constant=False, roots=rset)
        dp, D = trial('d', ['v'], 2)
        mp, M = trial('m', ['v'], 2)
        np_, N = trial('n', ['v'], 2)
        ODE = '(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)
        return dict(order=2, V=V, ODE=ODE, params=vp + dp + mp + np_, extra=[])

    if ansatz in (5.3, 10):
        # 2nd-order ODE, QUADRATIC coeffs, QUADRATIC inner variable.
        # (helium.sage 5.3 and 10 coincide in this formulation.)
        vp, V = trial('v', gens, 2, constant=False, roots=rset)
        dp, D = trial('d', ['v'], 2)
        mp, M = trial('m', ['v'], 2)
        np_, N = trial('n', ['v'], 2)
        ODE = '(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)
        return dict(order=2, V=V, ODE=ODE, params=vp + dp + mp + np_, extra=[])

    if ansatz == 8:
        # 1st-order ODE, linear coeffs, linear inner variable: M*DPsi - N*Psi = 0.
        vp, V = trial('v', gens, 1, constant=False, roots=rset)
        mp, M = trial('m', ['v'], 1)
        np_, N = trial('n', ['v'], 1)
        ODE = '(%s)*DPsi - (%s)*Psi' % (M, N)
        return dict(order=1, V=V, ODE=ODE, params=vp + mp + np_, extra=[])

    if ansatz == 9:
        # 1st-order ODE, constant rate: DPsi - n0*Psi = 0 (pure exponential in V).
        vp, V = trial('v', gens, 1, constant=False, roots=rset)
        np_, N = trial('n', ['v'], 0)          # N = n0
        ODE = 'DPsi - (%s)*Psi' % N
        return dict(order=1, V=V, ODE=ODE, params=vp + np_, extra=[])

    raise NotImplementedError(
        "ansatz %s not yet in the differential-algebra library.\n"
        "  product ansaetze (1,2,3: Psi=A(coords)*F(V))     -> add F as a jet, "
        "Psi-A*F as a defining eq;\n"
        "  rational argument (6,7: F(B/C))                  -> add inner w with "
        "C*w-B=0;\n"
        "  algebraic extension (11: gamma)                  -> add g with its "
        "minimal polynomial as an `extra` eq;\n"
        "  nested ODEs (12,13)                              -> two ODE functions, "
        "two inner variables." % ansatz)


# ==========================================================================
# the two Hamiltonians -> cleared native PDE
# ==========================================================================
# The Hamiltonian is written inline (clean physics, not helium.sage's ansatz
# code) and turned into a denominator-cleared differential polynomial in the
# coordinate jets Psi[c], Psi[c1,c2].

def _hamiltonian(pde_name):
    """Return (PsiF, expr) where expr = H(Psi)-E*Psi as a Symbolic expression
    in a genuine function PsiF of the coordinates."""
    Evar = var('E')
    if pde_name == 'hydrogen':
        x, y, z = var('x y z')
        # r is an independent ring symbol here; the relation r^2 = x^2+y^2+z^2 is
        # imposed separately as an ansatz equation (exactly as in the driver), so
        # the potential -1/r stays polynomial-after-clearing without a sqrt.
        r = var('r')
        PsiF = function('PsiF')(x, y, z)
        H = -1/2 * (diff(PsiF, x, 2) + diff(PsiF, y, 2) + diff(PsiF, z, 2)) \
            - (1/r) * PsiF
        rootsub = {}
    else:
        R1, R2, R12 = var('R1 R2 R12')
        PsiF = function('PsiF')(R1, R2, R12)
        H = (-1/2 * sum(diff(PsiF, Ri, 2) + 2/Ri * diff(PsiF, Ri) for Ri in (R1, R2))
             - (diff(PsiF, R12, 2) + 2/R12 * diff(PsiF, R12))
             - (R1**2 + R12**2 - R2**2)/(2*R1*R12) * diff(diff(PsiF, R12), R1)
             - (R2**2 + R12**2 - R1**2)/(2*R2*R12) * diff(diff(PsiF, R12), R2)
             - sum(2/Ri for Ri in (R1, R2)) * PsiF + 1/R12 * PsiF)
        rootsub = {}
    expr = H - Evar * PsiF
    for k, v_ in rootsub.items():
        expr = expr.subs({k: v_})
    return PsiF, expr


def build_pde_string(pde_name, coords):
    """H(Psi)-E*Psi, jetified and denominator-cleared, as a ring-parser string."""
    PsiF, expr = _hamiltonian(pde_name)
    cv = [var(c) for c in coords]
    order_key = {c: i for i, c in enumerate(coords)}   # derivation-declaration order

    # jetify: replace each derivative of PsiF (order <= 2, incl. mixed) and PsiF
    # itself with a fresh symbol carrying its jet string.  Derivative indices are
    # written in coordinate-declaration order (what the BLAD ranking expects).
    subst, jetmap = {}, {}
    n = 0
    for c1 in cv:
        for c2 in cv:
            d = diff(PsiF, c1, c2)
            if d in subst:
                continue
            t = var('JJ%d' % n); n += 1
            subst[d] = t
            idx = sorted([str(c1), str(c2)], key=lambda s: order_key[s])
            jetmap[str(t)] = 'Psi[%s]' % ','.join(idx)
    for c1 in cv:
        d = diff(PsiF, c1)
        t = var('JJ%d' % n); n += 1
        subst[d] = t
        jetmap[str(t)] = 'Psi[%s]' % str(c1)
    t0 = var('JJ%d' % n)
    subst[PsiF] = t0
    jetmap[str(t0)] = 'Psi'

    expr = expr.subs(subst)
    cleared = expr.numerator()               # clears coordinate/root denominators
    assert (cleared / expr).denominator() == 1 or expr.denominator() != 1
    s = str(cleared.expand())
    for k in sorted(jetmap, key=len, reverse=True):
        s = re.sub(r'(?<![\w\[])' + re.escape(k) + r'(?![\w\]])', jetmap[k], s)
    return s.replace('**', '^')


# ==========================================================================
# assemble the full problem
# ==========================================================================
def build_problem(pde_name, ansatz):
    coords, roots = coordinate_system(pde_name)
    spec = ansatz_spec(ansatz, coords, roots)
    order = spec['order']
    params = spec['params']

    chain, tower = _chain_rules(order, coords)
    root_eqs = ['%s^2 - (%s)' % (rn, rad) for rn, rad in roots]

    ansatz_eqs_str = (chain
                      + [spec['ODE']]
                      + ['v - (%s)' % spec['V']]
                      + root_eqs
                      + list(spec['extra']))

    # jets high->low: Psi-tower (DDPsi,DPsi,Psi), then v, then roots.
    jets = list(reversed(tower)) + ['v'] + [rn for rn, _ in roots]
    IVAR = coords
    DVAR = jets + ['E'] + params
    rk = dt.compute_ranking(IVAR, DVAR)
    R = rk.ring

    ansatz_eqs = [R(s) for s in ansatz_eqs_str]
    pconst = [R('%s[%s]' % (p, c)) for p in params for c in coords]
    pde = R(build_pde_string(pde_name, coords))

    # inner-variable coefficients: the params appearing in V.  If a solution
    # variety forces ALL of these to zero then v == 0 identically and the ansatz
    # has collapsed to a constant -- a degenerate (non-)solution.
    v_toks = set(re.findall(r'[A-Za-z]\w*', spec['V']))
    v_params = [p for p in params if p in v_toks]

    return dict(R=R, rk=rk, coords=coords, roots=roots, jets=jets,
                tower=tower, order=order, params=params, v_params=v_params,
                ansatz_eqs=ansatz_eqs, pconst=pconst, pde=pde,
                ansatz_eqs_str=ansatz_eqs_str)
