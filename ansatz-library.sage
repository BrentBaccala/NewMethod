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


def _log_relations(coords):
    """Cleared defining relations for the transcendental jet
    L = log(hyperradius^2); leaders L[c], initial = the radicand.  L has NO
    order-0 relation (log is transcendental; the free order-0 jet is exactly the
    additive integration constant of the log).  See the log-hyperradius ansaetze
    (17/17.1/18/19) and ~/project/reports/helium-new-ansatze.md.  NOTE: for
    helium L = log(R1^2+R2^2) = 2 log s (the triple-coalescence hyperradius); the
    hydrogen analogue L = log(x^2+y^2+z^2) = 2 log r is included so the ansaetze
    that are otherwise coordinate-agnostic can be smoke-tested on hydrogen."""
    if coords == ['R1', 'R2', 'R12']:
        rad, grad = 'R1^2 + R2^2', {'R1': '2*R1', 'R2': '2*R2', 'R12': '0'}
    else:                                          # hydrogen
        rad, grad = 'x^2 + y^2 + z^2', {'x': '2*x', 'y': '2*y', 'z': '2*z'}
    return ['(%s)*L[%s] - (%s)' % (rad, c, grad[c]) for c in coords]


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
# Entries are in NUMERICAL order by ansatz number.  Two return shapes:
#   Zeta / single-ODE-function family (5,5.1,5.2,5.3,8,9,10,19 and the
#     coeff-ring variants 14/15/16) -- dict(order, V, ODE, params, extra) where
#       V     : inner-variable polynomial string (defines jet `v`)
#       ODE   : the ODE as a differential polynomial in Psi/DPsi/DDPsi and v
#       params: ALL constant parameters introduced (v-coeffs then ODE-coeffs)
#       extra : any additional algebraic relations (e.g. algebraic extensions)
#   every other family (product 1/1.1/2/3, rational 6/7, nested 12, algext 13,
#     log 17/18, exp-Fock 20/20.1) -- dict(kind, jets_dep, equations, params,
#     v_params, amp_params); its differential-polynomial equations are inline.
# Each family's template is documented in a header comment at its first entry.

def ansatz_spec(ansatz, coords, roots):
    gens = coords + [rn for rn, _ in roots]
    rset = tuple(rn for rn, _ in roots)

    # ----- PRODUCT / EXPONENTIAL template: Psi = A(coords) * F(inner) ---------
    # The solution is a coordinate polynomial A times an ODE function F of an
    # inner variable.  Psi is defined by  Psi - A*F = 0  (leader Psi); F has its
    # own chain rule / ODE.  Unlike the Zeta(V) family the chain rule carries the
    # A-factor: differentiating Psi - A*F gives Psi[c] = A_c*F + A*F[c] (the
    # engine forms A_c since A is built from the independent coordinates).  The
    # inner-variable jet is named B (exp/Chi) or C (log); its coefficients are
    # v_params, so B==0 (exponential/log collapses to a constant) reads as the
    # degenerate case.  Homogenization is dropped (Thomas splits those loci).

    if ansatz == 1:
        # Psi = A * Phi,  Phi = exp(B):  Phi' = Phi, so the chain rule is
        # Phi[c] - Phi*B[c]  (no separate derivative jet).
        ap, A = trial('a', gens, 1, roots=rset)                    # A(coords)
        bp, B = trial('b', gens, 1, constant=False, roots=rset)    # inner B
        eqs = (['Psi - (%s)*Phi' % A]
               + ['Phi[%s] - Phi*B[%s]' % (c, c) for c in coords]
               + ['B - (%s)' % B])
        return dict(kind='product', jets_dep=['Psi', 'Phi', 'B'],
                    equations=eqs, params=ap + bp, v_params=bp, amp_params=ap)

    if ansatz == 1.1:
        # Hylleraas-type cusp with quadratic amplitude: Psi = A(deg 2)*exp(B).
        # No new template -- ansatz 1 with a quadratic A.  (report ansatz 1.1)
        ap, A = trial('a', gens, 2, roots=rset)
        bp, B = trial('b', gens, 1, constant=False, roots=rset)
        eqs = (['Psi - (%s)*Phi' % A]
               + ['Phi[%s] - Phi*B[%s]' % (c, c) for c in coords]
               + ['B - (%s)' % B])
        return dict(kind='product', jets_dep=['Psi', 'Phi', 'B'],
                    equations=eqs, params=ap + bp, v_params=bp, amp_params=ap)

    if ansatz == 2:
        # Psi = A * Xi,  Xi = log(C):  Xi' = 1/C, so the chain rule cleared of
        # the denominator is  C*Xi[c] - C[c].
        ap, A = trial('a', gens, 1, roots=rset)                    # A(coords)
        cp, C = trial('c', gens, 1, constant=False, roots=rset)    # inner C
        eqs = (['Psi - (%s)*Xi' % A]
               + ['C*Xi[%s] - C[%s]' % (c, c) for c in coords]
               + ['C - (%s)' % C])
        return dict(kind='product', jets_dep=['Psi', 'Xi', 'C'],
                    equations=eqs, params=ap + cp, v_params=cp, amp_params=ap)

    if ansatz == 3:
        # Psi = A * Chi,  Chi a 2nd-order ODE function of B with COORDINATE-
        # polynomial coefficients (helium.sage's "weird second-order mess"):
        #   pC*Chi'' - pD*Chi' - pF*Chi - pG = 0.
        ap, A = trial('a', gens, 1, roots=rset)
        bp, B = trial('b', gens, 1, constant=False, roots=rset)    # inner B
        cp, pC = trial('c', gens, 1, roots=rset)                   # ODE coeffs
        dp, pD = trial('d', gens, 1, roots=rset)                   # (coord polys)
        fp, pF = trial('f', gens, 1, roots=rset)
        gp, pG = trial('g', gens, 1, roots=rset)
        eqs = (['Psi - (%s)*Chi' % A]
               + ['Chi[%s] - DChi*B[%s]' % (c, c) for c in coords]
               + ['DChi[%s] - DDChi*B[%s]' % (c, c) for c in coords]
               + ['(%s)*DDChi - (%s)*DChi - (%s)*Chi - (%s)' % (pC, pD, pF, pG)]
               + ['B - (%s)' % B])
        return dict(kind='product', jets_dep=['Psi', 'DDChi', 'DChi', 'Chi', 'B'],
                    equations=eqs, params=ap + bp + cp + dp + fp + gp,
                    v_params=bp, amp_params=ap)

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

    # ----- RATIONAL-ARGUMENT template: Psi = Zeta(w), w = B/C -----------------
    # The inner argument is a rational function of the coordinates.  Rather than
    # carry a denominator, introduce w as a jet defined by the cleared relation
    # C*w - B = 0 (leader w, separant C -- so C != 0 is a Thomas inequation);
    # everything else is the Zeta(V) family with w in place of v.  B/C is
    # invariant under (B,C)->(lambda B, lambda C), an extra scaling dimension we
    # let Thomas carry (helium.sage suppressed it with homogenization).
    if ansatz in (6, 7):
        d_bc = 1 if ansatz == 6 else 2          # degree of B, C and of the ODE coeffs
        bp, B = trial('b', gens, d_bc, roots=rset)      # numerator   (with constant)
        cp, C = trial('c', gens, d_bc, roots=rset)      # denominator (with constant)
        dp, D = trial('d', ['w'], d_bc)                 # ODE coeffs in w = B/C
        mp, M = trial('m', ['w'], d_bc)
        np_, N = trial('n', ['w'], d_bc)
        eqs = (['Psi[%s] - DPsi*w[%s]' % (c, c) for c in coords]
               + ['DPsi[%s] - DDPsi*w[%s]' % (c, c) for c in coords]
               + ['(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)]
               + ['(%s)*w - (%s)' % (C, B)])            # C*w - B = 0  (w = B/C)
        return dict(kind='rational',
                    jets_dep=['DDPsi', 'DPsi', 'Psi', 'w'],
                    equations=eqs, params=bp + cp + dp + mp + np_,
                    v_params=bp)                         # B==0 -> w=0 -> degenerate

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

    # ----- NESTED-ODE template: Psi = Zeta(V), V depends on Theta = Theta(U) --
    # Two coupled ODE functions.  The inner function Theta solves an ODE in the
    # inner-inner variable U; the OUTER inner variable V is a polynomial in the
    # coordinates AND Theta, so the outer solution Psi = Zeta(V) is a function of
    # a function.  Psi/DPsi/DDPsi are the outer ODE function (= Zeta); Theta/
    # DTheta/DDTheta the inner one.  Reduction of the PDE cascades down both
    # towers: Psi[c] -> DPsi*V[c], V[c] -> (coords) + v_theta*Theta[c],
    # Theta[c] -> DTheta*U[c].
    if int(ansatz) == 12:
        deg = {12: (1, 1, 1, 1), 12.1: (2, 1, 1, 1), 12.2: (1, 2, 1, 1),
               12.3: (1, 1, 2, 1), 12.4: (1, 1, 1, 2)}[ansatz]
        maxdeg_v, maxdeg_u, ode_v, ode_u = deg

        up, U = trial('u', gens, maxdeg_u, constant=False, roots=rset)
        ap, A = trial('a', ['U'], ode_u)                # inner ODE coeffs in U
        bp, B = trial('b', ['U'], ode_u)
        cp, C = trial('c', ['U'], ode_u)
        # outer inner-variable V ranges over coords/roots AND the inner soln Theta
        vp, V = trial('v', gens + ['Theta'], maxdeg_v, constant=False, roots=rset)
        dp, D = trial('d', ['V'], ode_v)                # outer ODE coeffs in V
        mp, M = trial('m', ['V'], ode_v)
        np_, N = trial('n', ['V'], ode_v)
        eqs = (['Psi[%s] - DPsi*V[%s]' % (c, c) for c in coords]
               + ['DPsi[%s] - DDPsi*V[%s]' % (c, c) for c in coords]
               + ['(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)]
               + ['V - (%s)' % V]
               + ['Theta[%s] - DTheta*U[%s]' % (c, c) for c in coords]
               + ['DTheta[%s] - DDTheta*U[%s]' % (c, c) for c in coords]
               + ['(%s)*DDTheta - (%s)*DTheta - (%s)*Theta' % (A, B, C)]
               + ['U - (%s)' % U])
        return dict(kind='nested',
                    jets_dep=['DDPsi', 'DPsi', 'Psi', 'V',
                              'DDTheta', 'DTheta', 'Theta', 'U'],
                    equations=eqs, params=up + ap + bp + cp + vp + dp + mp + np_,
                    v_params=vp)

    # ----- ALGEBRAIC-EXTENSION template (helium.sage calls it 13) -------------
    # Psi = Zeta(V) with a 2nd-order ODE whose coefficients live in a quadratic
    # algebraic extension: g satisfies A*g^2 + B*g + C = 0 (A,B,C polys in V), and
    # the ODE coeffs D,M,N are polynomials in V and g.  g is a differential
    # indeterminate with an algebraic (order-0) defining relation; its derivatives
    # follow from prolonging that relation (separant 2*A*g+B).  This is the same
    # template ansatz 11 needs.
    if int(ansatz) == 13:
        deg = {13: (1, 1, 1), 13.1: (2, 1, 1), 13.2: (1, 2, 1), 13.3: (1, 1, 2),
               13.4: (2, 2, 1), 13.5: (2, 1, 2), 13.6: (2, 2, 2)}[ansatz]
        maxdeg_v, ode_deg, alg_deg = deg

        vp, V = trial('v', gens, maxdeg_v, constant=False, roots=rset)
        ap, A = trial('a', ['V'], alg_deg)              # min-poly coeffs in V
        bp, B = trial('b', ['V'], alg_deg)
        cp, C = trial('c', ['V'], alg_deg)
        dp, D = trial('d', ['V', 'g'], ode_deg)         # ODE coeffs in V and g
        mp, M = trial('m', ['V', 'g'], ode_deg)
        np_, N = trial('n', ['V', 'g'], ode_deg)
        eqs = (['Psi[%s] - DPsi*V[%s]' % (c, c) for c in coords]
               + ['DPsi[%s] - DDPsi*V[%s]' % (c, c) for c in coords]
               + ['(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)]
               + ['V - (%s)' % V]
               + ['(%s)*g^2 + (%s)*g + (%s)' % (A, B, C)])   # gamma minimal poly
        return dict(kind='algext',
                    jets_dep=['DDPsi', 'DPsi', 'Psi', 'g', 'V'],
                    equations=eqs, params=vp + dp + mp + np_ + ap + bp + cp,
                    v_params=vp)

    # ----- COEFFICIENT-RING EXTENSION templates (NewSol.tex ansatz collection) -
    # Psi = Zeta(v) with a 2nd-order ODE  D*DDPsi - M*DPsi - N*Psi = 0  whose
    # coefficients D,M,N live in an EXTENSION of the coefficient ring.  Ansatz 13
    # (already above) is the ALGEBRAIC case (g root of a quadratic).  14/15/16
    # below add the remaining cases from NewSol.tex's ansatz collection (the
    # section removed in Papers commit fc28cfca, Oct 2024).  These are my
    # interpretation of the schematic paper diagrams -- base cases only (degree 1
    # where possible), implemented faithfully; get-them-to-parse is the goal.

    if int(ansatz) == 14:
        # NewSol.tex: EXPONENTIAL element in the coefficient ring.  An exponential
        # t = exp(e0*v) sits in the ODE coefficient ring (t multiplies the
        # coefficients D,M,N, NOT Psi), with defining relation t[c] - e0*t*v[c]
        # per coordinate -- the same shape as ansatz 1's Phi[c]-Phi*B[c].  ONE
        # rate parameter e0; D,M,N are degree-1 polys in v and t.
        vp, V = trial('v', gens, 1, constant=False, roots=rset)
        dp, D = trial('d', ['v', 't'], 1)
        mp, M = trial('m', ['v', 't'], 1)
        np_, N = trial('n', ['v', 't'], 1)
        eqs = (['Psi[%s] - DPsi*v[%s]' % (c, c) for c in coords]
               + ['DPsi[%s] - DDPsi*v[%s]' % (c, c) for c in coords]
               + ['t[%s] - e0*t*v[%s]' % (c, c) for c in coords]
               + ['(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)]
               + ['v - (%s)' % V])
        return dict(kind='expext',
                    jets_dep=['DDPsi', 'DPsi', 'Psi', 't', 'v'],
                    equations=eqs, params=vp + ['e0'] + dp + mp + np_,
                    v_params=vp, amp_params=[])

    if int(ansatz) == 15:
        # NewSol.tex: 2nd-order HOLONOMIC element in the coefficient ring.  A
        # holonomic element t (jets t,Dt,DDt) solves its OWN 2nd-order ODE in v,
        # and Psi's ODE coefficients D,M,N are degree-1 polys in v, t and Dt.
        # Heavy (7 dependent jets) -- inherent to a holonomic coefficient ring.
        vp, V = trial('v', gens, 1, constant=False, roots=rset)
        pp, P = trial('p', ['v'], 1)                    # t's ODE coeffs (in v)
        qp, Q = trial('q', ['v'], 1)
        sp, S = trial('s', ['v'], 1)
        dp, D = trial('d', ['v', 't', 'Dt'], 1)         # Psi's ODE coeffs
        mp, M = trial('m', ['v', 't', 'Dt'], 1)
        np_, N = trial('n', ['v', 't', 'Dt'], 1)
        eqs = (['Psi[%s] - DPsi*v[%s]' % (c, c) for c in coords]
               + ['DPsi[%s] - DDPsi*v[%s]' % (c, c) for c in coords]
               + ['t[%s] - Dt*v[%s]' % (c, c) for c in coords]
               + ['Dt[%s] - DDt*v[%s]' % (c, c) for c in coords]
               + ['(%s)*DDt - (%s)*Dt - (%s)*t' % (P, Q, S)]
               + ['(%s)*DDPsi - (%s)*DPsi - (%s)*Psi' % (D, M, N)]
               + ['v - (%s)' % V])
        return dict(kind='holoext',
                    jets_dep=['DDPsi', 'DPsi', 'Psi', 'DDt', 'Dt', 't', 'v'],
                    equations=eqs, params=vp + pp + qp + sp + dp + mp + np_,
                    v_params=vp, amp_params=[])

    if int(ansatz) == 16:
        # NewSol.tex / helium.sage coded ansatz 16: an algebraic root nested
        # BELOW the extension.  Psi = Zeta(v) with a 2nd-order ODE; the novelty is
        # in the BASE -- an algebraic element g = sqrt(radicand) with
        # g^2 - radicand = 0, and the inner variable v ranging over the
        # coordinates AND g.  The DECIMAL encodes three degree bounds
        # (radicand / v / ODE-coeffs), matching the NewSol.tex Ansatz 16 table.
        # In particular 16.3 = (2,1,1) and 16.6 = (2,2,2) use a QUADRATIC radicand
        # -- g is the square root of a degree-2 polynomial, i.e. a genuinely new
        # coordinate of the same shape as hydrogen's r = sqrt(x^2+y^2+z^2), the
        # ingredient that let ansatz 5 find the J0 Bessel solution.  (helium.sage's
        # homogenized 16.31/16.61 collapse onto 16.3/16.6 here -- homogenization
        # is dropped; Thomas splits those loci itself.)
        deg16 = {16: (1, 1, 1), 16.1: (1, 2, 1), 16.2: (1, 1, 2),
                 16.3: (2, 1, 1), 16.4: (1, 2, 2), 16.5: (2, 1, 2),
                 16.6: (2, 2, 2)}
        rad_deg, v_deg, ode_deg = deg16[ansatz]
        vp, V = trial('v', gens + ['g'], v_deg, constant=False, roots=rset)
        kp, RAD = trial('k', gens, rad_deg, roots=rset)   # radicand: g = sqrt(RAD)
        ap, A = trial('a', ['v'], ode_deg)
        bp, B = trial('b', ['v'], ode_deg)
        cp, C = trial('c', ['v'], ode_deg)
        ODE = '(%s)*DDPsi + (%s)*DPsi + (%s)*Psi' % (A, B, C)
        eqs = (['Psi[%s] - DPsi*v[%s]' % (c, c) for c in coords]
               + ['DPsi[%s] - DDPsi*v[%s]' % (c, c) for c in coords]
               + [ODE]
               + ['v - (%s)' % V]
               + ['g^2 - (%s)' % RAD])                  # gamma: g^2 = radicand
        return dict(kind='algbase',
                    jets_dep=['DDPsi', 'DPsi', 'Psi', 'v', 'g'],
                    equations=eqs, params=vp + kp + ap + bp + cp,
                    v_params=vp, amp_params=[])

    # ----- LOG-HYPERRADIUS templates (Bartlett-Fock log; report ansaetze) ------
    # Psi carrying the log-non-analyticity L = log(R1^2+R2^2) = 2 log s admitted
    # via a transcendental jet with rational (cleared) derivatives; see
    # ~/project/reports/helium-new-ansatze.md for the full derivation.  17/17.1
    # (loglin) and 18 (product+log) hard-code the helium hyperradius gradient in
    # _log_relations, so they are HELIUM-SPECIFIC (17.1 additionally hard-codes
    # (R1^2+R2^2) in the trial form); a hydrogen build-check is not expected to be
    # meaningful for these.  19 rides the Zeta family via the extra_jets hook.

    if int(ansatz) == 17:
        if ansatz == 17:
            # Fock-linear: Psi = A + B*L, A,B degree 1 -- the minimal form with a
            # genuine log s.  (report ansatz 14)
            ap, A = trial('a', gens, 1, roots=rset)
            bp, B = trial('b', gens, 1, roots=rset)
            eqs = (['Psi - ((%s) + (%s)*L)' % (A, B)]   # leader Psi
                   + _log_relations(coords))            # leaders L[c]
            return dict(kind='loglin', jets_dep=['Psi', 'L'],
                        equations=eqs, params=ap + bp,
                        v_params=bp,            # B == 0  -> log gone -> DEGENERATE
                        amp_params=ap + bp)     # A==B==0 -> Psi == 0 -> TRIVIAL
        else:                                           # ansatz 17.1
            # Fock-exact slot: Psi = A(deg 2) + b0*(R1^2+R2^2)*L -- the log pinned
            # to the exact O(s^2 log s) Fock slot.  (report ansatz 14.1)
            ap, A = trial('a', gens, 2, roots=rset)
            eqs = (['Psi - ((%s) + b0*(R1^2 + R2^2)*L)' % A]
                   + _log_relations(coords))
            return dict(kind='loglin', jets_dep=['Psi', 'L'],
                        equations=eqs, params=ap + ['b0'],
                        v_params=['b0'], amp_params=ap + ['b0'])

    if ansatz == 18:
        # Kato x Fock: Psi = exp(B)*(A0 + A1*L) -- cusp exponents times the first
        # Fock log.  Product template extended by the log jet.  (report ansatz 15)
        bp,  Bx = trial('b', gens, 1, constant=False, roots=rset)
        a0p, A0 = trial('a', gens, 1, roots=rset)
        a1p, A1 = trial('h', gens, 1, roots=rset)
        eqs = (['Psi - Phi*((%s) + (%s)*L)' % (A0, A1)]
               + ['Phi[%s] - Phi*B[%s]' % (c, c) for c in coords]
               + ['B - (%s)' % Bx]
               + _log_relations(coords))
        return dict(kind='product', jets_dep=['Psi', 'Phi', 'B', 'L'],
                    equations=eqs, params=bp + a0p + a1p,
                    v_params=a1p,             # A1 == 0 -> log gone (ansatz-1 land)
                    amp_params=a0p + a1p)     # A0==A1==0 -> Psi == 0

    if ansatz == 19:
        # Zeta family with a log-extended inner variable: v ranges over the coords
        # AND L = log(R1^2+R2^2), so Psi = Zeta(v1*R1+v2*R2+v3*R12+v4*L).  Uses the
        # extra_jets hook in build_problem's Zeta branch and passes the log
        # relations through the family's `extra` slot.  (report ansatz 16)
        vp, V = trial('v', gens + ['L'], 1, constant=False, roots=rset)
        ap, A = trial('a', ['v'], 1)
        bp, B = trial('b', ['v'], 1)
        cp, C = trial('c', ['v'], 1)
        ODE = '(%s)*DDPsi + (%s)*DPsi + (%s)*Psi' % (A, B, C)
        return dict(order=2, V=V, ODE=ODE, params=vp + ap + bp + cp,
                    extra=_log_relations(coords), extra_jets=['L'])

    if ansatz == 20:
        # EXPONENTIATED FOCK (Myers-Umrigar-Sethna-Morgan 1991, sec IV): the log
        # goes INSIDE the exponent.  Psi = exp(B + C*L) with L = log(R1^2+R2^2) =
        # 2 log s (the triple-coalescence hyperradius) and B, C polynomial in the
        # (R1,R2,R12) coordinates.  This is the first library entry to place the
        # log in the exponent: exp(C*L) = (R1^2+R2^2)^C is a variable-exponent
        # (x^x-type) factor whose Taylor series carries the UNBOUNDED (ln s)^p Fock
        # tower with a finite parameter set.  17/17.1/18 keep the log linear /
        # multiplicative and so can only ever represent (ln s)^1 -- which is why
        # they came back clean-negative: too holonomic to hold the tower.
        #   B (deg 2, non-const)  -- the log-free exponent: Kato cusp psi_{1,0}
        #     (linear -Z(R1+R2)+R12/2) plus the O(s^2) pieces psi_{2,0} - 1/2 psi_{1,0}^2.
        #   C (deg 2)             -- the Fock r^2 log-slot coefficient; the exact
        #     leading term is psi_{2,1} prop. (pi-2)/(3pi) Z (R1^2+R2^2-R12^2)
        #     = (pi-2)/(3pi) Z Y_{2,1}, so C must reach degree 2.
        # (R1,R2,R12) ARE the KS-rationalized coordinates -- the hyperspherical
        # sqrt(1-sin a cos t), sqrt(1+sin a) irrationalities of the angular Fock
        # coefficients are polynomial here, so L is the only transcendental jet.
        # Refs: MUSM 1991 eq (23); Liverts 2022 eq (8); Fournais et al. 2004/2009;
        # ~/project/reports/helium-new-ansatze.md.
        bp, Bx = trial('b', gens, 2, constant=False, roots=rset)   # log-free exponent
        cp, Cx = trial('c', gens, 2, roots=rset)                   # Fock-log coefficient
        eqs = (['Psi[%s] - Psi*(B[%s] + C[%s]*L + C*L[%s])' % (c, c, c, c)
                for c in coords]                       # Psi = exp(B + C*L)
               + ['B - (%s)' % Bx]
               + ['C - (%s)' % Cx]
               + _log_relations(coords))               # leaders L[c]
        return dict(kind='product', jets_dep=['Psi', 'B', 'C', 'L'],
                    equations=eqs, params=bp + cp,
                    v_params=cp,        # C == 0 -> log gone -> Kato-only (DEGENERATE)
                    amp_params=[])      # exp is never 0: no Psi==0 (TRIVIAL) mode

    if ansatz == 20.1:
        # EXPONENTIATED FOCK + polynomial amplitude: Psi = A * exp(B + C*L),
        # A degree 1 in the coordinates.  Ansatz 20 with a Hylleraas-style
        # prefactor (cf. 1 -> 1.1): the amplitude lets Psi carry a polynomial
        # node/bulk factor on top of the exp-Fock singular structure, the way
        # hydrogen 2s = (1 - Zr/2) e^{-Zr/2} carries its node.  With an amplitude
        # the Psi==0 (TRIVIAL) mode is back (A == 0).  Refs as ansatz 20.
        ap, Ax = trial('a', gens, 1, roots=rset)                   # polynomial amplitude
        bp, Bx = trial('b', gens, 2, constant=False, roots=rset)   # log-free exponent
        cp, Cx = trial('c', gens, 2, roots=rset)                   # Fock-log coefficient
        eqs = (['Psi - (%s)*Phi' % Ax]                             # Psi = A * exp(B+C*L)
               + ['Phi[%s] - Phi*(B[%s] + C[%s]*L + C*L[%s])' % (c, c, c, c)
                  for c in coords]
               + ['B - (%s)' % Bx]
               + ['C - (%s)' % Cx]
               + _log_relations(coords))
        return dict(kind='product', jets_dep=['Psi', 'Phi', 'B', 'C', 'L'],
                    equations=eqs, params=ap + bp + cp,
                    v_params=cp,        # C == 0 -> log gone -> DEGENERATE
                    amp_params=ap)      # A == 0 -> Psi == 0 -> TRIVIAL

    raise NotImplementedError(
        "ansatz %s not yet in the differential-algebra library.\n"
        "  algebraic extension (11: gamma)  -> same template as 13." % ansatz)


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
    params = spec['params']
    root_eqs = ['%s^2 - (%s)' % (rn, rad) for rn, rad in roots]

    if 'equations' in spec:
        # explicit-equation families (product, nested, algext): the spec has
        # already assembled the differential-polynomial equations and jet list.
        jets_dep = spec['jets_dep']
        ansatz_eqs_str = list(spec['equations']) + root_eqs
        v_params = spec['v_params']
        amp_params = spec.get('amp_params', [])
        tower = None
    else:
        # Zeta(V) single-ODE-function family.
        order = spec['order']
        chain, tower = _chain_rules(order, coords)
        ansatz_eqs_str = (chain
                          + [spec['ODE']]
                          + ['v - (%s)' % spec['V']]
                          + list(spec['extra']))
        ansatz_eqs_str += root_eqs
        # extra_jets hook: a Zeta-family ansatz may extend the inner variable with
        # an extra differential jet (e.g. ansatz 19's log jet L) whose defining
        # relations arrive via the `extra` slot; rank it just below v.
        jets_dep = list(reversed(tower)) + ['v'] + list(spec.get('extra_jets', []))
        # inner-variable coefficients: params appearing in V.
        v_toks = set(re.findall(r'[A-Za-z]\w*', spec['V']))
        v_params = [p for p in params if p in v_toks]
        amp_params = []          # Zeta family: Psi is the free jet (no amplitude)

    # jets high->low: dependent jets, then roots.
    jets = list(jets_dep) + [rn for rn, _ in roots]
    IVAR = coords
    DVAR = jets + ['E'] + params
    rk = dt.compute_ranking(IVAR, DVAR)
    R = rk.ring

    ansatz_eqs = [R(s) for s in ansatz_eqs_str]
    pconst = [R('%s[%s]' % (p, c)) for p in params for c in coords]
    pde = R(build_pde_string(pde_name, coords))

    return dict(R=R, rk=rk, coords=coords, roots=roots, jets=jets,
                tower=tower, order=spec.get('order'), params=params,
                v_params=v_params, amp_params=amp_params, ansatz_eqs=ansatz_eqs,
                pconst=pconst, pde=pde, ansatz_eqs_str=ansatz_eqs_str)
