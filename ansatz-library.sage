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
# Selected polynomials are NOT forced nonzero up front: differential Thomas
# splits on the relevant initials/separants itself, so each disjoint cell
# carries its own inequations and the degenerate loci are reported, not
# excluded.
#
# build_problem(pde_name, ansatz) returns everything the solver needs:
#   dict(R, rk, coords, roots, jets, params, ansatz_eqs, pconst,
#        pdes, pde_params)
# where `pdes` is a LIST of differential polynomials (a PDE *system*; one
# element for the single-PDE problems) and `pde_params` the PDE's own
# constants (E for the Schroedinger problems, rho/mu for Navier-Stokes).
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
    if pde_name in ('navier-stokes', 'navier-stokes-nd'):
        return ['x', 'y', 'z', 't'], []
    raise ValueError("unknown pde %r" % pde_name)


# ==========================================================================
# the ansatz listing
# ==========================================================================
# Entries are in NUMERICAL order by ansatz number.  Two return shapes:
#   Zeta / single-ODE-function family (5, 5.1, 5.2, 5.3, 8, 9, 19) --
#     dict(order, V, ODE, params, extra) where
#       V     : inner-variable polynomial string (defines jet `v`)
#       ODE   : the ODE as a differential polynomial in Psi/DPsi/DDPsi and v
#       params: ALL constant parameters introduced (v-coeffs then ODE-coeffs)
#       extra : any additional algebraic relations (e.g. algebraic extensions)
#     build_problem assembles the chain rules for these from `order`.
#   every other entry (product 1/1.1/2/3/18/20/20.1, rational 6/7,
#     product-of-two 10, nested 12, algext 13, coeff-ring 14/15/16, log 17/17.1)
#     -- dict(kind, jets_dep, equations, params, v_params, amp_params); its
#     differential-polynomial equations are listed inline.
#
# EVERY entry carries a comment block in the ansatz-10 format: a banner naming
# the template and the closed form, a prose paragraph placing it against its
# neighbours, an indented display of the defining equations, closing prose for
# the design choices and references, and a "Code jet names:" line mapping the
# math notation to the identifiers used below.  The displays are written in the
# HELIUM coordinates (R₁, R₂, R₁₂); hydrogen substitutes (x, y, z) plus the
# algebraic root r, which adds one term to every coordinate trial polynomial.

def ansatz_spec(ansatz, coords, roots):
    gens = coords + [rn for rn, _ in roots]
    rset = tuple(rn for rn, _ in roots)

    # ----- KATO-CUSP EXPONENTIAL template: Ψ = A·exp(B) ----------------------
    # A polynomial amplitude times the exponential of a linear form — the minimal
    # cusp-carrying shape, and the first of the PRODUCT entries (1, 1.1, 2, 3, 18,
    # 20.1).  In a product template Ψ itself is a leader, defined by Ψ − A·Φ = 0;
    # unlike the Zeta family the chain rule then carries the A-factor, since
    # differentiating Ψ − A·Φ gives Ψ[c] = A[c]·Φ + A·Φ[c] and the engine forms
    # A[c] itself (A is built from the independent coordinates).  exp is its own
    # derivative, so Φ = exp(B) needs no derivative jet of its own: the chain rule
    # closes on Φ, and the whole ansatz costs 3 dependent jets.
    #
    #   amplitude:     A = a₀ + a₁·R₁ + a₂·R₂ + a₃·R₁₂
    #   exponent:      B =      b₀·R₁ + b₁·R₂ + b₂·R₁₂        (no constant term —
    #                                                          it would only
    #                                                          rescale A)
    #   the product:   Ψ − A·Φ
    #   chain rule:    Φ[c] − Φ·B[c]                          (Φ′ = Φ)
    #
    # B's coefficients are the v_params, so B ≡ 0 — the exponential collapsing to
    # 1 and leaving a bare polynomial Ψ = A — is the degenerate case, which Thomas
    # splits off as its own cell rather than excluding.  A ≡ 0 (the amp_params) is
    # the trivial Ψ ≡ 0 face.
    # Code jet names: Φ→Phi.

    if ansatz == 1:
        ap, A = trial('a', gens, 1, roots=rset)                    # A(coords)
        bp, B = trial('b', gens, 1, constant=False, roots=rset)    # inner B
        eqs = (['Psi - (%s)*Phi' % A]
               + ['Phi[%s] - Phi*B[%s]' % (c, c) for c in coords]
               + ['B - (%s)' % B])
        return dict(kind='product', jets_dep=['Psi', 'Phi', 'B'],
                    equations=eqs, params=ap + bp, v_params=bp, amp_params=ap)

    # ----- QUADRATIC-AMPLITUDE CUSP template: Ψ = A(deg 2)·exp(B) ------------
    # Ansatz 1 with a Hylleraas-type quadratic amplitude — no new template, one
    # more degree on A.  The extra freedom lets Ψ carry a polynomial node/bulk
    # factor on top of the cusp, the way hydrogen's 2s = (1 − Zr/2)·e^{−Zr/2}
    # carries its node.  (Report ansatz 1.1.)
    #
    #   amplitude:     A = a₀ + a₁·R₁ + … + a₉·R₁₂²           (all deg ≤ 2)
    #   exponent:      B =      b₀·R₁ + b₁·R₂ + b₂·R₁₂
    #   the product:   Ψ − A·Φ
    #   chain rule:    Φ[c] − Φ·B[c]                          (Φ′ = Φ)
    #
    # Degenerate / trivial faces exactly as ansatz 1: B ≡ 0 kills the exponential,
    # A ≡ 0 kills Ψ.
    # Code jet names: Φ→Phi.

    if ansatz == 1.1:
        ap, A = trial('a', gens, 2, roots=rset)
        bp, B = trial('b', gens, 1, constant=False, roots=rset)
        eqs = (['Psi - (%s)*Phi' % A]
               + ['Phi[%s] - Phi*B[%s]' % (c, c) for c in coords]
               + ['B - (%s)' % B])
        return dict(kind='product', jets_dep=['Psi', 'Phi', 'B'],
                    equations=eqs, params=ap + bp, v_params=bp, amp_params=ap)

    # ----- LOGARITHM template: Ψ = A·log(C) ----------------------------------
    # A polynomial amplitude times the logarithm of a linear form.  Ξ = log(C) has
    # Ξ′ = 1/C, so the chain rule has to be cleared of its denominator; the
    # separant of the cleared relation is C, which Thomas turns into the
    # inequation C ≠ 0 — the branch point of the log — on every regular cell,
    # while reporting the C ≡ 0 locus as a cell of its own.
    #
    #   amplitude:     A = a₀ + a₁·R₁ + a₂·R₂ + a₃·R₁₂
    #   argument:      C =      c₀·R₁ + c₁·R₂ + c₂·R₁₂        (no constant term)
    #   the product:   Ψ − A·Ξ
    #   chain rule:    C·Ξ[c] − C[c]                          (Ξ′ = 1/C, cleared)
    #
    # C's coefficients are the v_params: C ≡ 0 kills the log outright.  This is
    # the MULTIPLICATIVE log; for a log admitted additively see 17/17.1, through
    # the inner variable see 19, and inside an exponent see 20/20.1.
    # Code jet names: Ξ→Xi.

    if ansatz == 2:
        ap, A = trial('a', gens, 1, roots=rset)                    # A(coords)
        cp, C = trial('c', gens, 1, constant=False, roots=rset)    # inner C
        eqs = (['Psi - (%s)*Xi' % A]
               + ['C*Xi[%s] - C[%s]' % (c, c) for c in coords]
               + ['C - (%s)' % C])
        return dict(kind='product', jets_dep=['Psi', 'Xi', 'C'],
                    equations=eqs, params=ap + cp, v_params=cp, amp_params=ap)

    # ----- COORDINATE-COEFFICIENT ODE template: Ψ = A·X(B) -------------------
    # helium.sage's "weird second-order mess": a polynomial amplitude times an
    # unknown X obeying a 2nd-order ODE whose coefficients are polynomials in the
    # COORDINATES, not in the inner variable.  That is the one structural break
    # from the Zeta family (5, 5.1–5.3), where the coefficients are polynomials in
    # v — here the relation is not an ODE in B at all, but a linear relation among
    # the jets with coordinate-polynomial coefficients, and an inhomogeneous one
    # at that (the pG term, which no other entry in the library carries).
    #
    #   amplitude:     A  = a₀ + a₁·R₁ + a₂·R₂ + a₃·R₁₂
    #   inner var:     B  =      b₀·R₁ + b₁·R₂ + b₂·R₁₂
    #   coefficients:  pC, pD, pF, pG — one linear coordinate polynomial apiece,
    #                  each with its own constant term (4 params each for helium)
    #   the product:   Ψ − A·X
    #   chain rules:   X[c] − X′·B[c],  X′[c] − X″·B[c]
    #   the relation:  pC·X″ − pD·X′ − pF·X − pG = 0
    #
    # B's coefficients are the v_params (B ≡ 0 collapses X to a constant), A's the
    # amp_params.  Six trial polynomials makes this the widest product entry.
    # Code jet names: X/X′/X″→Chi/DChi/DDChi; pC/pD/pF/pG params are c*/d*/f*/g*.

    if ansatz == 3:
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

    # ----- ZETA template: Ψ = Z(v), one ODE function of one inner variable ---
    # The flagship of the library and the base case of the whole Zeta family: one
    # unknown function Z of one floating linear inner variable v, obeying a
    # 2nd-order linear ODE with polynomial coefficients in v.  Every entry with
    # the `order`/`V`/`ODE` return shape (5, 5.1–5.3, 8, 9, 19) is this template
    # with the degrees or the generator list changed; build_problem assembles the
    # chain rules for them from `order`, so they are not listed in `equations`.
    #
    #   inner variable (floating):  v = v₁·R₁ + v₂·R₂ + v₃·R₁₂   (+ v₄·r for H)
    #   chain rules (Z = Z(v)):     Ψ[c] − Ψ′·v[c],  Ψ′[c] − Ψ″·v[c]
    #   the ODE:      A(v)·Ψ″ + B(v)·Ψ′ + C(v)·Ψ = 0            (deg ode_deg)
    #
    # The decimal selects (maxdeg_v, ode_deg) — the degree of the inner variable
    # in the coordinates and of the ODE coefficients in v:
    #   5   = (1,1)   5.1 = (2,1)   5.2 = (1,2)   5.3 = (2,2)
    # 5.1 lets v be a conic (R₁², R₁·R₂, …) rather than a plane, which is the
    # lowest-degree shape a genuinely new coordinate can take.  5.2 raises the ODE
    # instead: degree 2 in the leading coefficient is what a 2nd-order linear ODE
    # needs to reach two regular singular points and a third at infinity, so the
    # hypergeometric equation — and with it the Legendre / Gegenbauer / Jacobi
    # factors — first becomes expressible there.  5.3 is the union of both, and
    # the largest member of the single-Zeta family; helium.sage numbered it and
    # its 10 identically, a coincidence that no longer holds here because ansatz
    # 10 is repurposed as the product-of-two-ODE-functions family.
    #
    # v FLOATS — it is not hard-coded to a radius — so for hydrogen the search
    # returns both the radial solutions (v ∝ r, the 1s exponential) and the E = 0
    # J₀ Bessel family on a plane, as separate primes of the projected ideal.
    # Ansatz 10 subsumes ansatz 5 as its G ≡ const face.
    #
    # TWO NAMING CONVENTIONS live under this branch, and they are load-bearing:
    # ansatz 5 predates the variants, numbers its inner-variable parameters from
    # v₁, calls the ODE coefficients a/b/c and writes the ODE with plus signs —
    # which is how the ideals are labelled in the paper, so it must not drift.
    # 5.1–5.3 number from v₀, use d/m/n and minus signs.  Both are cosmetic (the
    # unknown coefficients absorb the sign), but the output is read by humans and
    # transcribed, so the merge below preserves each verbatim rather than
    # normalizing one into the other.
    # Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi.

    if int(ansatz) == 5:
        maxdeg_v, ode_deg = {5: (1, 1), 5.1: (2, 1),
                             5.2: (1, 2), 5.3: (2, 2)}[ansatz]
        # (coefficient bases, v-numbering origin, ODE sign) -- see above.
        bases, start, sgn = (('a', 'b', 'c'), 1, '+') if ansatz == 5 \
                       else (('d', 'm', 'n'), 0, '-')
        vp, V = trial('v', gens, maxdeg_v, constant=False, start=start,
                      roots=rset)
        ap, A = trial(bases[0], ['v'], ode_deg)
        bp, B = trial(bases[1], ['v'], ode_deg)
        cp, C = trial(bases[2], ['v'], ode_deg)
        ODE = '(%s)*DDPsi %s (%s)*DPsi %s (%s)*Psi' % (A, sgn, B, sgn, C)
        return dict(order=2, V=V, ODE=ODE, params=vp + ap + bp + cp, extra=[])

    # ----- RATIONAL-ARGUMENT template: Ψ = Z(w), w = B/C ---------------------
    # The Zeta family with a RATIONAL inner argument.  Rather than carry a
    # denominator, w enters as a jet with the cleared defining relation C·w − B;
    # its separant is C, so C ≠ 0 becomes a Thomas inequation on every regular
    # cell and the C ≡ 0 locus is reported as a cell rather than assumed away.
    # Ansatz 6 takes B, C and the ODE coefficients all linear; ansatz 7 takes all
    # three quadratic — the two share this one code path.
    #
    #   numerator, denominator:  B, C — coordinate polynomials of degree d, both
    #                            with a constant term (d = 1 for 6, d = 2 for 7)
    #   inner variable:          C·w − B = 0                   (i.e. w = B/C)
    #   chain rules:             Ψ[c] − Ψ′·w[c],  Ψ′[c] − Ψ″·w[c]
    #   the ODE:  D(w)·Ψ″ − M(w)·Ψ′ − N(w)·Ψ = 0               (each deg d in w)
    #
    # B/C is invariant under (B,C) → (λB, λC), an extra scaling dimension we let
    # Thomas carry rather than gauge-fix, so every prime arrives with that one
    # spurious degree of freedom.  B's coefficients are the v_params: B ≡ 0 gives
    # w = 0, the degenerate cell.
    # Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi.

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

    # ----- FIRST-ORDER ZETA template: Ψ = Z(v), M·Z′ = N·Z -------------------
    # The 1st-order face of the Zeta family: one chain rule instead of two, one
    # ODE of order 1 with coefficients linear in v.  Its solutions are the
    # exponentials of a rational integral, so it is the cheapest ansatz that can
    # still express a Kato cusp.  Ansatz 5 finds the same solutions on its
    # a₀ = a₁ = 0 face, but carries the Ψ″ jet the whole way to get there — so
    # this is the one to reach for when 5 is too heavy to close.
    #
    #   inner variable:  v = v₀·R₁ + v₁·R₂ + v₂·R₁₂
    #   chain rule:      Ψ[c] − Ψ′·v[c]                       (order 1: one rule)
    #   the ODE:         (m₀ + m₁·v)·Ψ′ − (n₀ + n₁·v)·Ψ = 0
    #
    # Code jet names: Z/Z′→Psi/DPsi.

    if ansatz == 8:
        vp, V = trial('v', gens, 1, constant=False, roots=rset)
        mp, M = trial('m', ['v'], 1)
        np_, N = trial('n', ['v'], 1)
        ODE = '(%s)*DPsi - (%s)*Psi' % (M, N)
        return dict(order=1, V=V, ODE=ODE, params=vp + mp + np_, extra=[])

    # ----- CONSTANT-RATE ZETA template: Ψ = Z(v), Z′ = n₀·Z ------------------
    # Ansatz 8 with the rate held constant.  The ODE integrates in closed form to
    # Ψ = const·exp(n₀·v), so this is ansatz 1 with a constant amplitude, reached
    # through the Zeta machinery instead of the product machinery — and that
    # overlap is the point: the two templates must return the same variety, which
    # makes this the standard cross-check on a new PDE, a new ranking, or a change
    # to the reduction pipeline.  4 parameters for helium, the smallest in the
    # library.
    #
    #   inner variable:  v = v₀·R₁ + v₁·R₂ + v₂·R₁₂
    #   chain rule:      Ψ[c] − Ψ′·v[c]
    #   the ODE:         Ψ′ − n₀·Ψ = 0
    #
    # Code jet names: Z/Z′→Psi/DPsi.

    if ansatz == 9:
        vp, V = trial('v', gens, 1, constant=False, roots=rset)
        np_, N = trial('n', ['v'], 0)          # N = n0
        ODE = 'DPsi - (%s)*Psi' % N
        return dict(order=1, V=V, ODE=ODE, params=vp + np_, extra=[])

    # ----- PRODUCT-OF-TWO-ODE-FUNCTIONS template: Ψ = F(ξ)·G(η) --------------
    # A product of two unknown functions, each of its OWN floating linear inner
    # variable, each obeying its own 2nd-order ODE with coefficients linear in
    # that inner variable.  The minimal product analogue of ansatz 5 — and
    # ansatz 5 is its G≡const face, so this SUBSUMES ansatz 5: a hydrogen run
    # must re-find ansatz 5's radial + E=0 J₀-Bessel varieties as the single-
    # factor faces, plus — if the two-factor machinery closes — the E≠0
    # parabolic product on top.
    #
    #   inner variables (floating):  ξ = p₀·R₁ + p₁·R₂ + p₂·R₁₂   (+ a radial
    #                                η = q₀·R₁ + q₁·R₂ + q₂·R₁₂    term r for H)
    #   chain rules (F=F(ξ), G=G(η)):  F[c] − F′·ξ[c],  F′[c] − F″·ξ[c]
    #                                  G[c] − G′·η[c],  G′[c] − G″·η[c]
    #   product:                       Ψ − F·G
    #   the two ODEs:
    #       (a₀ + a₁·ξ)·F″ + (b₀ + b₁·ξ)·F′ + (c₀ + c₁·ξ)·F  = 0
    #       (d₀ + d₁·η)·G″ + (e₀ + e₁·η)·G′ + (f₀ + f₁·η)·G  = 0
    #
    # The inner variables FLOAT (not hard-coded): for hydrogen this discovers the
    # parabolic coordinates ξ=r+z, η=r−z of Landau–Lifshitz vol.3 §37 (up to the
    # arbitrary field axis); for helium it discovers a Hylleraas-type product
    # F(linear)·G(linear) — no separation is known there, so this is a genuine
    # coordinate-discovering search.  18 params + E, 9 dependent jets: expect the
    # heavy band.  Kept fully symmetric — NO leading-coefficient gauge-fix, so the
    # a₁=0 / d₁=0 branches (incl. the pure-Bessel face) stay reachable.
    # Code jet names: F″/F′→DDF/DF, G″/G′→DDG/DG, ξ/η→xi/eta.

    if ansatz == 10:
        pp, XI  = trial('p', gens, 1, constant=False, roots=rset)
        qq, ETA = trial('q', gens, 1, constant=False, roots=rset)
        eqs = (['xi - (%s)' % XI, 'eta - (%s)' % ETA]
               + ['F[%s] - DF*xi[%s]'    % (c, c) for c in coords]
               + ['DF[%s] - DDF*xi[%s]'  % (c, c) for c in coords]
               + ['G[%s] - DG*eta[%s]'   % (c, c) for c in coords]
               + ['DG[%s] - DDG*eta[%s]' % (c, c) for c in coords]
               + ['Psi - F*G']
               + ['(a0 + a1*xi)*DDF + (b0 + b1*xi)*DF + (c0 + c1*xi)*F']
               + ['(d0 + d1*eta)*DDG + (e0 + e1*eta)*DG + (f0 + f1*eta)*G'])
        odep = ['a0', 'a1', 'b0', 'b1', 'c0', 'c1', 'd0', 'd1', 'e0', 'e1', 'f0', 'f1']
        return dict(kind='product2',
                    jets_dep=['Psi', 'DDF', 'DF', 'F', 'DDG', 'DG', 'G', 'xi', 'eta'],
                    equations=eqs, params=pp + qq + odep,
                    v_params=pp + qq, amp_params=[])

    # ----- NESTED-ODE template: Ψ = Z(V), V = V(coords, Θ), Θ = Θ(U) ---------
    # Two coupled ODE functions, one inside the other.  The inner Θ solves its own
    # 2nd-order ODE in an inner-inner variable U; the OUTER inner variable V is a
    # polynomial in the coordinates AND Θ, so Ψ = Z(V) is a function of a
    # function.  Reduction of the PDE cascades down both towers:
    # Ψ[c] → Ψ′·V[c],  V[c] → (coords) + v_Θ·Θ[c],  Θ[c] → Θ′·U[c].
    #
    #   inner-inner var:  U = u₀·R₁ + u₁·R₂ + u₂·R₁₂          (deg maxdeg_u)
    #   inner chain:      Θ[c] − Θ′·U[c],  Θ′[c] − Θ″·U[c]
    #   inner ODE:        A(U)·Θ″ − B(U)·Θ′ − C(U)·Θ = 0      (deg ode_u in U)
    #   outer var:        V = v₀·R₁ + v₁·R₂ + v₂·R₁₂ + v₃·Θ   (deg maxdeg_v)
    #   outer chain:      Ψ[c] − Ψ′·V[c],  Ψ′[c] − Ψ″·V[c]
    #   outer ODE:        D(V)·Ψ″ − M(V)·Ψ′ − N(V)·Ψ = 0      (deg ode_v in V)
    #
    # The decimal selects (maxdeg_v, maxdeg_u, ode_v, ode_u):
    #   12   = (1,1,1,1)   12.1 = (2,1,1,1)   12.2 = (1,2,1,1)
    #   12.3 = (1,1,2,1)   12.4 = (1,1,1,2)
    # 8 dependent jets even at the base degrees, second only to ansatz 10's 9 —
    # two full ODE towers is inherently what this costs.
    # Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi, Θ/Θ′/Θ″→Theta/DTheta/DDTheta.

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

    # ----- ALGEBRAIC-EXTENSION template: Ψ = Z(V), ODE coeffs over Q(V)[g] ---
    # (helium.sage calls this one 13.)  The Zeta family with the ODE coefficients
    # living in a quadratic algebraic extension: g is a differential indeterminate
    # pinned by an order-0 minimal polynomial whose own coefficients are
    # polynomials in V, and the ODE coefficients D, M, N are polynomials in V AND
    # g.  g's derivatives follow from prolonging the minimal polynomial, whose
    # separant is 2A·g + B — so the branch locus of the extension arrives as a
    # Thomas inequation rather than an assumption.  Ansatz 11 (γ) needs exactly
    # this template.
    #
    #   inner variable:  V = v₀·R₁ + v₁·R₂ + v₂·R₁₂           (deg maxdeg_v)
    #   minimal poly:    A(V)·g² + B(V)·g + C(V) = 0          (deg alg_deg in V)
    #   chain rules:     Ψ[c] − Ψ′·V[c],  Ψ′[c] − Ψ″·V[c]
    #   the ODE:  D(V,g)·Ψ″ − M(V,g)·Ψ′ − N(V,g)·Ψ = 0        (deg ode_deg)
    #
    # The decimal selects (maxdeg_v, ode_deg, alg_deg):
    #   13   = (1,1,1)   13.1 = (2,1,1)   13.2 = (1,2,1)   13.3 = (1,1,2)
    #   13.4 = (2,2,1)   13.5 = (2,1,2)   13.6 = (2,2,2)
    # NOTE at ode_deg = 2 the ODE coefficients carry a g² monomial; on the A ≠ 0
    # cells the minimal polynomial reduces it back to degree 1 in g, so that
    # freedom is largely redundant there.  Contrast ansatz 16, where the algebraic
    # element sits in the BASE rather than in the coefficient ring.
    # Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi.

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

    # ----- EXPONENTIAL COEFFICIENT-RING template: coeffs over Q(v)[exp(e₀v)] -
    # The first of the three COEFFICIENT-RING entries (14, 15, 16): Ψ = Z(v) with
    # a 2nd-order ODE whose coefficients D, M, N live in an EXTENSION of the
    # coefficient ring.  Ansatz 13 above is the algebraic case of the same idea;
    # 14 here is the exponential one.  An exponential t = exp(e₀·v) sits in the
    # coefficient ring — it multiplies D, M, N, NOT Ψ — with the same defining
    # relation shape as ansatz 1's Φ.  That distinction is the whole content of
    # the entry: in ansatz 1 the exponential is a factor of the solution, here it
    # is a coefficient, so the solution is a Ψ whose ODE has exponential
    # coefficients rather than a Ψ that is itself an exponential.
    #
    #   inner variable:   v = v₀·R₁ + v₁·R₂ + v₂·R₁₂
    #   the exponential:  t[c] − e₀·t·v[c]                    (t = exp(e₀·v))
    #   chain rules:      Ψ[c] − Ψ′·v[c],  Ψ′[c] − Ψ″·v[c]
    #   the ODE:  D(v,t)·Ψ″ − M(v,t)·Ψ′ − N(v,t)·Ψ = 0        (deg 1 in v and t)
    #
    # ONE rate parameter, e₀.  From NewSol.tex's ansatz collection — the section
    # removed in Papers commit fc28cfca (Oct 2024) — read off the schematic
    # diagrams and implemented as the base case only (degree 1 where possible);
    # getting these to parse and decompose at all is the goal.
    # Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi.

    if int(ansatz) == 14:
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

    # ----- HOLONOMIC COEFFICIENT-RING template: coeffs over Q(v)[t, t′] ------
    # NewSol.tex's 2nd-order HOLONOMIC coefficient-ring case.  A holonomic element
    # t (jets t, t′, t″) solves its OWN 2nd-order ODE in v, and Ψ's ODE
    # coefficients are degree-1 polynomials in v, t and t′.  This is the general
    # case of which ansatz 14 is the 1st-order specialization: any Bessel / Airy /
    # hypergeometric factor can now appear in the coefficient ring.
    #
    #   inner variable:   v = v₀·R₁ + v₁·R₂ + v₂·R₁₂
    #   t's chain rules:  t[c] − t′·v[c],  t′[c] − t″·v[c]
    #   t's own ODE:      P(v)·t″ − Q(v)·t′ − S(v)·t = 0      (deg 1 in v)
    #   chain rules:      Ψ[c] − Ψ′·v[c],  Ψ′[c] − Ψ″·v[c]
    #   the ODE:  D(v,t,t′)·Ψ″ − M(v,t,t′)·Ψ′ − N(v,t,t′)·Ψ = 0
    #
    # Heavy — 7 dependent jets, two full towers of them.  That is inherent to a
    # holonomic coefficient ring, not an artifact of the encoding.  Base case only
    # (degree 1 everywhere), same provenance as ansatz 14.
    # Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi, t′/t″→Dt/DDt.

    if int(ansatz) == 15:
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

    # ----- ALGEBRAIC-BASE template: Ψ = Z(v), v over the coordinates AND g ---
    # NewSol.tex / helium.sage's coded ansatz 16: an algebraic root nested BELOW
    # the extension.  The novelty is in the BASE, not the coefficient ring — g is
    # a square root, g² − RAD = 0 with RAD an unknown coordinate polynomial, and
    # the inner variable v ranges over the coordinates AND g.  At rad_deg = 2
    # (16.3, 16.5, 16.6) g is the root of a QUADRATIC: a genuinely new coordinate
    # of the same shape as hydrogen's r = √(x²+y²+z²), which is the ingredient
    # that let ansatz 5 find the J₀ Bessel solution.  The difference is that here
    # the radicand is searched for rather than supplied.
    #
    #   the root:        g² − RAD = 0,  RAD = k₀ + k₁·R₁ + …  (deg rad_deg)
    #   inner variable:  v = v₀·R₁ + v₁·R₂ + v₂·R₁₂ + v₃·g    (deg v_deg)
    #   chain rules:     Ψ[c] − Ψ′·v[c],  Ψ′[c] − Ψ″·v[c]
    #   the ODE:  A(v)·Ψ″ + B(v)·Ψ′ + C(v)·Ψ = 0              (deg ode_deg in v)
    #
    # The decimal selects (rad_deg, v_deg, ode_deg), matching the NewSol.tex
    # Ansatz 16 table:
    #   16   = (1,1,1)   16.1 = (1,2,1)   16.2 = (1,1,2)   16.3 = (2,1,1)
    #   16.4 = (1,2,2)   16.5 = (2,1,2)   16.6 = (2,2,2)
    # NOTE g is NOT in the trial's `roots` set, so at v_deg = 2 the inner variable
    # carries a g² monomial that g² − RAD immediately reduces back to the
    # radicand — one redundant parameter on 16.1, 16.4 and 16.6.
    # Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi; RAD's params are k₀….

    if int(ansatz) == 16:
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

    # ----- LOG-HYPERRADIUS LINEAR template: Ψ = A + B·L ----------------------
    # The first of the LOG-HYPERRADIUS entries (17, 17.1, 18, 19), which carry the
    # Bartlett–Fock log-non-analyticity L = log(R₁² + R₂²) = 2·log s, s being the
    # triple-coalescence hyperradius.  L enters as a transcendental jet with
    # rational (cleared) derivatives and NO order-0 relation — log is
    # transcendental, so its free order-0 jet is exactly the additive integration
    # constant.  This entry admits it ADDITIVELY and linearly, the minimal form
    # with a genuine log s.  HELIUM-SPECIFIC: _log_relations hard-codes the helium
    # hyperradius gradient, and 17.1 additionally hard-codes (R₁² + R₂²) in the
    # trial form, so a hydrogen run is a build check and nothing more.
    #
    #   the log jet:  (R₁² + R₂²)·L[c] − ∂c(R₁² + R₂²)        (leaders L[c])
    #   17    Ψ − (A + B·L),  A, B deg 1                      (report ansatz 14)
    #   17.1  Ψ − (A + b₀·(R₁² + R₂²)·L),  A deg 2            (report ansatz 14.1)
    #
    # 17.1 pins the log to the exact O(s² log s) Fock slot instead of letting it
    # float.  B ≡ 0 (resp. b₀ ≡ 0) is the v_params face — the log gone; A ≡ B ≡ 0
    # is the trivial Ψ ≡ 0 face.  Both keep the log LINEAR and multiplicative, so
    # neither can represent more than (ln s)¹, which is why both came back
    # clean-negative; ansatz 20 is the exponentiated form that holds the tower.
    # Refs: ~/project/reports/helium-new-ansatze.md.

    if int(ansatz) == 17:
        if ansatz == 17:
            ap, A = trial('a', gens, 1, roots=rset)
            bp, B = trial('b', gens, 1, roots=rset)
            eqs = (['Psi - ((%s) + (%s)*L)' % (A, B)]   # leader Psi
                   + _log_relations(coords))            # leaders L[c]
            return dict(kind='loglin', jets_dep=['Psi', 'L'],
                        equations=eqs, params=ap + bp,
                        v_params=bp,            # B == 0  -> log gone -> DEGENERATE
                        amp_params=ap + bp)     # A==B==0 -> Psi == 0 -> TRIVIAL
        else:                                           # ansatz 17.1
            ap, A = trial('a', gens, 2, roots=rset)
            eqs = (['Psi - ((%s) + b0*(R1^2 + R2^2)*L)' % A]
                   + _log_relations(coords))
            return dict(kind='loglin', jets_dep=['Psi', 'L'],
                        equations=eqs, params=ap + ['b0'],
                        v_params=['b0'], amp_params=ap + ['b0'])

    # ----- KATO × FOCK template: Ψ = exp(B)·(A₀ + A₁·L) ----------------------
    # The cusp exponential times the first Fock log — ansatz 1's product template
    # extended by the log jet of ansatz 17.  These are the two singular structures
    # of the helium wavefunction at the triple-coalescence point, carried
    # simultaneously and multiplicatively.  HELIUM-SPECIFIC, for the same
    # _log_relations reason as 17.  (Report ansatz 15.)
    #
    #   exponent:      B  =      b₀·R₁ + b₁·R₂ + b₂·R₁₂
    #   amplitudes:    A₀ = a₀ + a₁·R₁ + a₂·R₂ + a₃·R₁₂
    #                  A₁ = h₀ + h₁·R₁ + h₂·R₂ + h₃·R₁₂
    #   the product:   Ψ − Φ·(A₀ + A₁·L)
    #   chain rule:    Φ[c] − Φ·B[c]                          (Φ = exp(B))
    #   the log jet:   (R₁² + R₂²)·L[c] − ∂c(R₁² + R₂²)
    #
    # A₁ ≡ 0 (the v_params) drops the log and lands back in ansatz-1 territory;
    # A₀ ≡ A₁ ≡ 0 is the trivial Ψ ≡ 0 face.
    # Code jet names: Φ→Phi; A₁'s parameters are named h₀…, not a₄….

    if ansatz == 18:
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

    # ----- LOG-EXTENDED ZETA template: Ψ = Z(v), v over coordinates AND L ----
    # The Zeta family with the log admitted through the INNER VARIABLE rather than
    # through an amplitude: v ranges over the coordinates and L, so Ψ is an
    # unknown function of a coordinate-plus-log combination.  Rides the Zeta
    # branch of build_problem via the `extra_jets` hook, with the log relations
    # passed through the family's `extra` slot and L ranked just below v.
    #
    #   inner variable:  v = v₀·R₁ + v₁·R₂ + v₂·R₁₂ + v₃·L
    #   the log jet:     (R₁² + R₂²)·L[c] − ∂c(R₁² + R₂²)
    #   chain rules:     Ψ[c] − Ψ′·v[c],  Ψ′[c] − Ψ″·v[c]
    #   the ODE:  (a₀ + a₁·v)·Ψ″ + (b₀ + b₁·v)·Ψ′ + (c₀ + c₁·v)·Ψ = 0
    #
    # Unlike 17/17.1/18 the log is not confined to degree 1 in the answer: Z is an
    # unknown function, so Z(… + v₃·L) can carry powers of the log — but only in
    # that single combination v, which is the restriction ansatz 20 removes.
    # (Report ansatz 16.)  Code jet names: Z/Z′/Z″→Psi/DPsi/DDPsi.

    if ansatz == 19:
        vp, V = trial('v', gens + ['L'], 1, constant=False, roots=rset)
        ap, A = trial('a', ['v'], 1)
        bp, B = trial('b', ['v'], 1)
        cp, C = trial('c', ['v'], 1)
        ODE = '(%s)*DDPsi + (%s)*DPsi + (%s)*Psi' % (A, B, C)
        return dict(order=2, V=V, ODE=ODE, params=vp + ap + bp + cp,
                    extra=_log_relations(coords), extra_jets=['L'])

    # ----- EXPONENTIATED FOCK template: Ψ = exp(B + C·L) ---------------------
    # Myers–Umrigar–Sethna–Morgan 1991, sec. IV: the log goes INSIDE the exponent.
    # This is the first library entry to put it there, and the reason it matters
    # is that exp(C·L) = (R₁² + R₂²)^C is a variable-exponent (xˣ-type) factor
    # whose Taylor series carries the UNBOUNDED (ln s)^p Fock tower with a finite
    # parameter set.  17/17.1/18 keep the log linear and multiplicative and so can
    # only ever represent (ln s)¹ — which is why they came back clean-negative:
    # too holonomic to hold the tower.
    #
    #   log-free exponent:  B = b₀·R₁ + … + b₈·R₁₂²  (deg 2, no constant) — the
    #       Kato cusp ψ₁,₀ = −Z(R₁ + R₂) + R₁₂/2 plus the O(s²) pieces
    #       ψ₂,₀ − ½·ψ₁,₀².
    #   log coefficient:    C = c₀ + c₁·R₁ + … + c₉·R₁₂²  (deg 2) — the Fock r²
    #       log slot; the exact leading term is
    #       ψ₂,₁ ∝ (π−2)/(3π)·Z·(R₁² + R₂² − R₁₂²) = (π−2)/(3π)·Z·Y₂,₁, so C
    #       must reach degree 2.
    #   the exponential:    Ψ[c] − Ψ·(B[c] + C[c]·L + C·L[c])
    #   the log jet:        (R₁² + R₂²)·L[c] − ∂c(R₁² + R₂²)
    #
    # (R₁, R₂, R₁₂) ARE the KS-rationalized coordinates — the hyperspherical
    # √(1 − sin a cos t), √(1 + sin a) irrationalities of the angular Fock
    # coefficients are polynomial here, so L is the only transcendental jet.
    # C ≡ 0 (the v_params) drops the log and leaves Kato-only: DEGENERATE.  exp is
    # never zero, so there is no Ψ ≡ 0 (TRIVIAL) mode and amp_params is empty.
    # Refs: MUSM 1991 eq (23); Liverts 2022 eq (8); Fournais et al. 2004/2009;
    # ~/project/reports/helium-new-ansatze.md.

    if ansatz == 20:
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

    # ----- EXPONENTIATED FOCK + AMPLITUDE: Ψ = A·exp(B + C·L) ----------------
    # Ansatz 20 with a Hylleraas-style polynomial prefactor (cf. 1 → 1.1).  The
    # amplitude lets Ψ carry a polynomial node/bulk factor on top of the exp-Fock
    # singular structure, the way hydrogen's 2s = (1 − Zr/2)·e^{−Zr/2} carries its
    # node.  Introducing A also brings the TRIVIAL mode back, which ansatz 20 does
    # not have: A ≡ 0 gives Ψ ≡ 0.
    #
    #   amplitude:          A = a₀ + a₁·R₁ + a₂·R₂ + a₃·R₁₂    (deg 1)
    #   log-free exponent:  B = b₀·R₁ + … + b₈·R₁₂²            (deg 2, no const)
    #   log coefficient:    C = c₀ + c₁·R₁ + … + c₉·R₁₂²       (deg 2)
    #   the product:        Ψ − A·Φ
    #   the exponential:    Φ[c] − Φ·(B[c] + C[c]·L + C·L[c])
    #   the log jet:        (R₁² + R₂²)·L[c] − ∂c(R₁² + R₂²)
    #
    # C ≡ 0 (the v_params) drops the log: DEGENERATE.  Refs as ansatz 20.
    # Code jet names: Φ→Phi.

    if ansatz == 20.1:
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

    # ----- NONLINEAR SECOND-DEGREE SYSTEM template: (u,v,w,p) over Ψ = Z(s) --
    # Figure 8 / equation (25) of NewMethod.tex — hence the number — the
    # nonlinear, second-order, second-degree ansatz whose four top-level
    # polynomials u, v, w, p were introduced precisely to serve as the
    # dependent variables of a SYSTEM of PDEs (incompressible Navier–Stokes
    # being the target pairing, sec:NavierStokes).  One unknown profile Z of
    # one floating linear inner variable — called s here, not v, so as not to
    # collide with the second velocity component — obeying a 2nd-order ODE
    # that carries all six second-degree jet monomials alongside the three
    # linear ones; each field is a linear background in the coordinates plus
    # the ODE element with its own amplitude.
    #
    #   chain rules:     Ψ[c] − Ψ′·s[c],  Ψ′[c] − Ψ″·s[c]      (c = x, y, z, t)
    #   the ODE:   (a₀+a₁s)Ψ″² + (b₀+b₁s)Ψ″Ψ′ + (c₀+c₁s)Ψ″Ψ + (d₀+d₁s)Ψ′²
    #              + (e₀+e₁s)Ψ′Ψ + (f₀+f₁s)Ψ² + (g₀+g₁s)Ψ″ + (h₀+h₁s)Ψ′
    #              + (i₀+i₁s)Ψ = 0
    #   inner variable:  s = s₁x + s₂y + s₃z + s₄t
    #   the fields:      u = u₁x + u₂y + u₃z + u₄t + u₅Ψ       (v, w, p alike)
    #
    # The decimal selects (leader_deg, background):
    #   25   = (2, yes)  42 params      25.1 = (1, yes)  40 params
    #   25.2 = (2, no)   30 params      25.3 = (1, no)   28 params
    # leader_deg = 2 includes the (a₀+a₁s)Ψ″² monomial; leader_deg = 1 omits
    # it — the cell where the paper says the viscous solutions live, since a
    # momentum residual of degree 1 in Ψ″ cannot be pseudo-divided by an ODE
    # of degree 2 in Ψ″ (the "viscous block" of sec:NavierStokes).  Everything
    # else in the ODE is unchanged: the other five second-degree monomials
    # stay in both cases.  background = no strips u, v, w to their Ψ term
    # alone (u − u₅Ψ), but p KEEPS its full linear part, so the run can
    # confirm p₁ = p₂ = p₃ = 0 rather than assume it.
    #
    # u₅, v₅, w₅ are the v_params — the amplitude vector a of the paper.
    # a ≡ 0 removes the ODE element from the velocity field entirely, which is
    # the DEGENERATE case.  amp_params is EMPTY: Ψ ≡ 0 leaves the linear
    # background, which is a real (if dull) solution of the system — one of
    # the linear strain-and-rotation flows — not a collapse to a trivial
    # solution, so no prime is tagged TRIVIAL on that account.
    # Refs: NewMethod.tex Figure 8 (nonlinear ansatz figure) and eq. (25)
    # [label `ansatz 18`], sec:NavierStokes; oracle ns-reduction-check.py.
    # Code jet names: Ψ/Ψ′/Ψ″→Psi/DPsi/DDPsi; the ODE's independent variable
    # is the jet s.
    #
    # NORMALIZATION LADDER 25.31–25.34 — a cheapest-last ladder of
    # progressively more normalized variants of the 25.3 corner
    # (leader_deg = 1, background = no, 28 params), selected by a third
    # dispatch field `norm`.  Levels are CUMULATIVE (level 3 includes
    # everything in levels 1 and 2); each rung is a strictly smaller
    # parameter space than the one before.  Normalization removes
    # parameters, never equations: every rung keeps all 14 ansatz equations.
    #
    # norm ≥ 1 (25.31, 25 params) — ROTATION NORMAL FORM.  Incompressible
    # Navier–Stokes is invariant under rotations of (x,y,z), with (u,v,w)
    # rotating as a vector.  Two of the three rotations align the wave
    # vector with the x-axis: s₂ = s₃ = 0, so s = s₁x + s₄t.  The residual
    # rotation about that axis puts the amplitude vector in the x–y plane:
    # w₅ = 0, so w ≡ 0.  The exclusion is NOT vacuous over an algebraically
    # closed field: SO(3,ℂ) orbits on ℂ³ are classified by
    # σ = s₁²+s₂²+s₃², and a nonzero isotropic vector (σ = 0, s ≠ 0) cannot
    # be rotated to the form (s₁,0,0) — that would force s₁ = 0.  Task 471
    # found exactly such a component (one of 25.2's eight primes contains
    # s₁²+s₂²+s₃²).  Over ℝ it is vacuous; over ℚ̄ it is a real component
    # that rungs 1–4 drop — the "punctured null cone" `excludes` entry.
    #
    # norm ≥ 2 (25.32, 22 params) — SCALING CHARTS.  Three commuting
    # scalings act on the parameter space and must each be pinned by one
    # normalization: κ (multiply the whole ODE by a constant), λ (Ψ → λΨ),
    # ν (s → νs).  Pins: s₁ = 1 (s = x + s₄t), v₅ = 1 (v = Ψ), g₀ = 1
    # (Ψ″ term becomes (1 + g₁s)Ψ″).  Their weights (κ,λ,ν) — the exponent
    # triple in p ⟼ κᵃλᵇνᶜ·p — are (0,0,1), (0,−1,0) and (1,−1,2), whose
    # 3×3 determinant is 1, so together they rigidify the torus with no
    # residual finite stabilizer.  Do NOT normalize u₅ instead of v₅:
    # continuity is (a·s)Ψ′ = u₅s₁Ψ′ once s₂ = s₃ = 0, so in the s₁ ≠ 0
    # chart it forces u₅ = 0 (modulo the degenerate Ψ′ ≡ 0); normalizing
    # u₅ = 1 would silently select the s₁ = 0 branch instead.  v₅ is the
    # transverse amplitude — the Stokes-layer direction.  v_params is EMPTY
    # from this rung on (see the code comment at the v_params assignment).
    #
    # norm ≥ 3 (25.33, 19 params) — PRESSURE TRIM: p₁ = p₂ = p₃ = 0, so
    # p = p₄t + p₅Ψ.  The one rung that is NOT a symmetry quotient.  Task
    # 471's generic-cell GTZ found p₁,p₂,p₃ in all 8 of 25.2's minimal
    # primes and in 25.3's genuine prime — with background = no every
    # inertial and viscous term carries a factor of Ψ, so p₁ is the sole
    # jet-degree-0 term in momentum-x and must vanish on its own.  But that
    # argument runs through the pseudo-division cofactor
    # (b₀+b₁s)Ψ′ + (c₀+c₁s)Ψ + (g₀+g₁s), so it is a GENERIC-CELL result;
    # special cells where the cofactor vanishes are not covered — hence the
    # `excludes` entry.
    #
    # norm = 4 (25.34, 9 params) — ODE TRUNCATION: drop the Ψ-quadratic
    # terms, leaving the linear ODE (1+g₁s)Ψ″ + (h₀+h₁s)Ψ′ + (i₀+i₁s)Ψ = 0
    # (g₀ already normalized to 1 by level 2).  This is the family the
    # Rayleigh–Stokes shear layer lives in: task 471's genuine prime
    # contains f₀,f₁,i₀,i₁ and 2×2 minors locking (d,e,h) ∝ (b,c,g), so the
    # quadratic block is inert there.  Do NOT additionally impose u₅ = 0
    # even though continuity implies it in this chart — it is a derived
    # consequence, not a normalization; folding it into the ansatz
    # definition would add a second excluded branch (Ψ′ ≡ 0) and make the
    # result harder to state.  Let the decomposition derive it.  Hence 9
    # parameters, not 8.
    #
    # TYPE DISTINCTION (for the paper): rungs 1–3 are CHARTS of ansatz 25 —
    # same equations, restricted parameter space — so their results lift
    # back to 25.3 modulo the recorded excluded locus (the `excludes` key,
    # carried as data and forwarded by build_problem).  Rung 4 is a
    # DIFFERENT ANSATZ FAMILY and its results do NOT lift; its `excludes`
    # is empty for that reason, not because nothing is given up.
    #
    #   25.31 = norm 1  25 params      25.32 = norm 2  22 params
    #   25.33 = norm 3  19 params      25.34 = norm 4   9 params

    if int(ansatz) == 25:
        assert coords == ['x', 'y', 'z', 't'], \
            "ansatz 25 is the Navier-Stokes system ansatz: it requires " \
            "coords ['x','y','z','t'], got %r" % (coords,)
        leader_deg, background, norm = {
            25:    (2, True,  0), 25.1:  (1, True,  0),
            25.2:  (2, False, 0), 25.3:  (1, False, 0),
            25.31: (1, False, 1), 25.32: (1, False, 2),
            25.33: (1, False, 3), 25.34: (1, False, 4),
        }[ansatz]
        if norm == 4:                       # linear ODE (quadratic block gone)
            ode_terms = ['(1 + g1*s)*DDPsi', '(h0 + h1*s)*DPsi',
                         '(i0 + i1*s)*Psi']
            odep = ['g1', 'h0', 'h1', 'i0', 'i1']
        else:
            ode_terms = (['(a0 + a1*s)*DDPsi^2'] if leader_deg == 2 else [])
            ode_terms += ['(b0 + b1*s)*DDPsi*DPsi', '(c0 + c1*s)*DDPsi*Psi',
                          '(d0 + d1*s)*DPsi^2', '(e0 + e1*s)*DPsi*Psi',
                          '(f0 + f1*s)*Psi^2',
                          '(%s + g1*s)*DDPsi' % ('1' if norm >= 2 else 'g0'),
                          '(h0 + h1*s)*DPsi', '(i0 + i1*s)*Psi']
            odep = ((['a0', 'a1'] if leader_deg == 2 else [])
                    + ['b0', 'b1', 'c0', 'c1', 'd0', 'd1', 'e0', 'e1',
                       'f0', 'f1']
                    + ([] if norm >= 2 else ['g0'])
                    + ['g1', 'h0', 'h1', 'i0', 'i1'])
        if norm >= 2:                       # s1 = 1 (ν-scaling chart)
            s_def, sp_ = 's - (x + s4*t)', ['s4']
        elif norm == 1:                     # s2 = s3 = 0 (rotation normal form)
            s_def, sp_ = 's - (s1*x + s4*t)', ['s1', 's4']
        else:
            s_def, sp_ = 's - (s1*x + s2*y + s3*z + s4*t)', ['s1', 's2', 's3', 's4']
        field_eqs, fparams = [], []
        for f in ('u', 'v', 'w', 'p'):
            if f == 'p' and norm >= 3:      # pressure trim: p1 = p2 = p3 = 0
                field_eqs.append('p - (p4*t + p5*Psi)')
                fparams += ['p4', 'p5']
            elif background or f == 'p':
                field_eqs.append('%s - (%s1*x + %s2*y + %s3*z + %s4*t + %s5*Psi)'
                                 % (f, f, f, f, f, f))
                fparams += ['%s%d' % (f, i) for i in range(1, 6)]
            elif f == 'w' and norm >= 1:    # w5 = 0 -> w ≡ 0
                field_eqs.append('w')
            elif f == 'v' and norm >= 2:    # v5 = 1 (λ-scaling chart)
                field_eqs.append('v - Psi')
            else:
                field_eqs.append('%s - %s5*Psi' % (f, f))
                fparams.append('%s5' % f)
        eqs = (['Psi[%s] - DPsi*s[%s]' % (c, c) for c in coords]
               + ['DPsi[%s] - DDPsi*s[%s]' % (c, c) for c in coords]
               + [' + '.join(ode_terms)]
               + [s_def]
               + field_eqs)
        # v_params (a ≡ 0 -> DEGENERATE) shrink with the amplitude vector and
        # are EMPTY from norm 2 on: with v5 = 1 the amplitude vector can never
        # vanish, so the a ≡ 0 test is vacuous — leaving it in would silently
        # change what a DEGENERATE verdict means.
        v_params = ([] if norm >= 2
                    else ['u5', 'v5'] if norm == 1
                    else ['u5', 'v5', 'w5'])
        # excluded locus, as data (forwarded by build_problem): what each
        # chart rung gives up relative to 25.3.  Rung 4 records none — it is
        # a different ansatz family, not a chart, so nothing "lifts" anyway.
        excludes = []
        if norm in (1, 2, 3):
            excludes += ['sigma = s1^2+s2^2+s3^2 = 0, s != 0 '
                         '(punctured null cone)']
        if norm in (2, 3):
            excludes += ['s1 = 0', 'v5 = 0', 'g0 = 0']
        if norm == 3:
            excludes += ['special cells where the 25.3 pseudo-division '
                         'cofactor (b0+b1*s)*DPsi + (c0+c1*s)*Psi + (g0+g1*s) '
                         'vanishes']
        return dict(kind='nssystem',
                    jets_dep=['u', 'v', 'w', 'p', 'DDPsi', 'DPsi', 'Psi', 's'],
                    equations=eqs, params=sp_ + fparams + odep,
                    v_params=v_params,
                    amp_params=[],
                    excludes=excludes)

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
    elif pde_name == 'helium':
        R1, R2, R12 = var('R1 R2 R12')
        PsiF = function('PsiF')(R1, R2, R12)
        H = (-1/2 * sum(diff(PsiF, Ri, 2) + 2/Ri * diff(PsiF, Ri) for Ri in (R1, R2))
             - (diff(PsiF, R12, 2) + 2/R12 * diff(PsiF, R12))
             - (R1**2 + R12**2 - R2**2)/(2*R1*R12) * diff(diff(PsiF, R12), R1)
             - (R2**2 + R12**2 - R1**2)/(2*R2*R12) * diff(diff(PsiF, R12), R2)
             - sum(2/Ri for Ri in (R1, R2)) * PsiF + 1/R12 * PsiF)
        rootsub = {}
    else:
        raise ValueError("unknown pde %r" % pde_name)
    expr = H - Evar * PsiF
    for k, v_ in rootsub.items():
        expr = expr.subs({k: v_})
    return PsiF, expr


def pde_params(pde_name):
    """The PDE's OWN constants (not ansatz parameters): ranked with the
    parameters in the lowest block, but NOT subject to the pconst constancy
    equations -- a PDE constant never appears in the ansatz, so it is never
    differentiated and needs no constancy relation."""
    if pde_name in ('hydrogen', 'helium'):
        return ['E']
    if pde_name == 'navier-stokes':
        return ['rho', 'mu']
    if pde_name == 'navier-stokes-nd':               # nondimensional: rho = mu = 1
        return []
    raise ValueError("unknown pde %r" % pde_name)


def pde_system(pde_name, coords):
    """The PDE system as a LIST of ring-parser strings.  The Schroedinger
    problems are one-element systems built through the sympy Hamiltonian
    path (build_pde_string)."""
    if pde_name in ('hydrogen', 'helium'):
        return [build_pde_string(pde_name, coords)]
    if pde_name in ('navier-stokes', 'navier-stokes-nd'):
        # Incompressible Navier-Stokes: continuity + three momentum equations.
        # These are already differential polynomials in the field jets u, v, w,
        # p -- no denominators to clear, so no sympy Hamiltonian detour.
        # 'navier-stokes-nd' is the nondimensional variant (rho = mu = 1).
        rho, mu = ('rho', 'mu') if pde_name == 'navier-stokes' else ('1', '1')
        mom = ('%(r)s*(%(f)s[t] + u*%(f)s[x] + v*%(f)s[y] + w*%(f)s[z])'
               ' + p[%(c)s] - %(m)s*(%(f)s[x,x] + %(f)s[y,y] + %(f)s[z,z])')
        return (['u[x] + v[y] + w[z]']
                + [mom % dict(r=rho, m=mu, f=f, c=c)
                   for f, c in (('u', 'x'), ('v', 'y'), ('w', 'z'))])
    raise ValueError("unknown pde %r" % pde_name)


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
def build_problem(pde_name, ansatz, ranking='orderly'):
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
    pparams = pde_params(pde_name)
    if ranking in ('elimination', 'block', 'elim'):
        # Block ranking: each jet its own block (highest first), the PDE's own
        # constants + params in a final low block.  The block comparison
        # dominates coordinate-order, so a high jet like DDPsi outranks ALL
        # lower-jet derivatives (e.g. DPsi[R1]).  The chain rule
        # DPsi[c]-DDPsi*v[c] then eliminates DDPsi directly instead of forcing
        # the cleared high-degree prolongation the orderly ranking does.
        DVAR = [[j] for j in jets] + [pparams + params]
    else:                                          # 'orderly' -- degrevlex
        DVAR = jets + pparams + params
    rk = dt.compute_ranking(IVAR, DVAR)
    R = rk.ring

    ansatz_eqs = [R(s) for s in ansatz_eqs_str]
    pconst = [R('%s[%s]' % (p, c)) for p in params for c in coords]
    pdes = [R(s) for s in pde_system(pde_name, coords)]

    return dict(R=R, rk=rk, coords=coords, roots=roots, jets=jets,
                tower=tower, order=spec.get('order'), params=params,
                v_params=v_params, amp_params=amp_params, ansatz_eqs=ansatz_eqs,
                pconst=pconst, pdes=pdes, pde_params=pparams,
                ansatz_eqs_str=ansatz_eqs_str,
                excludes=spec.get('excludes', []))
