# -*- mode: python -*-
#
# Sage script for the paper in the Journal of Computational Algebra.
#
# This is the differential-Thomas variant of joca-rg.sage.  Instead of running
# the paper's algorithm once against the single (generic) ansatz, it loads the
# disjoint simple systems ("cells") of the differential Thomas decomposition of
# the hydrogen ansatz -- computed by Maple's DifferentialThomas package and
# logged in ~/thomas-experiments/hydrogen_thomas.log -- and runs the algorithm
# once per cell:
#
#     for each Thomas cell:
#         specialize the ansatz to the cell's parameter stratum
#         reduce the PDE modulo that ansatz        (Ritt full reduction)
#         gather the remainder's coefficients in the parameters   (projection)
#         compute minimal associated primes        (GTZ)
#         drop primes excluded by the cell's inequations          (saturation)
#     union the surviving primes over all cells
#
# Each cell's parameter stratum (verified empirically across all 73 cells) is a
# subset of the constants {a0,a1,b0,b1,c0,c1,v1,v2,v3,v4} set to zero; so the
# cell-specialized ansatz is just the original ansatz with those constants set
# to 0 (the chain-rule relations survive intact, so Ritt reduction reduces
# fully).  Because the cells are *disjoint* and each is a squarefree regular
# simple system, the per-cell saturation is valid pointwise -- and the union
# over the cells is the complete decomposition.  This is the paper's
# completeness-by-composition thesis made executable: differential Thomas as
# the completeness oracle, the basic algorithm run cell-by-cell.
#
# Author: Brent Baccala
# Date: June 2026
# Tested on Ubuntu 24 with Sage 10.7 and DifferentialAlgebra 1.14.
#
#   sage joca-thomas.sage [--log PATH] [--cells 1,45,46] [--verbose-remainder]

import os, sys

import sympy

try:
    import DifferentialAlgebra
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(ex.msg + "\nInstall it with '%pip install DifferentialAlgebra'")

# --- command-line options -------------------------------------------------
def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default

LOG = _argval('--log', os.path.expanduser('~/thomas-experiments/hydrogen_thomas.log'))
CELLS_ARG = _argval('--cells')                       # e.g. "1,45,46"; default all
WANT_CELLS = set(int(n) for n in CELLS_ARG.split(',')) if CELLS_ARG else None
VERBOSE_REM = '--verbose-remainder' in sys.argv

# --- fixed-chain prolongation (the existential route) ----------------------
# --prolong K   run the fixed-chain prolongation per cell to order <= K,
#               instead of (well: in addition to) the universal projection.
# --prolong-stage {swell,closure,eliminate}
#               how far to take it.  'swell' only prolongs and reports sizes
#               (cheap, always terminates); 'closure' adds the D-closure test
#               (needs a Groebner basis); 'eliminate' adds the elimination to
#               K[c] (may be expensive -- this is the known cliff, at per-cell
#               scale).  Default: swell.
# --prolong-modp P
#               do the *commutative* stages (closure/elimination) over GF(P)
#               instead of QQ.  The differential reduction stays over QQ.
PROLONG = _argval('--prolong')
PROLONG_ORDER = int(PROLONG) if PROLONG else 0
PROLONG_STAGE = _argval('--prolong-stage', 'swell')
PROLONG_MODP = int(_argval('--prolong-modp', 0))
assert PROLONG_STAGE in ('swell', 'closure', 'eliminate'), \
    "--prolong-stage must be swell, closure or eliminate"

import functools, time
print = functools.partial(print, flush=True)

# the cell parser (thomas_cells.py) ships next to this script; also look in the
# log's directory so the parser can live with the log if preferred
_here = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv and sys.argv[0] else ''
for _d in (_here, os.path.dirname(LOG)):
    if _d and _d not in sys.path:
        sys.path.insert(0, _d)
import thomas_cells as tc

# Claude Sonnet 4's solution to printing "DPsi" as "\Psi'" in LaTeX
def patch_latex_varify():
    from sage.misc.latex import latex_varify
    import sage.misc.latex
    original = latex_varify
    def custom(a, is_fname=False):
        return r"\Psi'" if a == "DPsi" else original(a, is_fname=is_fname)
    sage.misc.latex.latex_varify = custom
patch_latex_varify()

# --- ring / PDE / ansatz (identical to joca.sage) -------------------------
x, y, z = sympy.var('x,y,z')
E = sympy.var('E')
v1, v2, v3, v4 = sympy.var('v1,v2,v3,v4')
a0, a1, b0, b1, c0, c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E, v1, v2, v3, v4, a0, a1, b0, b1, c0, c1]

Psi, DPsi, DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')

DiffRing = DifferentialAlgebra.DifferentialRing(derivations=[x, y, z],
                                                blocks=[[DDPsi, DPsi, Psi, v, r], constants],
                                                parameters=constants,
                                                notation='jet')

PDE = -int(1)/int(2)*(Psi[x, x] + Psi[y, y] + Psi[z, z])*r - Psi - E*r*Psi
print("PDE:", PDE)

ansatz0 = [Psi[x] - DPsi*v[x],
           Psi[y] - DPsi*v[y],
           Psi[z] - DPsi*v[z],
           DPsi[x] - DDPsi*v[x],
           DPsi[y] - DDPsi*v[y],
           DPsi[z] - DDPsi*v[z],
           (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
           v - (v1*x + v2*y + v3*z + v4*r),
           r**2 - x**2 - y**2 - z**2]
ansatz0 = list(map(sympy.expand, ansatz0))

# --- Sage polynomial ring (identical to joca.sage) ------------------------
PolyRing = PolynomialRing(QQ, names=[str(indet) for indet in DiffRing.indets(selection='all')])
PolyRing_constants = list(map(PolyRing, constants))


def build_system_of_equations(eqn, constants):
    """Factor each term into constant * non-constant parts, gather like
    non-constant monomials; return the (constant-only) coefficient polynomials.
    Identical to joca.sage / joca-rg.sage."""
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


# --- the differential-Thomas cells ----------------------------------------
jetbases = {'DDPsi': DDPsi, 'DPsi': DPsi, 'Psi': Psi, 'v': v, 'rho_': r}
const_syms = {'v1': v1, 'v2': v2, 'v3': v3, 'v4': v4,
              'a0': a0, 'a1': a1, 'b0': b0, 'b1': b1, 'c0': c0, 'c1': c1}
parse = tc.make_translator(jetbases, const_syms, {'x': x, 'y': y, 'z': z})

with open(LOG) as fh:
    cells = tc.split_cells(fh.read())
print(f"\nLoaded {len(cells)} Thomas cells from {LOG}\n" + "=" * 72)


def prime_key(P):
    """A canonical, sortable string for a prime ideal (for dedup / union)."""
    return tuple(sorted(str(g) for g in P.gens()))


def _ideal_subset(J, I):
    """True iff ideal J ⊆ ideal I (every generator of J lies in I)."""
    return all(g in I for g in J.gens())


def keep_minimal(d):
    """Keep only the ⊆-minimal ideals in a {key: (ideal, cells)} dict.

    We report the largest varieties, i.e. the smallest ideals: if I_a ⊊ I_b then
    V(I_a) ⊋ V(I_b), so V(I_b) is already contained in V(I_a) and is redundant in
    the union.  Drop every ideal that strictly contains another in the list (the
    superset ideals); the survivors are each presented in their simplest form.
    Returns (kept_dict, dropped) with dropped = [(key, P, cells, absorber_key)].
    """
    items = list(d.items())
    kept, dropped = {}, []
    for i, (ki, (Pi, ci)) in enumerate(items):
        absorber = None
        for j, (kj, (Pj, cj)) in enumerate(items):
            if i == j:
                continue
            sub_ji = _ideal_subset(Pj, Pi)
            if sub_ji and not _ideal_subset(Pi, Pj):          # Pj ⊊ Pi  -> drop Pi
                absorber = kj; break
            if sub_ji and _ideal_subset(Pi, Pj) and j < i:    # exact dup -> keep earlier
                absorber = kj; break
        if absorber is None:
            kept[ki] = (Pi, ci)
        else:
            dropped.append((ki, Pi, ci, absorber))
    return kept, dropped


# --- the fixed-chain prolongation (the existential route) -------------------
#
# The universal (membership) route above stops at the Ritt remainder rho and
# collects its coefficients.  The existential route needs the *differential
# closure* of T_j u {P}, and rho alone does not carry it: the constants-cutting
# conditions are the reduced derivatives D^k rho.  (Toy example: A: Psi''-c*Psi,
# P: Psi'-Psi.  rho = Psi'-Psi projects universally to the empty set and
# existentially to all of C -- but D(rho) reduced is (c-1)*Psi, so Ex = {c=1}.
# Verified: the fixed-chain code below reproduces exactly that.)
#
# The point of the *fixed-chain* form is that rho is never handed back to the
# differential engine as a new equation.  Adjoining it would promote a parametric
# jet to a leader, un-passivate the chain, and force a full re-completion
# (Delta-polynomials, splitting, a growing jet space) -- which is exactly what
# exhausts memory on cells 1, 2, 25, 26, 28.  Instead the cell is left exactly as
# Thomas produced it and used purely as a rewriting system:
#
#   * T_j is passive, so d/dx, d/dy, d/dz extend to derivations of the FIXED
#     algebra  R_j = QQ[c][x,y,z,r]/(r^2-x^2-y^2-z^2)[parametric jets].  Every
#     denominator that shows up (the ODE's initial a0+a1*v, the separant r) is an
#     inequation of the cell, so the divisions are legal.  We build the derivation
#     TABLE once per cell -- reduce D_d(u) for each single parametric jet u, which
#     is tiny -- and thereafter never call the differential engine again.  (BLAD in
#     the prolongation loop measured ~2000x slower: 200s vs 0.1s for one step.)
#
#   * so each D^k rho reduces back into the SAME few non-constant generators
#     (x,y,z,r,Psi,DPsi for the generic cell).  No jet space to expand.
#
# The decisive structural fact for THIS ansatz: the ansatz and the PDE are linear
# homogeneous in Psi, so rho and every D^k rho are LINEAR forms in the parametric
# jets:   gen_i = sum_j M_ij(c,x,y,z,r) * u_j.   Hence, with m = #parametric jets:
#
#   MEMBERSHIP (universal):  M == 0 identically  -- every entry vanishes.
#                            (This is precisely the paper's coefficient collection.)
#   EXISTENCE  (existential): rank M <= m-1 at SOME base point, i.e. all m x m
#                            minors of M vanish somewhere on the quadric.  The
#                            minors involve NO jets -- only x,y,z,r and c -- so the
#                            elimination drops the jet variables for free.
#
# Mem <= Ex is then immediate (M == 0 => all minors vanish everywhere).  The two
# readings are told apart by *how* the minors vanish: identically (membership) vs
# on an isolated base-point locus (a partial stratum).

def parametric_jets(spec):
    """Which of Psi, DPsi, DDPsi survive Ritt reduction against this cell?
    Those are the cell's free Taylor data.  (Generic cell: Psi, DPsi.  On the
    a0=a1=0 cells the ODE degenerates to first order, DPsi becomes a leader, and
    the parametric set collapses to {Psi} alone.)"""
    out = []
    for jet in (Psi, DPsi, DDPsi):
        _, red = DiffRing.differential_prem(jet, spec)
        if sympy.expand(red - jet) == 0:
            out.append(jet)
    return out


RING_SYMS = set(str(i) for i in DiffRing.indets(selection='all'))


def derivation_table(spec, params):
    """D_d(w) for every generator w of R_j, as (numerator, denominator) pairs.
    Computed by reducing the derivative of a SINGLE symbol -- cheap (ms).

    Also the DECISIVE VALIDITY CHECK for the whole fixed-chain construction.
    Lemma 2 says a *passive* cell's derivations close on the fixed ring R_j: the
    derivative of a parametric jet is either principal (rewritten by the chain)
    or parametric (already in R_j).  If reduction instead throws up a jet that is
    NOT a generator of R_j, the chain is not passive / not of finite type here,
    the "parametric" set is really infinite, and the construction is invalid.

    This is exactly the case on the degenerate strata.  joca-thomas.sage builds
    each cell's chain by zeroing the cell's constants in the RAW ansatz, which is
    the completed Thomas chain only where the ansatz was already passive (the
    generic cell).  On cells 26/27 (a0=a1=0) the ODE loses its DDPsi term
    altogether, so DDPsi is unconstrained -- and so are DDPsi[x], DDPsi[xx], ...:
      reduce(D_x(DDPsi))  ->  ... DDPsi[x] ...      <-- escapes the ring
    The real Thomas cell carries the extra integrability equations that close
    this; the re-specialized ansatz does not.  Harmless for the universal route
    (Ritt reduction still reduces the PDE fully -- which is why the recorded
    29-cell projection is sound), fatal for the existential one.  Refuse.
    """
    tbl = {}
    for d in (x, y, z):
        for w in [r] + list(params):
            h, num = DiffRing.differential_prem(DiffRing.differentiate(w, d), spec)
            escaped = sorted(s for s in map(str, sympy.sympify(num).free_symbols)
                             if s not in RING_SYMS)
            if escaped:
                raise ValueError(
                    f"derivation escapes the fixed ring: reduce(D_{d}({w})) "
                    f"involves {escaped}, which are not generators of R_j.  The "
                    f"cell chain is not passive (it is the re-specialized raw "
                    f"ansatz, not the completed Thomas cell), so the parametric "
                    f"set is misidentified and the existential route is invalid "
                    f"here.  Feed the cell's own equations "
                    f"(thomas_cells.cell_polys(parse_diff=True)) instead.")
            tbl[(w, d)] = (sympy.expand(num), sympy.expand(h))
    return tbl


def prolong_cell(spec, rem, Z, ineq_polys, order, stage, modp):
    """Fixed-chain prolongation of `rem` against the passive cell `spec`.
    Never modifies `spec`; never calls the differential engine after the table.

    IMPORTANT: `spec` must be the cell's *passive* chain.  joca-thomas.sage
    currently builds it by substituting the cell's zeroed constants into the raw
    ansatz, which is the completed chain ONLY on cells where the ansatz was
    already passive (the generic cell).  On degenerate strata the real Thomas
    cell carries extra integrability equations, and the specialized ansatz is not
    passive -- so the free data is misidentified and the rank test is invalid.
    We refuse to report in that case rather than emit a wrong answer."""
    t_all = time.time()
    out = dict(orders=[], params=None, closed=None, elim=None, elim_primes=None,
               passive=None, error=None)
    if rem == 0:
        out['closed'] = True                  # membership cell: nothing to prolong
        return out
    try:
        params = parametric_jets(spec)
        out['params'] = [str(p) for p in params]
        tbl = derivation_table(spec, params)

        K = QQ if not modp else GF(modp)
        S = PolynomialRing(K, names=PolyRing.variable_names(), order='degrevlex')
        g = {str(n): v for n, v in zip(S.variable_names(), S.gens())}
        def toS(p):
            return S(str(sympy.expand(p)).replace('**', '^'))
        X = [g['x'], g['y'], g['z']]
        rS = g['r']
        Q = rS**2 - g['x']**2 - g['y']**2 - g['z']**2
        pj = [g[str(p)] for p in params]

        # the derivation, as a polynomial operator (denominators cleared)
        num = {k: toS(v[0]) for k, v in tbl.items()}
        den = {k: toS(v[1]) for k, v in tbl.items()}
        def D(p, d, di):
            m = S.one()
            for w in [r] + list(params):
                m = m.lcm(den[(w, d)]) if den[(w, d)] != 0 else m
            acc = m * p.derivative(X[di])
            for w in [r] + list(params):
                dw = num[(w, d)] * (m // den[(w, d)])
                acc += dw * p.derivative(g[str(w)])
            return acc.reduce([Q])

        def st(p):
            return (len(p.monomials()), p.total_degree()) if p != 0 else (0, 0)

        # ---- PASSIVITY TEST ------------------------------------------------
        # A cell is passive iff its derivations commute.  With the scaled
        # (denominator-cleared) operator Dt = m*D, algebra gives the polynomial
        # identity
        #   m^3 * [D_i,D_j](u) = m*(Dt_i Dt_j u - Dt_j Dt_i u)
        #                        - Dt_i(m)*Dt_j(u) + Dt_j(m)*Dt_i(u),
        # so the RHS vanishes iff the cell is passive (m is a unit on the cell).
        # This is the Delta-polynomial condition in fixed-chain language, and it
        # costs a handful of reductions of single symbols.
        #
        # It matters: joca-thomas.sage builds `spec` by zeroing the cell's
        # constants in the RAW ansatz.  That is the completed Thomas chain only
        # where the ansatz was already passive.  On degenerate strata the real
        # cell carries extra integrability equations, `spec` is NOT passive, the
        # "parametric" jets are not free, and the rank test below is invalid.
        # Refuse rather than report a wrong answer.
        def mult(di):
            m = S.one()
            for w in [r] + list(params):
                if den[(w, (x, y, z)[di])] != 0:
                    m = m.lcm(den[(w, (x, y, z)[di])])
            return m
        passive, why = True, "derivations commute"
        for (i, j) in ((0, 1), (0, 2), (1, 2)):
            di, dj = (x, y, z)[i], (x, y, z)[j]
            mi, mj = mult(i), mult(j)
            if mi != mj:
                passive, why = False, "scaled derivations have different multipliers"
                break
            m = mi
            for u in pj:
                lhs = (m * (D(D(u, dj, j), di, i) - D(D(u, di, i), dj, j))
                       - D(m, di, i) * D(u, dj, j)
                       + D(m, dj, j) * D(u, di, i)).reduce([Q])
                if lhs != 0:
                    passive, why = False, f"[D_{di},D_{dj}]({u}) != 0"
                    break
            if not passive:
                break
        out['passive'] = (passive, why)
        if not passive:
            out['error'] = ("cell chain is NOT passive (" + why + "): `spec` is "
                            "the specialized raw ansatz, not the completed Thomas "
                            "chain, so the parametric jets are misidentified and "
                            "the rank test would be invalid.  Wire the cell's own "
                            "equations (thomas_cells.cell_polys(parse_diff=True)) "
                            "before trusting the existential route on this cell.")
            out['secs'] = time.time() - t_all
            return out

        rho = toS(rem)
        gens = [rho]; cur = [rho]
        n0, d0 = st(rho)
        out['orders'].append(dict(k=0, ngens=1, terms=[n0], maxdeg=d0, secs=0.0))
        for k in range(1, order + 1):
            t0 = time.time(); new = []
            for gg in cur:
                for di, d in enumerate((x, y, z)):
                    p = D(gg, d, di)
                    if p != 0:
                        new.append(p)
            if not new:
                out['closed'] = True; break
            ss = [st(p) for p in new]
            out['orders'].append(dict(k=k, ngens=len(new), terms=[s[0] for s in ss],
                                      maxdeg=max(s[1] for s in ss),
                                      secs=time.time() - t0))
            gens += new; cur = new

        # --- the coefficient matrix M: gens are linear forms in the param jets
        M = []
        for p in gens:
            rowv = []
            for u in pj:
                rowv.append(p.coefficient({u: 1, **{w: 0 for w in pj if w != u}}))
            lin = sum(c * u for c, u in zip(rowv, pj))
            if lin != p:
                out['error'] = "generators are NOT linear in the parametric jets"
                out['secs'] = time.time() - t_all
                return out
            M.append(rowv)
        out['nrows'] = len(M); out['ncols'] = len(pj)

        if stage == 'swell':
            out['secs'] = time.time() - t_all
            return out

        # --- m x m minors: the existence condition (rank <= m-1 somewhere) -----
        m = len(pj)
        t0 = time.time()
        minors = []
        for combo in Combinations(range(len(M)), m):
            det = matrix(S, [M[i] for i in combo]).determinant().reduce([Q])
            if det != 0:
                minors.append(det)
        out['nminors'] = len(minors)
        out['minor_secs'] = time.time() - t0
        out['identically_zero'] = (len(minors) == 0)
        # membership is M == 0 (rank 0), NOT "all minors vanish" (rank <= m-1)
        out['M_is_zero'] = all(e == 0 for rowv in M for e in rowv)

        if stage == 'eliminate':
            t0 = time.time()
            base = [g['x'], g['y'], g['z'], rS]
            I = S.ideal(minors + [Q]) if minors else S.ideal(Q)
            for s in [rS] + [toS(f) for f in ineq_polys]:
                if s != 0 and not s.is_unit():
                    I = I.saturation(S.ideal(s))[0]
            J = I.elimination_ideal(base)
            out['elim'] = J
            out['elim_secs'] = time.time() - t0
            if not modp:
                out['elim_primes'] = J.minimal_associated_primes()
    except Exception as ex:
        out['error'] = f"{type(ex).__name__}: {ex}"
    out['secs'] = time.time() - t_all
    return out


# Cache reduction+decomposition by parameter stratum (cells that zero the same
# constants share an identical specialized ansatz, hence identical reduction and
# primes; only their inequations -- and thus the pruning -- differ).
strata_cache = {}
union_primes = {}            # prime_key -> (ideal, contributing cells)  -- nontrivial
trivial_primes = {}          # prime_key -> (ideal, contributing cells)  -- Psi==0 strata

# The Maple base name of the order-0 wavefunction jet (Ps -> Psi).
_PSI_MAPLE_BASE = next(k for k, v in tc.JET_NAMES.items() if v == 'Psi')

def cell_forces_psi_zero(cell):
    """True iff the cell's defining equations include  <Psi>(x,y,z) = 0 -- i.e.
    DifferentialThomas ITSELF forced the wavefunction's order-0 jet to vanish (a
    Psi==0 stratum).  This is authoritative and immune to the Ritt-multiplier
    artifact that fools the reduction-based test below: on a stratum where v==0
    the PDE's multiplier h = 64*(a0+a1*v)*r^2 vanishes, so differential_prem
    returns a false remainder==0 AND the Psi-vs-ansatz reduction keeps Psi a live
    indeterminate -- yet the PDE forces Psi=0 (it reads -(1+E*r)*Psi there).
    Cell 14 (the spurious prime TH_E) is exactly this case; Thomas records it as
    `Ps(x,y,z) = 0`, which we honor here."""
    target = f"{_PSI_MAPLE_BASE}(x,y,z)"
    for eq in tc._maple_list(cell['eqs_str']):
        if '=' not in eq:
            continue
        lhs, rhs = (s.replace(' ', '') for s in eq.split('=', 1))
        if rhs == '0' and lhs == target:
            return True
    return False

for cell in cells:
    num = cell['num']
    if WANT_CELLS is not None and num not in WANT_CELLS:
        continue
    cp = tc.cell_polys(cell, parse, const_syms)
    Z = sorted((p for p in cp['param_eqs']), key=str)        # constants set to 0
    Zkey = tuple(map(str, Z))

    if Zkey not in strata_cache:
        sub = {p: 0 for p in Z}
        spec = [e for e in (sympy.expand(a.subs(sub)) for a in ansatz0) if e != 0]
        h, rem = DiffRing.differential_prem(PDE, spec)
        # Does the specialized ansatz force the wavefunction to vanish?  If Psi
        # reduces to 0 modulo the stratum's ansatz, the only "solution" there is
        # the trivial Psi == 0 -- a degenerate stratum the paper's generic
        # reduction never reports (it keeps Psi a live indeterminate).
        _, psi_rem = DiffRing.differential_prem(Psi, spec)
        # v ≡ 0 stratum: if all of v1..v4 vanish, the ansatz forces Psi spatially
        # constant (Psi_x = DPsi*v_x = 0, ...), and the PDE on a constant Psi reads
        # -(1+E*r)*Psi, which has NO nonzero solution -> Psi ≡ 0.  This is the
        # degenerate locus where the Ritt multiplier h = 64*(a0+a1*v)*r^2 vanishes
        # through v=0, so differential_prem returns a false remainder==0 with Psi
        # kept live; the reduction-based test below cannot see it.  (Cells 14/15/16
        # are this case -> the spurious primes TH_E and its sub-strata.)
        Zideal = ideal([PolyRing(p) for p in Z]) if Z else ideal(PolyRing.zero())
        v_zero = all(PolyRing(vv) in Zideal for vv in (v1, v2, v3, v4))
        trivial = (psi_rem == 0) or v_zero
        if rem == 0:
            eqns = ()
        else:
            eqns = build_system_of_equations(PolyRing(rem), PolyRing_constants)
        # ideal = projected coefficient conditions  +  the stratum (Z = 0)
        gens = list(eqns) + [PolyRing(p) for p in Z]
        I = ideal(gens) if gens else ideal(PolyRing.zero())
        primes = I.minimal_associated_primes()
        strata_cache[Zkey] = dict(spec_len=len(spec), rem=rem, eqns=eqns,
                                  primes=primes, trivial=trivial)

    sc = strata_cache[Zkey]
    # Per-cell triviality: the cached reduction-based test (per parameter stratum)
    # OR the cell's own Ps(x,y,z)=0 equation.  The latter catches the v==0 strata
    # where the multiplier artifact makes the reduction test wrongly report the
    # stratum as a nontrivial solution (e.g. cell 14 -> the spurious TH_E).
    cell_trivial = sc['trivial'] or cell_forces_psi_zero(cell)

    # cell inequations: a prime is excluded when an inequation polynomial lies
    # in it (then the inequation vanishes identically on that component -- the
    # cell forbids it).  param_ineqs are in QQ[constants][x,y,z]; jet_ineqs
    # (e.g. Ps != 0) are "nontrivial solution" conditions, not parameter ones.
    ineq_polys = [PolyRing(f) for f in cp['param_ineqs']]
    survivors = []
    for P in sc['primes']:
        if P.is_one():                       # unit ideal => empty variety
            continue
        if any(g in P for g in ineq_polys):  # inequation vanishes identically here
            continue
        survivors.append(P)

    tag = "TRIVIAL (Psi==0 forced)" if cell_trivial else "nontrivial"
    print(f"\n--- cell {num}: zero {{{', '.join(Zkey) or '(none, generic)'}}}; "
          f"ansatz {sc['spec_len']} eqs; "
          f"{len(cp['param_ineqs'])} param-ineqs, {len(cp['jet_ineqs'])} jet-ineqs; {tag} ---")
    if VERBOSE_REM:
        print("  remainder:", sc['rem'])
    if sc['rem'] == 0 and not cell_trivial:
        print("  PDE reduces to 0: the whole stratum solves the PDE (nontrivially)")
    print("  system of equations:", *sc['eqns'], sep="\n    " if sc['eqns'] else " (none)")
    if not survivors:
        print("  surviving solution varieties: NONE (all pruned by inequations / empty)")
    bucket = trivial_primes if cell_trivial else union_primes
    for P in survivors:
        print("   V:", P, "" if not cell_trivial else "   [trivial Psi==0]")
        key = prime_key(P)
        bucket.setdefault(key, (P, []))[1].append(num)

    # --- the existential route on this cell -------------------------------
    if PROLONG_ORDER:
        sub = {p: 0 for p in Z}
        spec = [e for e in (sympy.expand(a.subs(sub)) for a in ansatz0) if e != 0]
        pr = prolong_cell(spec, sc['rem'], Z, cp['param_ineqs'],
                          PROLONG_ORDER, PROLONG_STAGE, PROLONG_MODP)
        print(f"  [prolong] fixed-chain, order<={PROLONG_ORDER}, "
              f"stage={PROLONG_STAGE}"
              + (f", GF({PROLONG_MODP})" if PROLONG_MODP else ", QQ"))
        if pr['params'] is not None:
            print(f"    parametric jets of this cell: {pr['params']}")
        if pr.get('passive') is not None:
            ok, why = pr['passive']
            print(f"    passivity ([D_i,D_j] on the parametric jets): "
                  f"{'PASSIVE' if ok else 'NOT PASSIVE'}  ({why})")
        for o in pr['orders']:
            tt = o['terms']
            print(f"    k={o['k']}: {o['ngens']} gens; "
                  f"terms {min(tt)}..{max(tt)}; maxdeg {o['maxdeg']}; {o['secs']:.1f}s")
        if pr.get('nrows'):
            print(f"    coefficient matrix M: {pr['nrows']} rows x {pr['ncols']} "
                  f"cols (generators are linear in the parametric jets)")
        if pr.get('nminors') is not None:
            m = pr['ncols']
            if pr['identically_zero']:
                # rank M <= m-1 EVERYWHERE.  NOT membership -- membership is
                # M == 0 (rank 0).  Rank m-1 leaves a 1-dim kernel, i.e. a proper
                # subfamily solving at every c: a partial stratum over all of C.
                print(f"    all {m}x{m} minors vanish IDENTICALLY => rank M <= {m-1} "
                      f"everywhere => a proper subfamily solves for EVERY c")
                print(f"    (this is NOT membership: membership needs M == 0, i.e. "
                      f"rank 0.  M == 0? {pr['M_is_zero']})")
            else:
                print(f"    {pr['nminors']} nonzero {m}x{m} minors "
                      f"({pr['minor_secs']:.1f}s) => existence is a base-point "
                      f"condition; eliminate to get the stratum")
        if pr['closed']:
            print("    D-closed (derivatives vanish)")
        if pr['elim'] is not None:
            print(f"    ELIMINATION IDEAL in K[c] ({pr['elim_secs']:.1f}s):")
            for gg in pr['elim'].gens():
                print("       ", gg)
        if pr['elim_primes'] is not None:
            for P in pr['elim_primes']:
                print("    E:", P)
        if pr['error']:
            print("    [prolong] ERROR:", pr['error'])
        print(f"    [prolong] total {pr['secs']:.1f}s")

# --- union over all cells --------------------------------------------------
def dump_union(title, d):
    # Refine to the largest varieties: drop ideals that strictly contain another
    # in the list (their variety is already covered -> redundant in the union).
    kept, dropped = keep_minimal(d)
    note = f"; {len(dropped)} redundant superset ideal(s) dropped" if dropped else ""
    print(f"\n{title} ({len(kept)} maximal varieties{note}):\n")
    for key, (P, cells_for) in sorted(kept.items(), key=lambda kv: str(kv[0])):
        print("  V:", P)
        print("       (from cells:", ", ".join(map(str, sorted(set(cells_for)))) + ")")
    if dropped:
        print("\n  dropped (superset of a kept ideal => smaller, redundant variety):")
        kbygens = {k: P for k, (P, _) in d.items()}
        for key, P, cells_for, absorber in sorted(dropped, key=lambda t: str(t[0])):
            print("    V:", P, "  (from cells:",
                  ", ".join(map(str, sorted(set(cells_for)))) + ")")
            print("        ⊇", kbygens[absorber])

print("\n" + "=" * 72)
scope = 'selected' if WANT_CELLS else 'all'
dump_union(f"NONTRIVIAL solution varieties over {scope} cells", union_primes)
print("\n" + "-" * 72)
dump_union(f"TRIVIAL (Psi==0) strata over {scope} cells", trivial_primes)
