#!/usr/bin/env python3
"""
NewMethod paper: the xi ("first") separated equation, in SymPy.

    d/dxi ( xi df1/dxi ) + [ (1/2) E xi - (1/4) m^2/xi + beta1 ] f1 = 0

Companion to f1-ode.wls (Wolfram) and f1-ode-sage.sage (Sage/Maxima).
See cas-comparison.tex for the write-up.

Bottom line: SymPy solves ONLY the beta1 = 0 reduction (Bessel).  For
beta1 != 0 it has no applicable solver -- it has no Whittaker functions and
no confluent-hypergeometric ODE solver -- and falls back to a truncated
Frobenius series.  Run with any Python that has sympy + mpmath.
"""

import sympy as sp
import mpmath as mp

xi = sp.symbols("xi", positive=True)          # xi = r + z >= 0 (parabolic coord.)
EE, beta1 = sp.symbols("E beta1", real=True)  # constants ...
m = sp.symbols("m", positive=True)            # ... everything except xi and f1
f1 = sp.Function("f1")


def ode(mval=m, bval=beta1):
    """The equation, in the divergence form printed in the paper."""
    return sp.Eq(
        sp.diff(xi * sp.diff(f1(xi), xi), xi)
        + (sp.Rational(1, 2) * EE * xi - sp.Rational(1, 4) * mval**2 / xi + bval)
        * f1(xi),
        0,
    )


def rule(label, mval, bval):
    """What solvers does SymPy think apply?  This is the whole story."""
    try:
        hints = sp.classify_ode(ode(mval, bval), f1(xi))
    except Exception as exc:                              # noqa: BLE001
        hints = ("classify_ode raised %s" % type(exc).__name__,)
    print(f"\nclassify_ode [{label}]:\n    {hints}")
    return hints


def solve(label, mval, bval):
    print(f"\ndsolve [{label}]:")
    try:
        print("   ", sp.dsolve(ode(mval, bval), f1(xi)).rhs)
    except Exception as exc:                              # noqa: BLE001
        print(f"    FAILED: {type(exc).__name__}: {str(exc)[:120]}")


print("sympy", sp.__version__, "/ mpmath", mp.__version__)
print("\nODE:", ode())

# ---------------------------------------------------------------- classify
# beta1 = 0 gets '2nd_linear_bessel'; beta1 != 0 gets only the series solver.
rule("general, beta1 symbolic", m, beta1)
rule("beta1 = 0, m symbolic", m, sp.Integer(0))
rule("beta1 = 0, m = 2", sp.Integer(2), sp.Integer(0))

# ------------------------------------------------------------------ solve
solve("beta1 = 0, m symbolic", m, sp.Integer(0))          # -> Bessel, correct
solve("beta1 = 0, m = 2", sp.Integer(2), sp.Integer(0))   # -> Bessel, correct
solve("general, m = 0", sp.Integer(0), beta1)             # -> series + O(xi**6)
solve("general, m symbolic", m, beta1)                    # -> TypeError

# ---------------------------------------------------- no Whittaker at all
print("\nWhittaker functions in sympy:",
      [n for n in dir(sp) if "whit" in n.lower()] or "NONE")
print("Kummer U / hyperu in sympy:",
      [n for n in dir(sp) if n.lower() in ("hyperu", "kummeru")] or "NONE")

# ------------------------------------------- can it VERIFY a closed form?
# Hand SymPy the Kummer-M solution it cannot derive.  checkodesol returns
# False -- not because the expression is wrong, but because SymPy cannot
# collapse the contiguous-relation combination of hyper() terms.
k = sp.sqrt(EE / 2)
a = (1 + m) / 2 + sp.I * beta1 / (2 * k)
closed_form = xi ** (m / 2) * sp.exp(-sp.I * k * xi) * sp.hyper((a,), (m + 1,), 2 * sp.I * k * xi)
ok, residual = sp.checkodesol(ode(), sp.Eq(f1(xi), closed_form))
print(f"\ncheckodesol on the Kummer-M closed form: {ok}"
      f"  (residual has {sp.count_ops(residual)} ops -- unsimplified, not wrong)")

# ------------------------------------------- mpmath settles it numerically
print("\nmpmath residual of that same closed form (40 dps):")
mp.mp.dps = 40


def residual_at(E, b1, mv, x=mp.mpf("1.3")):
    kk = mp.sqrt(mp.mpf(E) / 2)
    aa = (1 + mv) / mp.mpf(2) + 1j * b1 / (2 * kk)
    g = lambda t: t ** (mv / mp.mpf(2)) * mp.e ** (-1j * kk * t) * mp.hyp1f1(aa, mv + 1, 2j * kk * t)
    return abs(mp.diff(lambda t: t * mp.diff(g, t), x)
               + (E * x / 2 - mp.mpf(mv) ** 2 / (4 * x) + b1) * g(x))


for E in (0.7, 3.0):
    for b1 in (0.0, 1.4):
        for mv in (0, 2, 3):
            print(f"    E={E} beta1={b1} m={mv}: {mp.nstr(residual_at(E, b1, mv), 4)}")
