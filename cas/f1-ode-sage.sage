#!/usr/bin/env sage
# NewMethod paper: the xi ("first") separated equation, in Sage / Maxima.
#
#     d/dxi ( xi df1/dxi ) + [ (1/2) E xi - (1/4) m^2/xi + beta1 ] f1 = 0
#
# Companion to f1-ode.wls (Wolfram) and f1-ode-sympy.py (SymPy).
# See cas-comparison.tex for the write-up.
#
# Bottom line: Sage's desolve(..., contrib_ode=True) -- i.e. Maxima's
# contrib_ode package -- SOLVES the general case in closed form, in terms of
# hypergeometric_M / hypergeometric_U.  Two caveats, both checked below:
#   (1) plain desolve (Maxima's ode2 alone) does NOT solve it;
#   (2) Maxima returns the xi^(-m/2) Frobenius branch, whose hypergeometric_M
#       half is undefined for positive integer m (second parameter 1-m <= 0)
#       -- exactly the physically relevant case.  Wolfram returns the
#       xi^(+m/2) branch, which has no such degeneracy.
#
# Run:  ~/miniforge3/envs/sage/bin/sage f1-ode-sage.sage

import mpmath as mp

xi = var('xi')                       # parabolic coordinate, xi = r + z >= 0
EE = var('EE'); b1 = var('b1'); m = var('m')   # constants
f = function('f1')(xi)

ode = diff(xi*diff(f, xi), xi) + (EE*xi/2 - m^2/(4*xi) + b1)*f == 0
print("Sage", sage.version.version)
print("Maxima", maxima.version())
print("\nODE:", ode)

# NOTE: ivar=xi is mandatory here.  Without it desolve raises
# "Unable to determine independent variable" -- the ODE carries four symbols
# and Sage will not guess which is the independent one.

print("\n=== plain desolve (Maxima ode2 only) ===")
try:
    print("   ", desolve(ode, f, ivar=xi))
except Exception as ex:
    print("    FAILED: %s: %s" % (type(ex).__name__, str(ex)[:160]))

print("\n=== desolve(contrib_ode=True), general case ===")
sol = desolve(ode, f, ivar=xi, contrib_ode=True)
for s in sol:
    print("   ", s)

print("\n=== desolve(contrib_ode=True), beta1 = 0 ===")
for s in desolve(ode.subs(b1 == 0), f, ivar=xi, contrib_ode=True):
    print("   ", s)

print("\n=== desolve(contrib_ode=True), beta1 = 0, m = 2 ===")
for s in desolve(ode.subs(b1 == 0).subs(m == 2), f, ivar=xi, contrib_ode=True):
    print("   ", s)

# ---------------------------------------------------------------------------
# Numeric verification of Maxima's general solution.
#
# Done in mpmath rather than Sage: differentiating Sage's hypergeometric_M
# symbolically hits "symbolic division by zero" in its _derivative_ method
# (it forms a/b with b = 1-m).  mpmath has no such trouble.
#
#   f = exp(s*xi/2) * F(a, 1-m, -s*xi) / xi^(m/2),   s = sqrt(2)*sqrt(-E)
#   a = -(s*beta1 + E*m - E)/(2E)
# ---------------------------------------------------------------------------
mp.mp.dps = 30

def maxima_branch(E, beta1, mv, which):
    s = mp.sqrt(2)*mp.sqrt(mp.mpc(-E))
    a = -(s*beta1 + E*mv - E)/(2*E)
    F = mp.hyp1f1 if which == 'M' else mp.hyperu
    return lambda t: mp.e**(s*t/2) * F(a, 1 - mv, -s*t) / t**(mp.mpf(mv)/2)

def residual(E, beta1, mv, which, x=mp.mpf('1.3')):
    g = maxima_branch(E, beta1, mv, which)
    return abs(mp.diff(lambda t: t*mp.diff(g, t), x)
               + (E*x/2 - mp.mpf(mv)^2/(4*x) + beta1)*g(x))

print("\n=== numeric residual of Maxima's branch (mpmath, 30 dps) ===")
for which in ('M', 'U'):
    print("  hypergeometric_%s:" % which)
    for E in (-1, 3):
        for beta1 in (0, 1.4):
            for mv in (0, 0.7, 2, 3):
                try:
                    r = mp.nstr(residual(E, beta1, mv, which), 4)
                except Exception as ex:
                    r = "UNDEFINED (%s)" % type(ex).__name__
                print("    E=%-2s beta1=%-4s m=%-3s -> %s" % (E, beta1, mv, r))
