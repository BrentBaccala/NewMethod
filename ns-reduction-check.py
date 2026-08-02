#!/usr/bin/env python3
# ns-reduction-check.py
# ------------------------------------------------------------------------
# Independent sympy check of the Navier-Stokes reduction written up in
# NewMethod.tex, section "Example: Incompressible Navier-Stokes".
#
# This does NOT use ansatz-library.sage or the Thomas machinery.  It
# substitutes Figure 8's ansatz into the incompressible Navier-Stokes
# equations by hand -- applying the chain rules
#
#     Psi_c = Psi' s_c,   Psi'_c = Psi'' s_c
#
# as a custom total-derivative operator -- and prints the residuals
# collected on the jets, then checks each coefficient the paper quotes.
# It is the oracle for the reduced system: anything the Sage pipeline
# produces for ansatz 25 must agree with it, cell-by-cell reductions aside.
#
# Run:  python3 ns-reduction-check.py     (needs only sympy)
#
# Author: Brent Baccala (AI assistant: Claude).  August 2026.

import sympy as sp

x,y,z,t = sp.symbols('x y z t')
rho,mu  = sp.symbols('rho mu')
s1,s2,s3,s4 = sp.symbols('s1 s2 s3 s4')
u1,u2,u3,u4,u5 = sp.symbols('u1 u2 u3 u4 u5')
v1,v2,v3,v4,v5 = sp.symbols('v1 v2 v3 v4 v5')
w1,w2,w3,w4,w5 = sp.symbols('w1 w2 w3 w4 w5')
p1,p2,p3,p4,p5 = sp.symbols('p1 p2 p3 p4 p5')
Ps,Pp,Ppp,Pppp = sp.symbols('Psi Psip Psipp Psippp')

s = s1*x + s2*y + s3*z + s4*t

def D(e, c):
    """total derivative w.r.t. coordinate c, using the ansatz chain rules
       Psi_c = Psi' s_c,  Psi'_c = Psi'' s_c,  Psi''_c = Psi''' s_c"""
    sc = sp.diff(s, c)
    return (sp.diff(e, c) + sp.diff(e, Ps)*Pp*sc
                          + sp.diff(e, Pp)*Ppp*sc
                          + sp.diff(e, Ppp)*Pppp*sc)

U = u1*x+u2*y+u3*z+u4*t; V = v1*x+v2*y+v3*z+v4*t
W = w1*x+w2*y+w3*z+w4*t; L = p1*x+p2*y+p3*z+p4*t
u = U + u5*Ps; v = V + v5*Ps; w = W + w5*Ps; p = L + p5*Ps

def lap(f): return sum(D(D(f,c),c) for c in (x,y,z))
def conv(f): return u*D(f,x) + v*D(f,y) + w*D(f,z)

cont = sp.expand(D(u,x)+D(v,y)+D(w,z))
momx = sp.expand(rho*(D(u,t)+conv(u)) + D(p,x) - mu*lap(u))
momy = sp.expand(rho*(D(v,t)+conv(v)) + D(p,y) - mu*lap(v))
momz = sp.expand(rho*(D(w,t)+conv(w)) + D(p,z) - mu*lap(w))

sigma = s1**2+s2**2+s3**2
adots = u5*s1 + v5*s2 + w5*s3

def show(name, e):
    print("== %s ==" % name)
    P_ = sp.Poly(e, Ps, Pp, Ppp, Pppp)
    for mono, coef in sorted(zip(P_.monoms(), P_.coeffs())):
        lbl = '*'.join(f for f,k in zip(['Psi','Psip','Psipp','Psippp'],mono) for _ in range(k)) or '1'
        print("  %-12s %s" % (lbl, sp.factor(sp.expand(coef))))
    print()

show('continuity', cont); show('x-momentum', momx)

def cf(e, mono):
    return sp.Poly(e, Ps, Pp, Ppp, Pppp).coeff_monomial(mono)
chk = [
 ("cont  Psip      == a.s",              cf(cont,Pp) - adots),
 ("momx  Psi*Psip  == rho*u5*(a.s)",     cf(momx,Ps*Pp) - rho*u5*adots),
 ("momy  Psi*Psip  == rho*v5*(a.s)",     cf(momy,Ps*Pp) - rho*v5*adots),
 ("momz  Psi*Psip  == rho*w5*(a.s)",     cf(momz,Ps*Pp) - rho*w5*adots),
 ("momx  Psipp     == -mu*u5*sigma",     cf(momx,Ppp) + mu*u5*sigma),
 ("momx  Psi       == rho*(u1u5+u2v5+u3w5)", cf(momx,Ps) - rho*(u1*u5+u2*v5+u3*w5)),
 ("momx  Psip      == rho*u5*s4+p5*s1+rho*u5*(s1U+s2V+s3W)",
      cf(momx,Pp) - (rho*u5*s4 + p5*s1 + rho*u5*(s1*U+s2*V+s3*W))),
 ("momx  jet-free  == rho*u4+p1+rho*(U*u1+V*u2+W*u3)",
      cf(momx,sp.Integer(1)) - (rho*u4+p1+rho*(U*u1+V*u2+W*u3))),
]
print("== checks ==")
for lbl, d in chk: print("  %-58s %s" % (lbl, sp.simplify(sp.expand(d)) == 0))

allowed = {(0,0,0,0),(1,0,0,0),(0,1,0,0),(1,1,0,0),(0,0,1,0)}
print("\n== jet monomials present ==")
for name,e in (('cont',cont),('momx',momx),('momy',momy),('momz',momz)):
    got = set(sp.Poly(e,Ps,Pp,Ppp,Pppp).monoms())
    print("  %-5s beyond {1,Psi,Psip,Psi*Psip,Psipp}: %s" % (name, sorted(got-allowed) or 'none'))
