# Ansatz-5 (hydrogen) partial-stratum search.
#
# Ansatz 5:  Psi = zeta(v),  v = v1 x + v2 y + v3 z + v4 r,  r^2 = x^2+y^2+z^2,
#            2nd-order ODE  L = A2 d^2 + L1 d + L0,  d = d/dv,
#            A2 = a0+a1 v,  L1 = b0+b1 v,  L0 = c0+c1 v.
#
# Reducing (H - E) Psi mod the ansatz and splitting on the independent radical r
# (K := v1^2+v2^2+v3^2 - v4^2) gives TWO first-order operators on zeta:
#   R1 = A1 d + B1 ,  A1 = -K L1 ,          B1 = 2E A2 - K L0
#   R0 = A0 d + B0 ,  A0 = 2 v4 (A2 - v L1), B0 = 2 (A2 - v4 v L0)
# multiplier h = A2 = a0+a1 v.
#
# A partial stratum (h != 0) needs a common nonzero solution zeta of L, R1, R0,
# i.e. ord GCRD(L,R1,R0) = 1.  With zeta' = rho*zeta, rho = -B1/A1:
#   proportionality  A1 B0 - A0 B1 = 0     (R0 selects the same mode)
#   Riccati          A2(rho'+rho^2) + L1 rho + L0 = 0   (that mode solves L)
# each an identity in v -> vanishing of every v-coefficient.
# We ask whether the resulting variety has a point with h not identically 0
# (i.e. (a0,a1) != (0,0)) and base-splitting genericity K != 0, v4 != 0.

Kp = PolynomialRing(QQ, 'a0,a1,b0,b1,c0,c1,K,v4,E')
a0,a1,b0,b1,c0,c1,K,v4,E = Kp.gens()
Rv.<v> = PolynomialRing(Kp)

A2 = a0 + a1*v
L1 = b0 + b1*v
L0 = c0 + c1*v
A1 = -K*L1
B1 = 2*E*A2 - K*L0
A0 = 2*v4*(A2 - v*L1)
B0 = 2*(A2 - v4*v*L0)

# derivative wrt v of a poly in Rv (coeffs are constants)
def dv(p):
    return p.derivative(v)

# proportionality:  A1 B0 - A0 B1  (poly in v)
prop = A1*B0 - A0*B1

# Riccati cleared by A1^2:  rho = -B1/A1
#   A1^2 (rho'+rho^2) = B1^2 - (B1' A1 - B1 A1')      [since rho'=-(B1'A1-B1A1')/A1^2]
#   A1^2 * Riccati = A2(B1^2 - B1' A1 + B1 A1') - L1 B1 A1 + L0 A1^2
ric = A2*(B1^2 - dv(B1)*A1 + B1*dv(A1)) - L1*B1*A1 + L0*A1^2

# collect all v-coefficients as generators of the ideal in Kp
gens = []
for p in (prop, ric):
    gens += list(p.coefficients())      # coefficients are elements of Kp
gens = [Kp(g) for g in gens if g != 0]

I = Kp.ideal(gens)
print("number of generators:", len(gens))
print("ideal dimension (affine, all params free):", I.dimension())

# Saturate to enforce h not identically zero: remove the component a0=a1=0.
# and genericity K!=0, v4!=0.  Saturate by K, v4, and by the ideal (a0,a1).
Isat = I
for f in (K, v4):
    Isat = Isat.saturation(Kp.ideal(f))[0]
# remove {a0=a1=0}: saturate by the irrelevant-ish ideal (a0,a1) via elimination trick
Isat_h = Isat.saturation(Kp.ideal(a0,a1))[0]

print("\nAfter saturating by K, v4 (genericity):")
print("  dim:", Isat.dimension(), " is <1>? ", Isat == Kp.ideal(Kp.one()))
print("After ALSO saturating by (a0,a1)  [h not identically 0]:")
print("  dim:", Isat_h.dimension(), " is <1>? ", Isat_h == Kp.ideal(Kp.one()))

if Isat_h != Kp.ideal(Kp.one()):
    print("\n  --> NONEMPTY: h!=0 partial strata CAN exist. Groebner basis:")
    for g in Isat_h.groebner_basis():
        print("     ", g)
else:
    print("\n  --> EMPTY: no h!=0 partial stratum for ansatz 5 on Schrodinger",
          "(off K=0, v4=0).")

# Also report the E-content: eliminate everything but E to see if strata force E.
print("\nElimination to E only (from Isat, genericity-saturated):")
elimE = Isat.elimination_ideal([a0,a1,b0,b1,c0,c1,K,v4])
print("  ", elimE.gens())
