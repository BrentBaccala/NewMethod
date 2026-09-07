# Ansatz 5 vs hydrogen ON THE RADIAL LOCUS  v1=v2=v3=0  (v = v4 r).
#
# This is the locus that ansatz5_gcrd.sage / ansatz5_gcrd2.sage EXCLUDE
# (they saturate by K and v4 and assume |v_lin| != 0 so that v and r split).
# Here v and r are dependent, no splitting occurs, and the Schrodinger
# equation reduces to a single second-order operator -- the classical radial
# equation -- so the partial-stratum test is a plain two-operator GCRD.
#
# Backs Theorem "the gap on the radial locus" in NewMethod.tex
# (\S The Membership and Consistency Loci for this Example).
#
# Normalization: v4 = 1 (rescaling v only rescales the coefficients).
#   M = v d^2 + 2 d + (2 + 2E v)                     <- reduced PDE
#   L = (a0+a1 v) d^2 + (b0+b1 v) d + (c0+c1 v)      <- ansatz top element
# One right-Euclidean step of L by M, cleared by v, leaves  Rt = A d + B.
# ord GCRD(L,M) = 1  <=>  A != 0  and  rho = -B/A solves M's Riccati.

P.<a0,a1,b0,b1,c0,c1,E,v> = PolynomialRing(QQ, 8)
A = v*(b0 + b1*v) - 2*(a0 + a1*v)
B = v*(c0 + c1*v) - (a0 + a1*v)*(2 + 2*E*v)
d = lambda f: f.derivative(v)

# (A^2) * [ v(rho' + rho^2) + 2 rho + 2 + 2E v ]  with rho = -B/A
F = v*(-d(B)*A + B*d(A) + B^2) - 2*B*A + (2 + 2*E*v)*A^2
J = P.ideal([F.coefficient({v: k}) for k in range(F.degree(v) + 1)])

IA = P.ideal([A.coefficient({v: k}) for k in range(3)])   # A != 0
Js = J.saturation(IA)[0].saturation(P.ideal([a0, a1]))[0] # L genuinely 2nd order

print("dimension:", Js.dimension())
el = Js.elimination_ideal([a0, a1, b0, b1, c0, c1, v])
print("possible E:", [g.factor() for g in el.gens()])
# ==>  (2E+1)(8E+1)(18E+1) :   E = -1/2, -1/8, -1/18 = -1/(2n^2), n = 1,2,3.
#
# The bound n <= 3 is structural: a common solution forces rho = -B/A, and the
# ns Coulomb mode has n-1 poles (the zeros of its degree-(n-1) polynomial
# factor), while first-degree ODE coefficients cap deg A at 2.
#
# The three families of L are obtained by back-substituting each root and
# solving the linear conditions L(zeta_n) = 0; see the theorem statement.
# For n >= 4 the only ansatz-5 operator annihilating the ns mode is M itself,
# which is membership, not a partial stratum.
