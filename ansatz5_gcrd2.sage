# Follow-up: (1) confirm the L1 != 0 branch is empty; (2) handle b0=b1=0 correctly.
Kp = PolynomialRing(QQ, 'a0,a1,b0,b1,c0,c1,K,v4,E')
a0,a1,b0,b1,c0,c1,K,v4,E = Kp.gens()
Rv.<v> = PolynomialRing(Kp)

A2 = a0 + a1*v; L1 = b0 + b1*v; L0 = c0 + c1*v
A1 = -K*L1;  B1 = 2*E*A2 - K*L0
A0 = 2*v4*(A2 - v*L1);  B0 = 2*(A2 - v4*v*L0)
def dv(p): return p.derivative(v)
prop = A1*B0 - A0*B1
ric  = A2*(B1^2 - dv(B1)*A1 + B1*dv(A1)) - L1*B1*A1 + L0*A1^2
gens = [Kp(g) for p in (prop,ric) for g in p.coefficients() if g!=0]
I = Kp.ideal(gens)
Isat = I.saturation(Kp.ideal(K))[0].saturation(Kp.ideal(v4))[0]
Isat_h = Isat.saturation(Kp.ideal(a0,a1))[0]

# (1) L1 != 0 branch: also saturate by (b0,b1).  Empty  ==>  no strata with L1!=0.
branch_L1 = Isat_h.saturation(Kp.ideal(b0,b1))[0]
print("L1 != 0 branch  is <1> (empty)? ", branch_L1 == Kp.ideal(Kp.one()))

# (2) b0=b1=0 case done correctly.  Then A1=0, R1=B1 (order 0): need B1=0.
#     B1 = 2E A2 - K L0 = 0 identically  ==>  L0 = (2E/K) A2.
#     Remaining: common solution of  L: A2 z'' + L0 z = 0   and R0: A0 z' + B0 z = 0.
#     With L0=(2E/K)A2 the ODE is z'' + (2E/K) z = 0.  R0 gives z'/z = rho0 = -B0/A0.
#     Common solution  <=>  rho0 solves the ODE's Riccati:  rho0' + rho0^2 + 2E/K = 0.
# Substitute b0=b1=0 and L0=(2E/K)A2 (clear K):  work in Kp2 without b's.
Kp2 = PolynomialRing(QQ, 'a0,a1,K,v4,E'); a0b,a1b,Kb,v4b,Eb = Kp2.gens()
Rv2.<w> = PolynomialRing(Kp2)
A2b = a0b + a1b*w
# L0 = (2E/K)A2  ->  clear K: work with K*Riccati etc.  R0 with b=0:
#   A0 = 2 v4 A2 ,  B0 = 2(A2 - v4 w L0) = 2 A2 (1 - v4 w (2E/K)) ; rho0 = -B0/A0
#   = -(1 - (2E v4/K) w)/v4 = -1/v4 + (2E/K) w .
# Riccati (clear K^2):  K^2(rho0' + rho0^2) + 2 E K = 0.
rho0_num = -Kb + 2*Eb*v4b*w          # = K*v4*rho0 ; so rho0 = rho0_num/(K v4)
# rho0 = rho0_num/(K v4).  rho0' = (2 E v4)/(K v4)=2E/K ; rho0^2 = rho0_num^2/(K v4)^2
# (K v4)^2 (rho0'+rho0^2+2E/K) = (K v4)^2*(2E/K) + rho0_num^2 + (K v4)^2*(2E/K)
#   = 2 rho0_num^2-free... just build it directly:
lhs = (Kb*v4b)^2*( (2*Eb*Kb)/Kb )  # placeholder; do it cleanly below
# clean: multiply Riccati rho0'+rho0^2+2E/K = 0 by (K v4)^2:
#   (K v4)^2 rho0' = (K v4)^2 * 2E/K = 2 E K v4^2
#   (K v4)^2 rho0^2 = rho0_num^2
#   (K v4)^2 * 2E/K = 2 E K v4^2
ric0 = 2*Eb*Kb*v4b^2 + rho0_num^2 + 2*Eb*Kb*v4b^2
print("\nb0=b1=0 case: Riccati*(Kv4)^2 as poly in w:", ric0)
sys0 = [Kp2(g) for g in ric0.coefficients() if g!=0]
J0 = Kp2.ideal(sys0).saturation(Kp2.ideal(Kb))[0].saturation(Kp2.ideal(v4b))[0]
J0h = J0.saturation(Kp2.ideal(a0b,a1b))[0]
print("b0=b1=0 branch, h!=0, is <1> (empty)? ", J0h == Kp2.ideal(Kp2.one()))
print("  (its groebner basis:", list(J0h.groebner_basis()), ")")
