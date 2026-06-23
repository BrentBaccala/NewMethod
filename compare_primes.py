from sage.all import PolynomialRing, QQ
V = ['x','y','z','DDPsi','DPsi','Psi','v','r','E','v1','v2','v3','v4','a0','a1','b0','b1','c0','c1']
R = PolynomialRing(QQ, V)
g = {n: R(n) for n in V}
def I(*s): return R.ideal([R(e) for e in s])
# joca.sage (NewMethod generic) primes
NM = {
 'NM1': I('a0','v4','v3','v2','v1'),
 'NM2': I('c1','c0','b1','b0','a1','a0'),
 'NM3': I('v1^2+v2^2+v3^2','a1','a0','v4'),
 'NM4': I('v4*c0-b0','E*c0-v4*c1','-v4^2*c1+E*b0','b1','a1-1/2*b0','a0','v3','v2','v1'),
 'NM5': I('v4*c0-b0','v1^2+v2^2+v3^2-v4^2','c1','b1','a1-b0','a0','E'),
}
# joca-thomas (refined) primes
TH = {
 'TH_A': I('-b1*c0+b0*c1','v4*c1-b1','v4*c0-b0','a1','a0','v3','v2','v1','E+1/2'),
 'TH_B': I('v4*c0-b0','E*c0-v4*c1','-v4^2*c1+E*b0','b1','a1-1/2*b0','a0','v3','v2','v1'),
 'TH_C': I('c0^2+4*b0*c1','4*v4*c1+c0','v4*c0-b0','b1+1/2*c0','a1','a0','v3','v2','v1','E+1/8'),
 'TH_D': I('v4*c0-b0','v1^2+v2^2+v3^2-v4^2','c1','b1','a1-b0','a0','E'),
 'TH_E': I('b0','a1','a0','v4','v3','v2','v1'),
}
def sub(J,K): return all(p in K for p in J.gens())   # J ⊆ K
print("Exact matches (J==K):")
for nn,J in NM.items():
    for tn,K in TH.items():
        if sub(J,K) and sub(K,J): print(f"  {nn} == {tn}")
print("\nContainments (NM ⊊ TH  => NM is a LARGER variety):")
for nn,J in NM.items():
    for tn,K in TH.items():
        if sub(J,K) and not sub(K,J): print(f"  {nn} ⊊ {tn}   (V({nn}) ⊋ V({tn}))")
print("\nContainments (TH ⊊ NM  => TH is a LARGER variety):")
for tn,K in TH.items():
    for nn,J in NM.items():
        if sub(K,J) and not sub(J,K): print(f"  {tn} ⊊ {nn}   (V({tn}) ⊋ V({nn}))")
print("\nNM primes with NO containment relation to any TH prime:")
for nn,J in NM.items():
    rel = any(sub(J,K) or sub(K,J) for K in TH.values())
    if not rel: print(f"  {nn}: {J.gens()}")
print("TH primes with NO containment relation to any NM prime:")
for tn,K in TH.items():
    rel = any(sub(K,J) or sub(J,K) for J in NM.values())
    if not rel: print(f"  {tn}: {K.gens()}")

print("\n=== UNION (total solution set) comparison ===")
from functools import reduce
Jnm = reduce(lambda a,b: a.intersection(b), NM.values())
Jth = reduce(lambda a,b: a.intersection(b), TH.values())
# V(Jnm) = U V(NM_i);  V(Jnm) ⊆ V(Jth)  iff  Jth ⊆ rad(Jnm)
def variety_subset(Jbig, Jsmall):   # V(Jsmall) ⊆ V(Jbig)  iff  Jbig ⊆ rad(Jsmall)
    return all(g in Jsmall.radical() for g in Jbig.gens())
nm_in_th = variety_subset(Jth, Jnm)   # U V(NM) ⊆ U V(TH) ?
th_in_nm = variety_subset(Jnm, Jth)   # U V(TH) ⊆ U V(NM) ?
print(f"  U V(NM) ⊆ U V(TH): {nm_in_th}")
print(f"  U V(TH) ⊆ U V(NM): {th_in_nm}")
print(f"  unions equal as varieties: {nm_in_th and th_in_nm}")
