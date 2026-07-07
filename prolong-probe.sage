# Pointwise existence probe for the bounded-prolongation route, v2.
# Build (A u {P})^[NMAX] ONCE (generators tagged by prolongation order, so
# every smaller N is a subset); then for each test point c* and each N,
# decide jet existence by 1 in J(c*) -- a parameter-free GB over GF(p).
# Unbuffered, timed at every stage.

import sys
import time

p = 32003
K = GF(p)
NMAX = int(sys.argv[sys.argv.index('--nmax')+1]) if '--nmax' in sys.argv else 4

T0 = time.time()
def stage(msg):
    print(f"[{time.time()-T0:8.1f}s] {msg}", flush=True)

const_names = ['E', 'v1', 'v2', 'v3', 'v4', 'a0', 'a1', 'b0', 'b1', 'c0', 'c1']
deps = ['Psi', 'DPsi', 'DDPsi', 'v']

def multi_indices(n):
    return [(i, j, t - i - j)
            for t in range(n + 1)
            for i in range(t + 1)
            for j in range(t - i + 1)]

midx = multi_indices(NMAX)

def jetname(w, m):
    s = 'x' * m[0] + 'y' * m[1] + 'z' * m[2]
    return w if not s else w + '_' + s

jet_list = [(w, m) for w in deps for m in midx]
names = ([jetname(w, m) for (w, m) in jet_list]
         + ['r', 'rinv', 'x', 'y', 'z'] + const_names)
R = PolynomialRing(K, names, order='degrevlex')
g = dict(zip(names, R.gens()))
jet = {(w, m): g[jetname(w, m)] for (w, m) in jet_list}
stage(f"ring built: {R.ngens()} variables (NMAX={NMAX})")

X, Y, Z = g['x'], g['y'], g['z']
E = g['E']
v1, v2, v3, v4 = g['v1'], g['v2'], g['v3'], g['v4']
a0, a1, b0, b1, c0, c1 = (g[s] for s in ['a0','a1','b0','b1','c0','c1'])
Psi, DPsi, DDPsi, v = (g[s] for s in deps)
r, rinv = g['r'], g['rinv']

succ = {}
for (w, m) in jet_list:
    succ[jet[(w, m)]] = [jet.get((w, tuple(m[i] + (1 if i == d else 0)
                                           for i in range(3))))
                         for d in range(3)]
succ[X] = [R.one(), R.zero(), R.zero()]
succ[Y] = [R.zero(), R.one(), R.zero()]
succ[Z] = [R.zero(), R.zero(), R.one()]
succ[r] = [X*rinv, Y*rinv, Z*rinv]
succ[rinv] = [-X*rinv**3, -Y*rinv**3, -Z*rinv**3]
for s in const_names:
    succ[g[s]] = [R.zero()] * 3

def Dtot(poly, d):
    out = R.zero()
    for w in poly.variables():
        s = succ[w][d]
        if s is None:
            raise ValueError(f"jet overflow: {w}")
        if s:
            out += poly.derivative(w) * s
    return out

def prolongations(eq, maxord):
    memo = {(0, 0, 0): eq}
    for m in multi_indices(maxord):
        if sum(m) == 0:
            continue
        d = next(i for i in range(3) if m[i] > 0)
        m0 = tuple(m[i] - (1 if i == d else 0) for i in range(3))
        memo[m] = Dtot(memo[m0], d)
    return memo

def e(w, *m):
    return jet[(w, tuple(m))]

solved = [
    (e('Psi',1,0,0) - DPsi*e('v',1,0,0), ('Psi',(1,0,0)), 1),
    (e('Psi',0,1,0) - DPsi*e('v',0,1,0), ('Psi',(0,1,0)), 1),
    (e('Psi',0,0,1) - DPsi*e('v',0,0,1), ('Psi',(0,0,1)), 1),
    (e('DPsi',1,0,0) - DDPsi*e('v',1,0,0), ('DPsi',(1,0,0)), 1),
    (e('DPsi',0,1,0) - DDPsi*e('v',0,1,0), ('DPsi',(0,1,0)), 1),
    (e('DPsi',0,0,1) - DDPsi*e('v',0,0,1), ('DPsi',(0,0,1)), 1),
    (v - (v1*X + v2*Y + v3*Z + v4*r), ('v',(0,0,0)), 0),
]
ODE = (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi
r_rel = r**2 - X**2 - Y**2 - Z**2
r_inv = r*rinv - 1
PDE = (e('Psi',2,0,0) + e('Psi',0,2,0) + e('Psi',0,0,2))*r + 2*Psi + 2*E*r*Psi
kept = [(ODE, 0), (r_rel, 0), (r_inv, 0), (PDE, 2)]

stage("prolonging ...")
rules = {}            # jet -> (rhs, order)
tagged = []           # (generator, order tag)
for eq, (w, m0), o in solved:
    for m, peq in prolongations(eq, NMAX - o).items():
        L = jet[(w, tuple(m0[i] + m[i] for i in range(3)))]
        ordL = sum(m0) + sum(m)
        rhs = L - peq
        if L in rules:
            d = rules[L][0] - rhs
            if d:
                tagged.append((d, ordL))
        else:
            rules[L] = (rhs, ordL)
for eq, o in kept:
    for m, peq in prolongations(eq, NMAX - o).items():
        tagged.append((peq, o + sum(m)))
stage(f"prolonged: {len(rules)} rules, {len(tagged)} generators")

# Resolve rules in one pass, dependency order: v-jets reference only the
# base; DPsi-jets reference v-jets and DDPsi-jets; Psi-jets reference
# DPsi-jets and v-jets.  Then substitute into the kept generators.
fam = {'v': 0, 'DPsi': 1, 'Psi': 2}
order_of = {jet[(w, m)]: (fam[w], sum(m)) for (w, m) in jet_list if w in fam}
resolved = {}
for L in sorted(rules, key=lambda u: order_of[u]):
    rhs = rules[L][0]
    pend = [u for u in rhs.variables() if u in resolved]
    resolved[L] = rhs.subs({u: resolved[u] for u in pend}) if pend else rhs
for L in resolved:
    assert all(u not in resolved for u in resolved[L].variables()), L
stage("rules resolved (one pass)")

def substitute(q):
    pend = [u for u in q.variables() if u in resolved]
    return q.subs({u: resolved[u] for u in pend}) if pend else q

tagged = [(substitute(q), t) for q, t in tagged]
tagged = [(q.subs({Psi: R.one()}), t) for q, t in tagged if q]
tagged = [(q, t) for q, t in tagged if q]
for q, t in tagged:
    assert all(u not in resolved for u in q.variables())
stage(f"substituted: {len(tagged)} generators; "
      f"per-order counts: { {N: sum(1 for _, t in tagged if t <= N) for N in range(2, NMAX+1)} }")

# ------------------------------------------------------------ test points
half = K(1)/K(2)
pts = [
 ('generic random 1',      [K.random_element() for _ in range(11)], False),
 ('new E=0 (parabolic)',   [0, 1,0,0,1, 0,1, 1,0, 1,0],             True),
 ('classical (E=1)',       [1, 0,0,0,1, 0,1, 2,0, 2,2],             True),
 ('bad-locus e^-r E=-1/2', [-half, 0,0,0,1, 0,0, 1,0, 1,0],         True),
 ('bad-locus 2s E=-1/8',   [-K(1)/K(8), 0,0,0,1, 0,0, -2,1, -2,half], True),
 ('bad-locus non-solving', [7, 0,0,0,1, 0,0, 1,2, 3,5],             False),
 ('v==0 stratum (E=2)',    [2, 0,0,0,0, 0,1, 1,0, 1,0],             False),
 ('generic random 2',      [K.random_element() for _ in range(11)], False),
]

for name, vals, solving in pts:
    sub = {g[c]: K(vv) for c, vv in zip(const_names, vals)}
    row = []
    for N in range(2, NMAX + 1):
        t = time.time()
        J = R.ideal([q.subs(sub) for q, tag in tagged if tag <= N])
        exists = R.one() not in J
        row.append(f"N={N}:{'yes' if exists else 'NO'}({time.time()-t:.1f}s)")
        stage(f"{name:24s} " + '  '.join(row))
    expect = 'solution' if solving else 'non-solution'
    print(f"  -> {name}: {'  '.join(row)}   [{expect}]", flush=True)
