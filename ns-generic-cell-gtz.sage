# Generic-cell membership locus for ansatz 25.x: reduce the four PDEs against
# the ansatz equations + constancy relations (the generic cell), project onto
# the constants exactly as thomas-ansatz-solve.sage does, and take minimal
# associated primes.  This is what the driver would compute for the generic
# (no-parameter-equations) cell of the Thomas decomposition.
import sys, os, re
import sympy

load(os.path.expanduser('~/Papers/NewMethod/ansatz-library.sage'))

ANSATZ = float(sys.argv[-1]) if '.' in sys.argv[-1] else int(sys.argv[-1])
PDE_NAME = sys.argv[-2]
prob = build_problem(PDE_NAME, ANSATZ)
R = prob['R']
COORDS = prob['coords']
PARAMS = prob['params']
PDE_PARAMS = prob['pde_params']

reductors = list(prob['ansatz_eqs']) + list(prob['pconst'])

def full_prem(p, reductors, max_passes=64):
    r = R(p); h = R.one()
    for _ in range(max_passes):
        r2, h2 = R.differential_prem(r, reductors)
        h = h * h2
        if r2 == r:
            return r2, h
        r = r2
    return r, h

PolyRing = PolynomialRing(QQ, names=(COORDS + prob['jets'] + PDE_PARAMS + PARAMS))
PolyRing_constants = [PolyRing(p) for p in (PDE_PARAMS + PARAMS)]

def build_system_of_equations(eqn, constants):
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

_all = []
for P in prob['pdes']:
    rem, h = full_prem(P, reductors)
    rp = PolyRing(str(rem).replace('^', '**'))
    if not rp.is_zero():
        _all.extend(build_system_of_equations(rp, PolyRing_constants))
eqns = tuple(set(_all))
print("combined system: %d equations" % len(eqns))
I = ideal(list(eqns))
import time
t0 = time.time()
primes = I.minimal_associated_primes()
print("GTZ %.1fs -> %d minimal primes" % (time.time() - t0, len(primes)))
for P in primes:
    print("  V:", sorted([str(g) for g in P.gens()]))
