# Coherence (passivity) check for the generic cell that --generic-cell builds.
#
# The cell's equations are the ansatz equations + the constancy relations, taken
# as given.  A differential system only defines a normal form when it is
# PASSIVE: every integrability condition between two equations whose leaders are
# derivatives of the same dependent variable must reduce to zero modulo the
# system.  thomas-ansatz-solve.sage does not test this -- it only checks that
# leaders are distinct.  This script tests it.
#
# For equations p, q with leaders u_A, u_B on the same head u, and L = lcm(A,B),
# the Delta-polynomial is
#
#     Delta = sep(q) * D^{L-A}(p)  -  sep(p) * D^{L-B}(q)
#
# reduced to a fixpoint against the whole system.  Nonzero remainders are
# integrability conditions the generic cell is MISSING.
#
#   sage coherence-check.sage --pde navier-stokes --ansatz 25.3
import os, sys
sys.path.insert(0, os.path.expanduser('~/DifferentialThomas-sage'))
sys.path.insert(0, os.path.expanduser('~/sage-differential-polynomial/src'))
load(os.path.expanduser('~/Papers/NewMethod/ansatz-library.sage'))


def _argval(flag, default=None):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else default


PDE_NAME = _argval('--pde', 'navier-stokes')
A = _argval('--ansatz', '25.3')
ANSATZ = float(A) if '.' in str(A) else int(A)

prob = build_problem(PDE_NAME, ANSATZ)
R = prob['R']
COORDS = prob['coords']
eqs = list(prob['ansatz_eqs']) + list(prob['pconst'])
eqs = [e for e in eqs if not e.is_zero()]
print("%s / ansatz %s: %d equations, coords %s"
      % (PDE_NAME, ANSATZ, len(eqs), COORDS))


def full_prem(p, reductors, max_passes=64):
    r, h = p, R.one()
    for _ in range(max_passes):
        r2, h2 = R.differential_prem(r, reductors)
        h = h * h2
        if r2 == r:
            return r2
        r = r2
    return r


def split_leader(L):
    """'Psi[x,y]' -> ('Psi', ['x','y']);  'DDPsi' -> ('DDPsi', [])."""
    if L is None:
        return None, None
    if '[' not in L:
        return L, []
    head, rest = L.split('[', 1)
    return head, rest.rstrip(']').split(',')


# Group equations by the HEAD of their leader.  Only equations sharing a head
# can have an integrability condition between them.
by_head = {}
for e in eqs:
    head, idx = split_leader(e.leader())
    if head is None:
        continue
    by_head.setdefault(head, []).append((e, idx))

print("\nleader heads: %s"
      % ", ".join("%s(%d eq)" % (h, len(v)) for h, v in sorted(by_head.items())))

missing, checked = [], 0
for head, group in sorted(by_head.items()):
    for i in range(len(group)):
        for j in range(i + 1, len(group)):
            p, Aidx = group[i]
            q, Bidx = group[j]
            # multidegree of each leader, and their lcm
            from collections import Counter
            ca, cb = Counter(Aidx), Counter(Bidx)
            lcm = {c: max(ca.get(c, 0), cb.get(c, 0)) for c in set(ca) | set(cb)}
            # D^{L-A} p  and  D^{L-B} q
            dp, dq = p, q
            for c, n in lcm.items():
                for _ in range(n - ca.get(c, 0)):
                    dp = dp.differentiate(c)
                for _ in range(n - cb.get(c, 0)):
                    dq = dq.differentiate(c)
            delta = q.separant() * dp - p.separant() * dq
            rem = full_prem(delta, eqs)
            checked += 1
            if not rem.is_zero():
                missing.append((p.leader(), q.leader(), rem))

print("\nchecked %d integrability conditions" % checked)
if not missing:
    print("PASSIVE: every integrability condition reduces to zero.")
    print("The generic cell is coherent -- full_prem against it IS a normal form.")
else:
    print("NOT PASSIVE: %d integrability condition(s) do NOT reduce to zero."
          % len(missing))
    for la, lb, rem in missing:
        s = str(rem)
        print("\n  [%s , %s] ->\n    %s" % (la, lb, s[:400] + (" ..." if len(s) > 400 else "")))
    print("\nThese are equations the true generic cell has and this one does not,")
    print("so --generic-cell can UNDER-report: a PDE that would reduce to zero")
    print("modulo the completed cell may leave a nonzero remainder here, adding")
    print("spurious constraints and shrinking the membership locus.")
