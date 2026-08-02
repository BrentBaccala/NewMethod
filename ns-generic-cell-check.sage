# Generic-cell reduction check for ansatz 25.x vs the ns-reduction-check.py
# oracle: reduce the four Navier-Stokes PDEs against the ansatz equations +
# constancy relations (the GENERIC cell's triangular system, before any Thomas
# splitting) and print the remainders collected on the jets.
import sys, os, re
import sympy

load(os.path.expanduser('~/Papers/NewMethod/ansatz-library.sage'))

ANSATZ = float(sys.argv[-1]) if '.' in sys.argv[-1] else int(sys.argv[-1])
PDE_NAME = sys.argv[-2]
prob = build_problem(PDE_NAME, ANSATZ)
R = prob['R']

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

# sympy symbols for everything
names = (prob['coords'] + prob['jets'] + prob['pde_params'] + prob['params'])
S = {n: sympy.Symbol(n) for n in names}

def blad_to_sympy(e):
    s = str(e)
    assert '_' not in re.sub(r'\b(DPsi|DDPsi)\b', '', s), \
        "derivative jet survived reduction: %s" % s
    return sympy.expand(sympy.sympify(s.replace('^', '**'), locals=S))

jets3 = [S['Psi'], S['DPsi'], S['DDPsi']]
labels = ['continuity', 'x-momentum', 'y-momentum', 'z-momentum']
rems = []
for lbl, P in zip(labels, prob['pdes']):
    rem, h = full_prem(P, reductors)
    rems.append(rem)
    print("== %s ==  (cofactor h = %s)" % (lbl, h))
    e = blad_to_sympy(rem)
    P_ = sympy.Poly(e, *jets3)
    for mono, coef in sorted(zip(P_.monoms(), P_.coeffs())):
        lb = '*'.join(f for f, k in zip(['Psi', 'DPsi', 'DDPsi'], mono)
                      for _ in range(k)) or '1'
        print("  %-14s %s" % (lb, sympy.factor(sympy.expand(coef))))
    print()

# jet-monomial support check (prediction 1), leader_deg 2 only
allowed = {(0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1)}
print("== jet monomials beyond {1, Psi, DPsi, Psi*DPsi, DDPsi} ==")
for lbl, rem in zip(labels, rems):
    got = set(sympy.Poly(blad_to_sympy(rem), *jets3).monoms())
    print("  %-12s %s" % (lbl, sorted(got - allowed) or 'none'))
