# -*- mode: python -*-
#
# Parse the simple systems (cells) emitted by the DifferentialThomas run in
# hydrogen_thomas.log and translate each cell's equations / inequations from
# Maple function-notation into François Boulier's DifferentialAlgebra jet
# notation, so the paper's algorithm (reduce / project / decompose) can be run
# against each cell.  Shared by the prototype (plain python) and joca-thomas.sage.
#
# Maple <-> Sage/DifferentialAlgebra name map:
#   jets:   DDPs->DDPsi  DPs->DPsi  Ps->Psi  Vf->v  rho->rho_   (rho_ avoids the
#                                                                'r' remainder name)
#   params: V1..V4 -> v1..v4 ; a0,a1,b0,b1,c0,c1 unchanged ; E only in the PDE.
#
# Maple prints derivatives as diff(Ps(x,y,z),x), diff(diff(Ps(x,y,z),x),y), ...
# and the zero-derivative parameters as a0(x,y,z) etc.  Parameter-constancy
# equations  diff(<param>(x,y,z),v) = 0  are dropped (the DifferentialRing
# already declares the constants as parameters, so their derivatives vanish).

import re
import sympy

# --- the names, as they appear in the Maple log ---------------------------
JET_NAMES = {'DDPs': 'DDPsi', 'DPs': 'DPsi', 'Ps': 'Psi', 'Vf': 'v', 'rho': 'rho_'}
PARAM_NAMES = {'V1': 'v1', 'V2': 'v2', 'V3': 'v3', 'V4': 'v4',
               'a0': 'a0', 'a1': 'a1', 'b0': 'b0', 'b1': 'b1', 'c0': 'c0', 'c1': 'c1'}
IVARS = ('x', 'y', 'z')


def split_cells(logtext):
    """Return list of dicts {num, neqs, nineqs, eqs_str, ineqs_str} from the log."""
    blocks = re.split(r'--- cell (\d+) : (\d+) eqs, (\d+) ineqs ---', logtext)
    cells = []
    for i in range(1, len(blocks), 4):
        num, neqs, nineqs, body = (int(blocks[i]), int(blocks[i + 1]),
                                   int(blocks[i + 2]), blocks[i + 3])
        m_eq = re.search(r'EQS:\s*(\[.*?\])\s*\n\s*INEQS:', body, re.S)
        m_in = re.search(r'INEQS:\s*(\[.*?\])\s*(?:\n|$)', body, re.S)
        cells.append(dict(num=num, neqs=neqs, nineqs=nineqs,
                          eqs_str=m_eq.group(1) if m_eq else '[]',
                          ineqs_str=m_in.group(1) if m_in else '[]'))
    return cells


def _maple_list(s):
    """Split a Maple list string '[a, b, c]' into top-level elements."""
    s = s.strip()
    assert s[0] == '[' and s[-1] == ']', s[:40]
    s = s[1:-1]
    out, depth, cur = [], 0, []
    for ch in s:
        if ch in '([': depth += 1
        elif ch in ')]': depth -= 1
        if ch == ',' and depth == 0:
            out.append(''.join(cur)); cur = []
        else:
            cur.append(ch)
    if ''.join(cur).strip():
        out.append(''.join(cur))
    return [e.strip() for e in out]


def make_translator(jetbases, const_syms, ivar_syms):
    """
    Build a function maple_expr_str -> sympy expression in DifferentialAlgebra
    jet notation.  `jetbases` maps Sage jet name -> DifferentialAlgebra indexedbase
    object; `const_syms` maps param Sage name -> constant sympy symbol; `ivar_syms`
    maps 'x'/'y'/'z' -> sympy symbol.
    """
    x, y, z = ivar_syms['x'], ivar_syms['y'], ivar_syms['z']
    # undefined sympy Functions for every Maple base name
    funcs = {nm: sympy.Function(nm) for nm in list(JET_NAMES) + list(PARAM_NAMES)}

    def parse(expr_str):
        s = expr_str.replace('^', '**')
        ns = dict(funcs); ns.update(x=x, y=y, z=z, diff=sympy.diff)
        e = sympy.sympify(s, locals=ns, evaluate=True)
        # build xreplace map over the atoms actually present
        repl = {}
        for f in e.atoms(sympy.Function):
            if isinstance(f, sympy.Derivative):
                continue
            nm = f.func.__name__
            if nm in JET_NAMES:
                repl[f] = jetbases[JET_NAMES[nm]]            # order-0 jet
            elif nm in PARAM_NAMES:
                repl[f] = const_syms[PARAM_NAMES[nm]]        # constant
        for d in e.atoms(sympy.Derivative):
            base = d.expr.func.__name__
            # flatten (var, count) pairs into a flat list of derivation vars
            dvars = []
            for v, n in d.variable_count:
                dvars.extend([v] * int(n))
            if base in JET_NAMES:
                repl[d] = jetbases[JET_NAMES[base]][tuple(dvars)]
            elif base in PARAM_NAMES:
                repl[d] = sympy.Integer(0)                   # params are constant
        return sympy.expand(e.xreplace(repl))

    return parse


# Any raw Maple equation string mentioning one of these base names involves a
# jet (the unknown function or its derivatives); everything else is a pure
# parameter relation.  Used to skip the huge solved-form jet equations of the
# generic cells without sympy-parsing them.
_JET_RX = re.compile(r'\b(?:DDPs|DPs|Ps|Vf|rho)\b')


def cell_polys(cell, parse, const_syms, parse_diff=False):
    """
    Translate one cell.  By default only the (small) parameter relations and
    inequations are parsed; the huge solved-form jet equations of the generic
    cells are skipped (pass parse_diff=True to also translate them).  Returns:

      param_eqs   : sympy polynomials in the constants only (the cell's
                    parameter stratum: each must vanish).
      param_ineqs : constant-only inequation polynomials (must be != 0).
      jet_ineqs   : raw strings of inequations that involve jets (e.g. Ps != 0),
                    i.e. "nontrivial solution" conditions, recorded for report.
      diff_eqs    : (only if parse_diff) the jet-leader differential polynomials,
                    denominators cleared.
    """
    const_set = set(const_syms.values())
    param_eqs, param_ineqs, jet_ineqs, diff_eqs = [], [], [], []

    for eq in _maple_list(cell['eqs_str']):
        lhs_s, rhs_s = eq.split('=', 1)
        is_jet = bool(_JET_RX.search(eq))
        if not is_jet:
            # parameter relation; drop constancy eqs  diff(<param>,v) = 0
            if lhs_s.strip().startswith('diff(') and rhs_s.strip() == '0':
                continue
            num, _ = sympy.fraction(sympy.together(parse(lhs_s) - parse(rhs_s)))
            param_eqs.append(sympy.expand(num))
        elif parse_diff:
            num, _ = sympy.fraction(sympy.together(parse(lhs_s) - parse(rhs_s)))
            diff_eqs.append(sympy.expand(num))

    for ineq in _maple_list(cell['ineqs_str']):
        lhs_s = ineq.split('<>')[0].strip()
        if _JET_RX.search(ineq):
            jet_ineqs.append(lhs_s)
        else:
            num, _ = sympy.fraction(sympy.together(parse(lhs_s)))
            param_ineqs.append(sympy.expand(num))

    return dict(param_eqs=param_eqs, param_ineqs=param_ineqs,
                jet_ineqs=jet_ineqs, diff_eqs=diff_eqs)
