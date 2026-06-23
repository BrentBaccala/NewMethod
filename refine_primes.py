#!/usr/bin/env python
# refine_primes.py <joca-thomas output file>
#
# Post-process a joca-thomas.sage run WITHOUT re-running it: among the reported
# solution varieties, drop every ideal that is a strict SUPERSET of another in
# the same list.  We want the largest varieties = the smallest ideals, so if
# I_a ⊊ I_b then V(I_a) ⊋ V(I_b): V(I_b) is already contained in V(I_a) and is
# redundant in the union.  Keep the ⊆-minimal ideals (each presented in its
# simplest form).  This is the same refinement now built into joca-thomas.sage;
# this script applies it to an already-produced output.

import sys, re
from sage.all import PolynomialRing, QQ

OUT = sys.argv[1]
text = open(OUT).read()

IDEAL_RX = re.compile(
    r'V:\s*Ideal \((?P<gens>.*?)\) of Multivariate Polynomial Ring in '
    r'(?P<vars>.+?) over Rational Field')
CELLS_RX = re.compile(r'\(from cells:\s*(?P<cells>[^)]*)\)')

def split_top(s):
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

def parse_block(block):
    """Return list of dicts {ideal, gens_str, cells} for each V: line in a block."""
    out = []
    lines = block.splitlines()
    for i, ln in enumerate(lines):
        m = IDEAL_RX.search(ln)
        if not m:
            continue
        varnames = [v.strip() for v in m.group('vars').split(',')]
        R = PolynomialRing(QQ, varnames)
        gens = [R(g) for g in split_top(m.group('gens'))]
        I = R.ideal(gens)
        cells = ''
        for ln2 in lines[i + 1:i + 3]:
            cm = CELLS_RX.search(ln2)
            if cm:
                cells = cm.group('cells').strip(); break
        out.append(dict(ideal=I, gens_str=m.group('gens'), cells=cells))
    return out

def subset(J, I):
    """True iff ideal J ⊆ ideal I (every generator of J lies in I)."""
    return all(g in I for g in J.gens())

def minimal(items):
    """Keep the ⊆-minimal ideals; return (kept, dropped) where dropped notes the
    minimal ideal that absorbs each redundant one."""
    kept, dropped = [], []
    for i, a in enumerate(items):
        redundant_under = None
        for j, b in enumerate(items):
            if i == j:
                continue
            # b ⊊ a  =>  V(a) ⊊ V(b), so a is redundant
            if subset(b['ideal'], a['ideal']) and not subset(a['ideal'], b['ideal']):
                redundant_under = b; break
            # exact duplicate: keep the earlier index, drop the later
            if (subset(b['ideal'], a['ideal']) and subset(a['ideal'], b['ideal'])
                    and j < i):
                redundant_under = b; break
        (dropped if redundant_under else kept).append((a, redundant_under))
    return [a for a, _ in kept], [(a, b) for a, b in dropped]

def refine(title, block):
    items = parse_block(block)
    if not items:
        return
    kept, dropped = minimal(items)
    print("=" * 72)
    print(f"{title}: {len(items)} reported -> {len(kept)} after dropping "
          f"{len(dropped)} redundant (superset) ideal(s)\n")
    for a in kept:
        print(f"  V: Ideal ({a['gens_str']})")
        print(f"       (from cells: {a['cells']})")
    if dropped:
        print("\n  dropped (superset of a kept ideal -> smaller variety, redundant):")
        for a, b in dropped:
            print(f"    ({a['gens_str']})")
            print(f"        ⊇ ({b['gens_str']})  [from cells {a['cells']}]")
    print()

# split the file into the NONTRIVIAL and TRIVIAL sections
nt = re.search(r'(NONTRIVIAL solution varieties.*?)(?:\n-{5,}|\Z)', text, re.S)
tr = re.search(r'(TRIVIAL \(Psi==0\) strata.*?)(?:\n-{5,}|\Z)', text, re.S)
if nt: refine("NONTRIVIAL solution varieties", nt.group(1))
if tr: refine("TRIVIAL (Psi==0) strata", tr.group(1))
