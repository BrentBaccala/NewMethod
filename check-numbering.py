#!/usr/bin/env python3
"""Check that the numbering of the submitted version has not drifted.

The version submitted to JOCA on 2026-09-15 has fixed theorem, section,
equation and figure numbers: a referee's "Lemma 14" must still be Lemma 14
in the revision.  frozen-numbering.txt records label -> number for every
\\label in that version; this script diffs the current NewMethod.aux
against it.

    ./check-numbering.py            report drift against the frozen snapshot
    ./check-numbering.py --freeze   overwrite the snapshot from the .aux
                                    (only when a renumbering is intended)

New labels are fine and are reported separately.  A number that moved, or
a frozen label that vanished, is an error (exit 1).
"""
import re
import sys
import os

AUX = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'NewMethod.aux')
FROZEN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'frozen-numbering.txt')


def group(s, i):
    """Return (contents, index past) of the brace group starting at s[i]=='{'."""
    assert s[i] == '{'
    depth, j = 0, i
    while j < len(s):
        if s[j] == '{':
            depth += 1
        elif s[j] == '}':
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise ValueError('unbalanced braces at %d' % i)


def read_aux(path):
    """label -> (printed number, anchor kind) for every \\newlabel."""
    text = open(path, encoding='utf-8', errors='replace').read()
    out = {}
    for m in re.finditer(r'\\newlabel\s*(?=\{)', text):
        name, i = group(text, m.end())
        if name.endswith('@cref'):      # cleveref's shadow entries
            continue
        body, _ = group(text, i)
        number, j = group(body, 0)
        kind = ''
        # hyperref's fifth field is the anchor, e.g. "theorem.14".
        fields = []
        while j < len(body) and len(fields) < 4:
            f, j = group(body, j)
            fields.append(f)
        if len(fields) >= 3:
            kind = fields[2].split('.', 1)[0]
        out[name] = (number, kind)
    return out


def write_frozen(entries, path):
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('# label\tnumber\tkind -- numbering of the submitted '
                 'version, regenerate only on purpose\n')
        for name in sorted(entries):
            number, kind = entries[name]
            fh.write('%s\t%s\t%s\n' % (name, number, kind))


def read_frozen(path):
    entries = {}
    for line in open(path, encoding='utf-8'):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        entries[parts[0]] = (parts[1], parts[2] if len(parts) > 2 else '')
    return entries


def main():
    if not os.path.exists(AUX):
        sys.exit('no %s -- build the paper first' % AUX)
    current = read_aux(AUX)

    if '--freeze' in sys.argv:
        write_frozen(current, FROZEN)
        print('froze %d labels into %s' % (len(current), FROZEN))
        return 0

    frozen = read_frozen(FROZEN)
    moved = [(n, frozen[n][0], current[n][0]) for n in sorted(frozen)
             if n in current and current[n][0] != frozen[n][0]]
    gone = [n for n in sorted(frozen) if n not in current]
    added = [n for n in sorted(current) if n not in frozen]

    for name, was, now in moved:
        print('MOVED   %-45s %s -> %s' % (name, was, now))
    for name in gone:
        print('MISSING %-45s was %s' % (name, frozen[name][0]))
    for name in added:
        print('new     %-45s %s (%s)' % (name, current[name][0],
                                         current[name][1]))
    if moved or gone:
        print('\n%d frozen number(s) drifted -- the submitted numbering is '
              'no longer intact.' % (len(moved) + len(gone)))
        return 1
    print('%d frozen numbers intact, %d new label(s).' % (len(frozen),
                                                          len(added)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
