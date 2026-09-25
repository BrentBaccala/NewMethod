#!/usr/bin/env python3
"""Copy the <N> fixed theorem numbers from the current NewMethod.tex onto the
matching environments of an older version (argv[1], rewritten in place), for
make-diff.sh.  The old version's numbered environments are the \\begin{...}
lines outside comment environments, the \\ifNSexample blocks and % comments;
they must match the current file's <N>-annotated ones type for type."""
import re, sys

ENVS = 'theorem|lemma|corollary|proposition|algorithm|example|exercise|definition'
old_path, new_path = sys.argv[1], sys.argv[2]

fixed = re.findall(r'^[^%\n]*\\begin\{(' + ENVS + r')\}<(\d+)>',
                   open(new_path).read(), re.M)

lines = open(old_path).read().split('\n')
skip = 0          # nesting depth of comment / \ifNSexample regions
sites = []
for i, line in enumerate(lines):
    code = re.sub(r'(?<!\\)%.*', '', line)
    if re.search(r'\\begin\{comment\}|^\s*\\ifNSexample\b', code):
        skip += 1
    if not skip:
        for m in re.finditer(r'\\begin\{(' + ENVS + r')\}(?!<)', code):
            sites.append((i, m.group(1)))
    if re.search(r'\\end\{comment\}', code) or (skip and re.match(r'\s*\\fi\b', code)):
        skip -= 1

if [e for _, e in sites] != [e for e, _ in fixed]:
    sys.exit('diff-number-old.py: environment sequences differ:\n  old %s\n  new %s'
             % ([e for _, e in sites], [e for e, _ in fixed]))
for (i, env), (_, n) in zip(sites, fixed):
    lines[i] = lines[i].replace('\\begin{%s}' % env, '\\begin{%s}<%s>' % (env, n), 1)
open(old_path, 'w').write('\n'.join(lines))
