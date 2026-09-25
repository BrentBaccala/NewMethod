#!/bin/sh
# Build marked-up "what changed since submission" PDFs with latexdiff.
#
#   ./make-diff.sh [REV]        (REV defaults to the tag "submitted")
#
# Produces, in this directory (all gitignored):
#   NewMethod-diff-bars.pdf   change bars in the margin; new text plain,
#                             removed text struck through
#   NewMethod-diff-color.pdf  new text green, removed text red + struck through
#
# latexdiff is not packaged on samsung; it's the upstream release unpacked
# under ~/opt/latexdiff (override with LATEXDIFF=...).
set -e
cd "$(dirname "$0")"
REV=${1:-submitted}
LATEXDIFF=${LATEXDIFF:-$HOME/opt/latexdiff/latexdiff-1.3.4/latexdiff-fast}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

git show "$REV:NewMethod.tex" > "$tmp/old.tex"
# The current source pins each theorem-like environment's number with <N>.
# Give the old source the same <N> annotations, so latexdiff sees no change
# there (and doesn't wrap a \DIFadd around the <N>, which would break it).
python3 diff-number-old.py "$tmp/old.tex" NewMethod.tex

for style in bars color; do
  perl "$LATEXDIFF" --encoding=utf8 --preamble="diff-preamble-$style.tex" \
 "$tmp/old.tex" NewMethod.tex \
       > "NewMethod-diff-$style.tex" 2> "$tmp/latexdiff.err" \
    || { cat "$tmp/latexdiff.err"; exit 1; }
  # (its "Wide character" warnings come from temp files used only for
  # matching, not from the output, so they are dropped)
  for pass in 1 2 3; do
    pdflatex -interaction=nonstopmode "NewMethod-diff-$style.tex" > /dev/null || true
  done
  echo "NewMethod-diff-$style.pdf: $(grep -ac '^!' NewMethod-diff-$style.log) LaTeX error(s)"
done
