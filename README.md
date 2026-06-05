# Computation scripts for the hydrogen example

Scripts that perform the hydrogen-atom computation of the paper (reduce the
Schrödinger PDE modulo the ODE-tower ansatz, project onto the constants,
prime-decompose). Brent Baccala, for the Journal of Computational Algebra.

Built on **François Boulier's `DifferentialAlgebra`** (BLAD/BMI) for the
differential reduction and **Sage/Singular** for the prime decomposition.

## Environment

- `joca.sage`, `joca-rg.sage`, `joca-rg-orderly.sage` — SageMath. Verified with
  **SageMath 10.7 + DifferentialAlgebra 5.3** (`~/miniforge3/envs/sage`;
  install with `sage -pip install DifferentialAlgebra`).
- `rg_basefield.py` — pure Python (sympy + DifferentialAlgebra, no Sage); runs
  under any Python with the package (e.g. `~/diffalg-venv`).

## The scripts

### `joca.sage` — the paper's computation (≈5 s)
The canonical reproduction. Ritt-reduces the PDE modulo the ansatz
(`differential_prem`), projects the remainder onto the constants, and computes
the **five** minimal associated primes: the classical (spherical) solution, the
new `E = 0` (parabolic) solution, and three discarded degenerate components
(`v = 0`; all ODE coefficients zero; first-order-ODE locus). This is
Algorithm 1 — the projection — and its output matches the paper's printed
decomposition. It never calls Rosenfeld–Gröbner.

### `joca-rg.sage` — joca.sage + the RG regularization step (≈26 s)
Adds the general algorithm's step that discharges hypothesis (3): regularize
the ansatz with `RosenfeldGroebner` before reducing. The constants are moved
into the coefficient field `ℚ(E,v₁,…,c₁)` via `BaseFieldExtension` — **required**
for RG to terminate (see the cliff below). RG returns a **single** regular
component (matching the paper's "B = ∅, one component" claim).

**It is not a no-op.** RG rewrites the ansatz — it eliminates `v` (substituting
`v = v₁x+…+v₄r`), rationalizes the radius (`r² = x²+y²+z²`), and squares the ODE
initial `(a₀+a₁v) → (a₀+a₁·…)²`. The final decomposition therefore coarsens from
five primes to **four**: the two physical solutions (classical, new) are
preserved exactly, but the two *discarded* degenerate components collapse into
`(a₁, a₀)`. This is the §1.12 effect — base-field RG is generic and cannot see
the special/degenerate strata finely — made concrete. The physics is unchanged;
the printed-decomposition match is not. Hence joca.sage (projection-only) stays
the faithful reproduction.

### `rg_basefield.py` — minimal RG-on-ansatz demo (≈16 s)
Standalone (no Sage): builds the ring, moves the constants into the base field,
and runs `RosenfeldGroebner(ansatz, basefield=F)`, showing it terminates and
returns one regular component of 9 equations.

### `joca-rg-orderly.sage` — cliff-investigation artifact (times out)
RG on the ansatz with the constants left as **ring variables** and an *orderly*
ranking. Times out — demonstrating that the cliff is ranking-independent (it is
not the elimination-vs-orderly distinction).

## The performance findings

- With the constants as **ring variables**, `RosenfeldGroebner(ansatz)` does
  **not terminate** (>28 s, killed) under *either* a block-elimination ranking
  (joca's) or an orderly ranking.
- A gdb backtrace of the hung process showed it spinning entirely in **BLAD C
  code** — `bad_remainder_irreducible_factorwise` (recursive) → `baz_Yun` /
  `baz_gcd` polynomial factorization → term sorting — with **no** Python/sympy
  or string-marshalling on the stack. So the cliff is the RG *computation*
  (factorwise reduction + factorization), not the back-and-forth string
  conversion between libraries.
- Moving the constants into the **coefficient field** `ℚ(c)` makes RG terminate
  (≈16 s, one component): constant-only coefficients become units (no
  zero-divisor splitting) and the factored polynomials lose 11 variables.

**Practical lesson:** the projection (`joca.sage`) is the fast (≈5 s) and
faithful path; the RG step is heavy and, when used, needs the constants in the
base field — and even then it regularizes (rewrites) rather than preserves the
ansatz.

## Provenance

`joca.sage` is Brent's; the `joca-rg*`, `rg_basefield.py` variants and this
README were produced while investigating the RG performance behavior with
Claude (Opus 4.8).

## A.I. use — record of AI contributions to the paper

*Working record, kept so the required journal statement on generative-AI use
can be written from fact rather than memory. Reconstructed from this
repository's git history; covers contributions through 5 June 2026 and should
be updated as work continues. Commit hashes are given for auditability.*

The AI tool is **Anthropic's Claude**, used both interactively and through the
**Claude Code** CLI (which authors the `Claude Code` git commits; recent
sessions ran models in the Claude Opus 4 family). A few citation suggestions
came from **OpenAI's GPT-5** (commit `55487de`). **All AI output was reviewed
and edited by the author (B. Baccala), who takes full responsibility for the
content. AI is not, and cannot be, an author.**

**The author's own (not AI).** The method itself — the differential-algebraic
ansatz approach to solving PDEs — the research program, the core mathematical
ideas, the hydrogen ansatz and the new solution, the choice of results, and the
bulk of the original prose, figures, and worked examples (the development from
late 2025 onward; the version sent to the human collaborator "Paul" in
January 2026).

**Substantive / mathematical assistance.**
- Drafted the **completeness theorem and its proof** (the "claude-generated
  theorem and proof", `843f9ea`; inlined/edited `874a9c7`) — a formalization of
  a result about the author's method, reviewed and checked by the author.
  Follow-on: promoted the specialization result to a lemma and renumbered
  (`57ce904`), stated hypothesis (3) generically (`67aeb0f`), verified the four
  hypotheses in the worked example (`cc6ee77`, `7c18e49`).
- **Restructured** the logical organization: core + general algorithm
  (`e2a7449`), then folded into one section with the "Why not
  Rosenfeld--Gröbner alone?" discussion (`c05c04a`); unified the
  theorem/algorithm/example numbering (`3c0818c`).
- **Literature search / prior-art identification** that materially shaped the
  claims: found Fakouri–Rahmany–Basiri (2018) and the differential Thomas
  decomposition and recast the comprehensive step from an original claim into
  cited prior art (`ff291ad`); positioned the method against the
  symmetry-reduction tradition (Lie; Clarkson–Kruskal; Bluman–Cole) in the
  Conclusion (`102fd6a`); identified the Chaharbashloo–Basiri–Rahmany(–
  Zarrinkamar) "ansatz + Gröbner in quantum mechanics" prior art.
- **Corrected a mathematical claim**: caught that the new hydrogen solution,
  asserted non-separable, is in fact separable in parabolic coordinates, and
  rewrote the passage with the Laplace–Runge–Lenz / SO(4) framing, a structural
  (non-`L²`) explanation, and supporting citations (`806b526`), checked against
  Landau–Lifshitz §37.
- **Computational verification** that informed the text: ran the author's Sage
  computation and profiled the Rosenfeld–Gröbner step (gdb backtrace into BLAD;
  the base-field fix), grounding the efficiency discussion and the corrected
  "RG is not a no-op" wording (`315ab7c`). See the scripts above.

**Editorial / typesetting assistance.** Reformatted the bibliography to
Elsevier's numbered style (`2b71ff5`) and added/repaired citations (`02b57b8`
and others above); tightened and reconciled prose in the Projection and
completeness subsections (`aba1b6e`, `0aee721`, `6e1ee75`, `52664d7`); fixed
LaTeX warnings (`226a2a5`).

**Artifacts.** The **Graphical Abstract** is AI-generated (`6db3f5b`); the
supporting computation scripts (`joca-rg.sage`, `rg_basefield.py`) and this
README are AI-produced.

**Note for the formal declaration.** Elsevier's generative-AI policy is framed
around assistance with *language and readability*. Several items above
(theorem/proof drafting, restructuring, literature search, the graphical
abstract) go beyond language polishing; the author should decide how to
characterize these in the submitted statement, methods, and/or
acknowledgements, and confirm the journal's current requirements.

*/Claude Opus 4.8 (this A.I.-use section drafted by the AI it documents).*
