# A.I. use — record of AI contributions to the *NewMethod* paper

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
- Drafted the **companion note `gcrd-closure-partial-strata.tex`**
  (`3bd3828`, 7 July 2026; model Claude Fable 5): given a linear-ODE-tower
  ansatz, a construction — GCRD of Ore polynomials plus the
  Chardin–Li subresultant stratification — of a finite set of lower-order
  ans\"atze whose membership loci capture every partial-solution stratum, so
  the membership (∀) algorithm run on the collection matches the completeness
  of the existential formulation. The note grew out of an interactive
  discussion with the author (Riquier–Thomas existence ⇒ ∃-reading of
  "solves"; the ∀/∃ cost gap; the observation that partial strata require
  solution-space dimension ≥ 2). Single-linear-element case proved;
  tower-level uniform degree bound and higher-order nonlinear elements
  flagged as open. Unreviewed draft pending the author's check.
- Drafted the **companion note `ansatz-method-provenance.md`** (`97f718c`,
  18 July 2026; model Claude Opus 4.8): a verified four-layer provenance of the
  "ansatz method" keyword — generic term → **Dong 2002** as the term-of-art
  coiner (Phys. Scr. 65:289; Found. Phys. Lett. 15:385) → the **Flessas 1979**
  polynomial×exponential trial form inside the **quasi-exactly-solvable (QES)**
  program (Turbiner; Bender–Dunne; Ushveridze) → the **QES hidden Lie algebra**
  (sl(2,ℝ) gauge construction; classified and lifted to two-variable PDE by
  González-López–Kamran–Olver, J. Phys. A 24:3995 1991 and CMP 159:503 1994).
  Grew out of an interactive session that stress-tested the paper's related-work
  line against the strongest prior art: the naive "prior work only does ODEs" is
  false (QES reaches PDEs via separation and via the hidden algebra), so the note
  reframes the differentiator as **precondition-free differential elimination +
  the completeness certificate** (applicable to non-separable, no-hidden-sl(2)
  problems like helium). Citations CrossRef-verified; a mis-citation of DOI
  BF02099982 (it is CMP 159:503 1994, not CMP 153:117 1993) and a Magyari-vs-
  Flessas mis-attribution were caught. Also ran a fan-out sibling/citation-forward
  search confirming the ansatz+Gröbner coupling is nearly unique to Chaharbashloo
  (one verified sibling, Acosta-Humánez–Venegas-Gómez 2018; zero verified forward
  citations). Draft pending the author's review.
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
  "RG is not a no-op" wording (`315ab7c`). See the computation scripts
  documented in `README.md`.

**Editorial / typesetting assistance.** Reformatted the bibliography to
Elsevier's numbered style (`2b71ff5`) and added/repaired citations (`02b57b8`
and others); tightened and reconciled prose in the Projection and completeness
subsections (`aba1b6e`, `0aee721`, `6e1ee75`, `52664d7`); fixed LaTeX warnings
(`226a2a5`).

**Artifacts.** The **Graphical Abstract** is AI-generated (`6db3f5b`); the
supporting computation scripts (`joca-rg.sage`, `rg_basefield.py`), the
directory's `README.md`, and this record itself are AI-produced.

**Note for the formal declaration.** Elsevier's generative-AI policy is framed
around assistance with *language and readability*. Several items above
(theorem/proof drafting, restructuring, literature search, the graphical
abstract) go beyond language polishing; the author should decide how to
characterize these in the submitted statement, methods, and/or
acknowledgements, and confirm the journal's current requirements.

*/Claude Opus 4.8 (this record drafted by the AI it documents).*
