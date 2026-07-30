# A.I. use — record of AI contributions to the *NewMethod* paper

*Working record, kept so the required journal statement on generative-AI use
can be written from fact rather than memory. Reconstructed from this
repository's git history; covers contributions through 27 July 2026 and should
be updated as work continues. Commit hashes are given for auditability.*

The AI tool is **Anthropic's Claude**, used both interactively and through the
**Claude Code** CLI (which authors the `Claude Code` git commits; sessions
through mid-2026 ran models in the Claude Opus 4 family, later ones Claude
Opus 5). A few citation suggestions
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
- **Repaired the containment proofs in "The Constant Loci"** (author's commit
  `92e77e5`, 27 July 2026; model Claude Opus 5). Found that the claimed
  containment V_∀ ⊆ V_{∃\Ψ} is **false** as stated — consistency of the ansatz
  is too weak a hypothesis. Supplied the missing condition (that the ansatz not
  force Ψ to vanish, 1 ∉ [A(c*)]:Ψ^∞), the two-line proof (equal ideals have
  equal saturations, so P ∈ [A(c*)] gives [A(c*),P]:Ψ^∞ = [A(c*)]:Ψ^∞), and a
  counterexample family A(c) = {Ψ″−Ψ, Ψ′−cΨ}, in which (c²−1)Ψ ∈ [A(c)] and so
  V_∀ = ℂ∖{−1} while V_{∃\Ψ} = {1}. Observed that the draft's homogeneity step
  (Ψ ∈ [A(c*)] ⟹ P ∈ [A(c*)]) is what *generates* such counterexamples rather
  than a step toward the theorem, and that the repaired proof needs no
  homogeneity hypothesis at all — homogeneity is required only upstream, to make
  V_{∃\Ψ} the object of interest. Also caught a ∨-for-∧ error in the V_∀ ⊆ V_∃
  proof. The author incorporated these; flagged and still open at the time of
  writing are a quantifier gap between the pointwise hypothesis and the global
  containment, the now-redundant weaker consistency assumption, and the missing
  prose justification for the third proof.
- **Located the differential-algebra definition of homogeneity, and the correct
  grading for the ansatz** (27 July 2026; model Claude Opus 5). Text-extracted
  and searched 201 PDFs (the `~/project/papers` collection and the author's
  *Books/Differential Algebra*), plus DjVu OCR of the classical texts, for a
  differential-algebra definition of "homogeneous differential equation".
  **Robertz, LNM 2121, Def. A.3.13 / Lemma A.3.14 / Remark A.3.16** is the only
  genuine one — homogeneous differential ideal with respect to the standard
  grading, the equivalence to closure of the solution set under scalar
  multiplication, and the caveat that homogeneity is not invariant under
  coordinate changes of the dependent variables; van der Put–Singer define the
  term for scalar *linear* equations only (b = 0). Ritt (1932, 1950), Kolchin
  *Differential Algebraic Groups*, Hubert I/II, Boulier–Lazard–Ollivier–Petitot,
  Bächler et al., the BLAD manual, Fakouri and Lange-Hegermann use "linear
  homogeneous" as an undefined adjective or never (Kolchin *DAAG* and Kaplansky
  are image-only scans and were not checked). Identified that the paper's ansatz
  ideals *are* homogeneous under the **Ψ-degree grading** — auxiliary
  indeterminates v, x, y, z, r in degree 0, a grading the derivations preserve —
  resolving the author's concern that the ansatz ideal need not be homogeneous,
  and flagged that the hypothesis "if P is homogenous (and the ansatz too)" is
  false under the paper's own definition of homogeneity, since the generator
  defining the independent variable is Ψ-free.
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
- **Drafted the 6a/6b split of the core algorithm's assembly step** (July 2026),
  at the author's request: worked out that replacing the intersection by a union
  of the per-component varieties, each cut down to its own cell, yields an
  intermediate locus V_{∃∀} with V_∀ ⊆ V_{∃∀} ⊆ V_∃. Supplied the quantifier
  table distinguishing the three loci (∀/∀, ∃/∀, ∃/∃ over components and over
  the solutions within a component), both containment arguments, the observation
  that the paper's own three-loci example makes the upper containment maximally
  strict (V_{∃∀} = ∅ against V_∃ = ℂ), the reason no rearrangement of the
  assembly step can compute V_∃ — the ∀-projection has already discarded the
  needed information — and the Ψ-reduction filter that places the intermediate
  locus inside V_{∃∖Ψ} in the homogeneous case. Also noted that, unlike V_∀ and
  V_∃, the intermediate locus is not intrinsic: it depends on which
  decomposition was computed.
- **Computational verification** that informed the text: ran the author's Sage
  computation and profiled the Rosenfeld–Gröbner step (gdb backtrace into BLAD;
  the base-field fix), grounding the efficiency discussion and the corrected
  "RG is not a no-op" wording (`315ab7c`). See the computation scripts
  documented in `README.md`.

**Editorial / typesetting assistance.** Reformatted the bibliography to
Elsevier's numbered style (`2b71ff5`) and added/repaired citations (`02b57b8`
and others); tightened and reconciled prose in the Projection and completeness
subsections (`aba1b6e`, `0aee721`, `6e1ee75`, `52664d7`); fixed LaTeX warnings
(`226a2a5`). Line-level review of the Constant Loci subsection (July 2026):
notation collisions (a set-builder colon abutting an ideal saturation; a
differential ring written with the independent variable and the differential
indeterminate interchanged), an independent variable named inconsistently
between the algebraic and analytic readings of the same example, an unresolved
author query left in the body text, a locus asserted in the summary line whose
announced computation is never carried out, and the paper-wide
"homogenous"/"homogeneous" spelling split.

**Artifacts.** The **Graphical Abstract** is AI-generated (`6db3f5b`); the
supporting computation scripts (`joca-rg.sage`, `rg_basefield.py`), the
directory's `README.md`, and this record itself are AI-produced.

**Note for the formal declaration.** Elsevier's generative-AI policy is framed
around assistance with *language and readability*. Several items above
(theorem/proof drafting, restructuring, literature search, the graphical
abstract) go beyond language polishing; the author should decide how to
characterize these in the submitted statement, methods, and/or
acknowledgements, and confirm the journal's current requirements.

*/Claude Opus 4.8 (this record drafted by the AI it documents); updated
27 July 2026 by /Claude Opus 5.*
