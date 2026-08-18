# A.I. use — record of AI contributions to the *NewMethod* paper

*Working record, kept so the required journal statement on generative-AI use
can be written from fact rather than memory. Reconstructed from two sources:
this repository's git history, which says what changed and when, and the
session transcripts, task reports and compactions kept alongside it, which say
whose idea it was. The git history alone cannot settle the second question, and
it is not sufficient even for the first — the author committed AI-drafted
material under his own name at least six times, and the AI committed the
author's proposals under its name many more. Every claim below of the form
"proposed by the author" or "drafted by the AI" rests on the transcripts, not
on the commit's author field. Covers contributions through 18 August 2026 and
should be updated as work continues. Commit hashes are given for auditability,
and the appendix tabulates every commit to the paper or its companion notes
that the AI touched.*

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
  `53c4eba`, 27 July 2026; model Claude Opus 5). Found that the claimed
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

**Companion notes and working documents.** Twelve `.tex`/`.md` documents sit
beside `NewMethod.tex` in this directory. **All twelve were drafted by the AI**;
none is part of the submitted paper, and their status ranges from checked and
absorbed into the paper to unreviewed draft. They are recorded here because
material has repeatedly moved from a note into the paper, and because several of
them are the source of statements the paper now makes. Two are already described
above (`gcrd-closure-partial-strata.tex` and `ansatz-method-provenance.md`); the
remaining ten are:

- **`rg-saturation-and-the-bad-locus.tex`** (`59ccab9`, `b659ce8`, 5 June 2026;
  model Claude Opus 4.8). Found that `joca-rg.sage` was computing a *wrong*
  decomposition — 4 primes with a spurious `(a₀,a₁)` — not merely a coarse one:
  a regular differential chain represents the *saturated* ideal [Ā]:H̄^∞, and
  the script was reducing against the chain's bare equations, dropping the
  saturation and re-admitting the bad locus as a component. The note writes up
  the mechanism, including that `normal_form` returns the reduction as a
  fraction whose denominator (a₀+a₁v)² vanishes exactly on the bad locus, so
  the bad locus is a *pole* of the saturated reduction — which is why the
  regularized chain cannot see the two genuine strata living there. This note's
  identification of (a₀,a₁) as "the bad locus B where an initial vanishes" is
  what later exposed a false claim in the paper itself (below).
- **`ck-direct-method-examples.tex`** (`243b96d`, 8 June 2026 — swept into the
  repository by an unrelated commit, hence the misleading message). A worked
  Clarkson–Kruskal direct-method companion (heat kernel and the hydrogen J₀),
  written because the 1989 original is paywalled and Clarkson–Mansfield 1994
  was substituted as the citable source (`cc09d29`).
- **`method-comparison.md`** (`3d39b0a`, `65fe337`, 11–12 June 2026). A 525-line
  comparison of the method against the constructible-set / comprehensive-Gröbner
  / GTZ toolchain, parametric Rosenfeld–Gröbner, and the differential Thomas
  decomposition, later extended with a worked Δ(Ψₓ,Ψ_y) critical-pair split as a
  concrete atlas-versus-image example. Not cited by the paper, but the survey it
  contains is upstream of the paper's eventual decision to lead with the Thomas
  decomposition.
- **`membership-vs-variety-partial-strata.tex`** (`901ce6a`, 18 June 2026;
  extended `72fe9f9`, 9 July, and `f3ba3f5`, 23 July). The note that names the
  **partial-solution strata** — the gap between the membership locus V_∀ and the
  projection of the joint variety — and so supplies the vocabulary the paper's
  three-loci discussion now uses. It grew out of an interactive session in which
  the distinction was worked out jointly with the author. Later additions: a PDE
  counterexample showing that autoreducedness does not preclude strata, and an
  explicit GCRD/Gröbner verification that ansatz 5 on the hydrogen equation has
  no genuine (h ≠ 0) partial strata at all — its whole V_∀-versus-V_{∃∖Ψ} gap
  lies on the ODE-degeneration loci the decomposition already isolates.
- **`bounded-prolongation-theorem.md`** (`5d56eb6`, 2 July 2026, 871 lines;
  corrected `bc7a8c4`, 5 July). Proposes replacing reduce-then-project by
  prolong to a finite involutive order, algebraize the jets, never saturate,
  eliminate to the constants. A verified deep-research pass (23 of 25 claims
  confirmed by adversarial re-checking) returned a **negative verdict on
  novelty**, and the note says so: the core pipeline is essentially
  Lange-Hegermann 2014 over the differential Thomas decomposition, and the
  prolong-then-algebraize construction appears verbatim in
  Golubitsky–Kondratieva–Ovchinnikov 2009. What survives as contribution is the
  explicit N₀ = 2 for the finite-type hydrogen ansatz, the saturation-free
  reading, and the algorithm-agnostic packaging. Two errors of the AI's were
  caught afterwards, both prompted by the author's questions: the elimination
  J ∩ K[c] is identically (0), because the trivial solution Ψ ≡ 0 lies over
  every c, and the correct form is (J : Ψ^∞) ∩ K[c] (`bc7a8c4`); and the N₀ = 2
  claim is false for the *combined* system A ∪ {P}, being the involutive order
  of the ansatz alone. This note is unreviewed and its route was not adopted by
  the paper.
- **`riquier-parametric-data.tex`** (`6f2d20a` as markdown, converted `c9f682a`,
  7 July 2026). Cauchy–Kovalevskaya, the Riquier existence theorem, and the
  ansatz-5 parametric data — written as background during a tutorial thread with
  the author that later converged on the GCRD-closure construction.
- **`prolongation-projection-algorithm.tex`** (`b8e2cec`, 10 July 2026, 995
  lines; extended `95d687c`, 13 July, 553 lines). **The algorithm developed here
  is the author's proposal** — prolong, track, split, project; the AI supplied
  the development, the semantics and the resolution of the puzzle the author had
  posed about the ansatz element v = v₁x + v₂y + v₃z. The resolution is that the
  quantifier discipline (independents universal, dependents and constants
  existential) is semantically right but must be implemented by two operations
  that do not commute: universal quantification is coefficient collection,
  existential quantification is Gröbner elimination, and the latter is sound only
  after prolongation to passivity. The naive ordering computes a locus
  **incomparable to both targets** — over-constraining on the v-equation, and
  under-constraining on a two-line toy example where it returns ℂ against the
  correct {c = 1}. The later extension locates precisely what the paper's staged
  route (decompose the ansatz alone, reduce, project) loses: a cell hides a
  stratum iff its reduced remainder carries a parametric jet. Delivered as a
  task-runner report; not incorporated into the paper.
- **`coherent-regular-simple-passive-hierarchy.tex`** (`4063c0c`, `fd8f64d`,
  20 July 2026). Pins *coherent*, *regular*, *simple* and *passive/involutive*
  to primary sources and establishes that they are **not** synonyms but a
  one-directional hierarchy. Written after an earlier informal claim of the
  AI's — that Seiler's AAECC 2009 papers state a coherence ⟺ involution
  *equivalence* — was checked and **found false**: those papers never mention
  coherence or Rosenfeld at all, and the sourced statement is a one-way
  implication. The note is the one companion document the paper leans on
  hardest; its two later corrections (`41a7102`) and its new §5 (`747002c`) are
  narrated in the arcs below.
- **`membership-existence-equivalence.tex`** (`f3fda04`, `50cf452`, 23 July
  2026). Proves V_∀ = V_{∃∖Ψ} ⟺ no partial strata ⟺, for linear-ODE ansätze,
  ord GCRD(L, {R_β}) ∈ {0, n} everywhere; coincidence proved for one-dimensional
  ansätze, with a follow-on remark on which hypotheses are actually necessary.
- **`regularity-counterexamples.tex`** (author's commit `933ef24`, 23 July
  2026). This file is **entirely AI-drafted material**, and it is the one case in
  the repository where the author preserved AI text by moving it *out* of the
  paper rather than by keeping it in. The counterexamples came into the paper
  with the completeness proof (`874a9c7`): a demonstration that hypothesis (4)
  needs squarefreeness and not merely a regular chain — A = {u³ − u²}, whose
  separant shares the factor u, so that u − 1 is fully reduced, nonzero, and
  still lies in (A) : S_A^∞. The author removed these subsections from
  \S*Algorithm* as too detailed for a journal paper and copied them here. The
  material is not in the submitted paper; it is recorded because the commit that
  moved it is one of the six in which the author committed AI prose under his
  own name.

**A correction to the paper arising from the notes** (12 June 2026; model
Claude Opus 4.8).

- **The bad locus is not empty** (`18ddc11`, `d965dbb`). The paper asserted in
  three places that the hydrogen ansatz has an empty bad locus B, resting on the
  claim that an element linear in its leader satisfies hypothesis (4)
  automatically. That conflates *linear-in-leader* with *monic*: the ODE element
  (a₀+a₁v)Ψ″ + (b₀+b₁v)Ψ′ + (c₀+c₁v)Ψ is linear in Ψ″ but non-monic, its initial
  equals its separant a₀+a₁v, and that vanishes identically on {a₀ = a₁ = 0},
  dropping the leader. So B is the codimension-2 subspace {a₀ = a₁ = 0} ⊆ ℂ¹¹ —
  and it *contains genuine decomposition strata*. This matched what the
  rg-saturation note had independently called the bad locus, which is how it was
  caught. Three passages corrected and a vague Hubert citation sharpened to
  Theorem 5.13 / Definition 5.5.

The bullets that follow record a single connected arc, 28–30 July 2026 (model
Claude Opus 5), carried out at the author's direction: shifting the paper's
foundations from Rosenfeld–Gröbner regular systems to the differentially simple
systems of the differential Thomas decomposition, and following the consequences
through the completeness proof and the algorithm.

- **Rewrote the preliminaries on simple systems** (`7532ee5`, `68af580`).
  Filled in the stubbed block: Janet division, multiplicative and
  non-multiplicative variables, non-multiplicative prolongation,
  Janet-involutive, algebraically simple, differentially simple (S1)–(S4), and
  the Thomas decomposition itself; dropped Boulier–Lazard–Ollivier–Petitot
  Def. 22, restated the specialization lemma for a Thomas component, and
  collapsed two theorem hypotheses into one. **Introduced a bug and caught
  it**: removing BLOP Def. 22 left the proof citing BLOP Thm. 23, whose
  hypothesis is literally "regular differential system". Repaired by stating
  Rosenfeld's Lemma inline with (C1)–(C3) unfolded and bridging from
  differential simplicity — (C3) from (S2), triangularity from (S1), (C2) from
  square-freeness via Bächler et al. footnote 3 — which the next commit then
  made moot.
- **Established that the "Janet-involutive ⟹ coherent" step the paper was
  leaning on has no published source, and found the bypass** (`d015a61`).
  Literature search: Bächler et al. §1's "hence they are coherent" has nothing
  behind it, Chen and Gao (2003) list the relation as an open problem, and
  Bächler et al. never attach an ideal to a simple system at all, so the ideal
  theory could not come from there either. The replacement is **Robertz,
  LNM 2121, Prop. 2.2.50** — one numbered equivalence whose only hypothesis is
  simplicity, delivering what coherence, Rosenfeld and Lazard were being
  assembled to deliver. Verified in the source rather than taken from the
  search report, which turned up **Prop. 2.2.72** alongside it — the
  decomposition identity √(E:q^∞) = ⋂_i (E^(i):(q^(i))^∞) — closing a gap that
  had been flagged as unsourced. Net −54 lines: critical pairs, Δ-polynomials,
  coherence, Rosenfeld's Lemma and the Hubert regular-chain block all left the
  preliminaries. Also answered the author's objection that Prop. 2.2.50 says
  nothing about the pseudo-remainder lying in the *algebraic* ideal: Ritt full
  reduction is partial reduction followed by algebraic reduction, the old proof
  stopped at the seam, and Prop. 2.2.50 states the composite (Li–Wang Thm. 8
  confirms by stopping at the same seam explicitly).
- **Read Li–Wang from page images and corrected the companion hierarchy note**
  (`41a7102`, `coherent-regular-simple-passive-hierarchy.tex`). The local PDF is
  an image-only scan; rendering it page by page established that Li–Wang **never
  define "passive"** — §4.1 defers to Wu (1989) — so their Thm. 3 is about
  Wu-passivity, not Janet-involutivity, and the note's implication lattice had
  been using it to justify a Janet-involutive arrow. Split into two arrows.
  Second correction to the same note: it described 𝕋 as a characteristic set of
  a **prime** differential ideal, where Li–Wang Cor. 2 says **radical** — a
  second-hand claim the note had carried unchecked.
- **Corrected the Robertz 2018 citation and the proposition numbers**
  (`93a031b`, `6617bdd`). It is a **published article** — *Les cours du CIRM*
  vol. 6 no. 1, Course no. III, pp. 1–37, doi:10.5802/ccirm.28 — not lecture
  notes; the English title in the PDF metadata is a leftover `\hypersetup` from
  the 2014 book. The published version inserts two definitions, so **the
  proposition numbers shift by two** from the preprint the paper had been citing
  (3.29/3.30 → 3.31/3.32); verified against the published PDF, which was
  downloaded and added to the local collection. Every Robertz citation site now
  carries both references, at the author's request, and the proposition was
  numbered into the paper's shared theorem sequence.
- **Dropped the denominator hypothesis** (`3d92072`). h(c*) ≠ 0 follows from c*
  lying in the cell, once "cell" is taken to mean the projection of the
  component's solution set — with the subtlety that h involves jet variables, so
  what is actually needed is h ≢ 0 as a polynomial. Algebraic simplicity written
  out as (A1)–(A3) in the same enumerated style as (S1)–(S4), at the author's
  request.
- **Reformulated the membership locus on the radical, and made the theorem an
  equivalence** (`d9c4c8d`). Proposed by the author (that the argument yields
  P ∈ √[A(c*)] for free) and sharpened here: the exact object is [A(c*)]:q^∞,
  which contains the radical, and every step is reversible, so
  c* ∈ 𝒱 ⟺ P ∈ [A(c*)]:q^∞ and the sufficiency disclaimer could go. Four
  coupled changes followed — V_∀ redefined on the radical, the standing
  assumption simplified to the transparent Ψ ∉ √[A(c*)], the ∃-projection
  dropped from the algorithm (kept on record as a cheaper, complete but
  inexact variant), and a new corollary reading Prop. 2.2.72 with q = 1, so that
  the per-component test decides the saturation while the intersection over
  components decides the radical.
- **Recorded how a component's cell is obtained** (`83a4465`) — read off the
  equations and inequations whose leaders are constants, justified by Bächler
  et al. Remark 2.3's extension property, with a warning against confusing this
  with the algorithm's projection step (applying the latter to the ansatz's own
  equations would return the empty cell). Then **corrected the algorithm's
  assembly** (`26ff25e`) from a union of components to
  ⋂_i (𝒱_i ∪ (ℂⁿ∖C_i)), the cells being quasi-affine and overlapping rather
  than partitioning.
- **Added a second, cheaper assembly step** (`1380cf6`, 30 July 2026), at the
  author's request: alongside step 6a (the intersection, which returns the
  membership locus exactly) a step 6b taking the union of the per-component
  varieties cut down to their own cells, yielding an intermediate locus V_{∃∀}
  with V_∀ ⊆ V_{∃∀} ⊆ V_∃. Supplied the quantifier table distinguishing the
  three loci (∀/∀, ∃/∀, ∃/∃, over components and over the solutions within a
  component), both containment arguments, the observation that the paper's own
  three-loci example makes the upper containment maximally strict (V_{∃∀} = ∅
  against V_∃ = ℂ), the reason no rearrangement of the assembly step can compute
  V_∃ — the ∀-projection has already discarded the needed information — and the
  Ψ-reduction filter that places the intermediate locus inside V_{∃∖Ψ} in the
  homogeneous case. Also noted that, unlike V_∀ and V_∃, the intermediate locus
  is not intrinsic: it depends on which decomposition was computed.
  **Restructured the algorithm** in a follow-on (`232cec8`), prompted by the
  author's question of how the new branch is actually computed: since
  𝒱_i ∪ (ℂⁿ∖C_i) = W_i ∪ (ℂⁿ∖C_i), both assemblies read 𝒱_i only through
  W_i = 𝒱_i ∩ C_i, so the prime decomposition was hoisted out of a trailing
  step and into the per-component loop, and the two assembly formulas became
  steps 7a and 7b. Named the ideals the loop actually builds (the coefficient
  ideal J_i, the cell's equation ideal E_i and inequation product h_i, the
  saturation (J_i + E_i):h_i^∞), and supplied three corrections: saturation
  returns the Zariski *closure*, so each surviving prime still carries h_i ≠ 0
  as a residual condition and it is this — not anything in the assembly — that
  makes the output constructible; the redundancy test needs the inequations
  (prime containment alone can fail when V(𝔭) lies inside V(h_j)); and the
  radical is never formed, an ideal and its radical having the same minimal
  primes, so the paragraph instructing the reader to take the radical first was
  wrong as an instruction and was rewritten as motivation. Finally supplied the
  missing recipe for step 7a (`07ab265`), which until then was a formula with no
  method behind it: writing each factor as a union of locally closed pieces
  (W_i, one open set per generator of E_i, and V(h_i)), the rule that two pieces
  meet in V(𝔞+𝔟) ∖ V(fg), the emptiness-and-discard test (the saturation being
  the unit ideal), and the observation that the combinatorial blow-up over
  choices of one piece per factor — not the ideal arithmetic — is what makes 7a
  the expensive branch, with two mitigations that keep it in hand. That
  hand-rolled recipe was then **replaced by the published one** (`ab1166c`)
  after a literature search established that set-theoretic computation with
  constructible sets is solved and implemented: the right operation is a
  disjoint refinement of the family W_1, C_1, …, W_r, C_r, after which each
  block lies inside or outside every W_i and C_i and both assemblies are a
  selection of blocks — irredundant by construction, so the containment test
  drafted for 7b became unnecessary and was removed. Identified the refinement
  as the set-theoretic coprime factorization problem of the already-cited CTD
  paper (whose difference algorithm works on regular systems — a regular chain
  plus an inequation, the same shape as step 6's output), its implementation as
  MakePairwiseDisjoint in ConstructibleSetTools, the algebraic Thomas
  decomposition as the cheaper route here since it is step 1's algorithm in the
  non-differential case and returns disjoint output already, and Brunat–Montes
  for a canonical rather than merely irredundant representation. Two references
  added and the previously uncited Montes bibitem put to use. One citation trap
  caught: search engines report the companion manuscript "Computing with
  Constructible Sets in Maple" as J. Symbolic Comput. 45(1) 124–149; those pages
  are in fact Boulier–Lemaire–Moreno Maza on differential characteristic sets,
  and the manuscript appears never to have been published, so the ICCSA 2008
  proceedings version is cited instead. Finally **reformatted the algorithm to
  the conventions of the papers it cites** (`4afff4d`) — Bächler et al. (JSC
  2012) and Fakouri et al. (2018): a named algorithm (MembershipLocus), Input
  and Output as pre- and postconditions rather than prose, numbered pseudocode
  lines with assignment arrows and explicit control flow, named subalgorithms
  called and cited the way those papers call their own, and the 7a/7b branch
  recast as a mode flag with two returns. Caught in passing a pre-existing
  notation collision — r was both the number of components and the remainder
  produced by reduction — and renamed the component count to s. Answering the
  author's question of whether step 1 can decide the standing consistency
  assumption, established that it can and that the assumption was therefore
  unnecessary as a precondition (`b585d77`): inconsistency cannot surface as a
  defective component, since every simple system has a solution (ThomasDecomp
  Remark 2.3) and the decomposition discards a branch the moment it produces
  c_= for a non-zero constant or 0_≠ — it surfaces as a gap in the covering by
  cells, the cells being projections of a partition of Sol(A(c)). Refining ℂⁿ
  along with the W_i and C_i then makes consistency a by-product, and the ∀
  assembly needs no hypothesis at all, a block lying in no cell satisfying its
  condition vacuously and being returned — which is correct, the radical being
  the unit ideal there. The ∃∀ assembly still needs it, as does the sandwich,
  since at an inconsistent c* the point lies in V_∀ but not V_∃. Consequently
  dropped the same hypothesis from the assembly corollary (`f4996bb`), the
  intersection over an empty index set being the whole ring and hence exactly
  the radical there; and, at the author's suggestion, moved the refinement
  inside the branch, observing further that the ∃∀ condition mentions neither
  ℂⁿ nor the cells, so that branch refines s sets rather than 2s+1. In the same
  pass **corrected my own Output specification**: a refinement block is a
  Boolean combination of the sets refined and an intersection of irreducible
  sets need not be irreducible, so the promised presentation by primes required
  a final normalization line, which had been missing. On the author's objection
  that an algorithm cannot assign "the ideal generated by …", an infinite
  object, reworked the algorithm to carry finite generating sets throughout
  (`3670235`), with the representation convention stated once before it and the
  saturations, the Output and the surrounding remarks adjusted to match. Added a
  table of the four subroutines the algorithm calls, with what each computes and
  its source (`bf45291`), cited by name rather than by number after finding that
  the local copy of the Thomas-decomposition paper is the arXiv preprint, not
  the published JSC version the bibliography names — the same preprint/published
  gap that shifted the Robertz proposition numbers.
- **Surveyed the differential-algebra literature for standard ring notation and
  migrated the paper to it** (July 2026; model Claude Opus 5, research delegated
  to a subagent over ~30 sources — the `~/project/papers` collection, an OCR of
  Ritt's first two chapters, the BLAD/BMI manuals and the live
  `DifferentialAlgebra` source; report at
  `~/project/reports/differential-algebra-ring-notation-survey.md`). The verdict
  was unanimous: braces hold differential indeterminates, and no source writes
  ℚ{x,y,z} for independent variables. Established that the paper's convention,
  attributed to Ritt, is not Ritt's — §3 defines 𝔉{y₁,…,yₙ} over dependent
  variables and §29 adjoins the independent variable, with z′ = 1, to the
  *coefficient field*, which is the author's own "associated indeterminate" idea
  on the other side of the braces. Anchor citation: Boulier's own docstring in
  `DifferentialRing.pyx`, "Following Ritt and Kolchin notation, R = Q[x,a,b]{z,y}
  with derivation d/dx", matching what `joca.sage` already does — only the prose
  diverged. Migrated 17 sites (`3efd8bd`), 11 inside tikz figures, and fixed
  three further defects the survey exposed: the paper's only full statement of
  its own convention sat inside a `\begin{comment}` block and was never typeset;
  ℂ(c){t}[Ψ] inverted the parameter c, forbidding the c\* = 0 specialisations the
  method depends on; and the quotient rings /(r²−x²−y²−z²) contradicted the
  paper's own text two pages earlier. The author supplied the decisive
  correction to the recommendation as first drafted: it had put the constants in
  the coefficient ring, where they are inert, when the algorithm needs them
  ranked in the lowest block so that the decomposition splits on them — which is
  the entire source of the cells.
- **Audited the algorithm's own citations, and found six defects — five of them
  mine** (30–31 July 2026). Each was caught by the author asking whether a claim
  was actually true, not by any verification pass of mine, and they are recorded
  individually because the pattern matters: the newly-formalized algorithm cited
  more precisely than it had checked. (i) `DisjointRefinement` was a subroutine
  name **I invented** and attributed to the already-cited CTD paper; renamed to
  `SMPD`, the paper's actual Algorithm 3 (`8de5053`). (ii) The prose credited
  `MakePairwiseDisjoint`, which is the *weaker* operation and does not appear in
  the ICCSA paper cited for it; rewritten to cite each source for what it
  documents. (iii) `Reduce` was cited as Bächler et al. Alg. 2.6 — the
  *algebraic* Reduce; the differential one is Alg. 3.3 (`260d919`). (iv) That
  correction was then superseded on the author's question of whether Ritt's
  reduction could simply be used throughout: the paper's introduction,
  preliminaries, projection section, hydrogen example and the completeness
  proof all use **Ritt full reduction**, and Proposition 1's criterion is a
  pseudo-remainder rather than Janet-reducibility, so line 3 became
  `FullReduce`, cited to Ritt §6 (`2553068`, `631e6e9`). (v) `minAss` was cited
  to GTZ as a whole; narrowed to Cor. 3.2(v) and §9 (`d35e9cd`). (vi) Numbered
  subroutine citations and the Example 2 notation were reconciled with the
  algorithm's lines 4–6 (`dd996e2`). A seventh defect was **not** mine and
  predates the arc: the cell-reading paragraph rested on the premise that the
  constants are ranked below every other indeterminate, whereas the paper's own
  hydrogen ranking display has read "… > constants > x,y,z" all along, so prose
  and ranking had disagreed for some time; the algorithm's lines 5–6 now select
  membership in ℚ[c₁,…,cₙ] directly rather than relying on the leader.
- **Notation and exposition follow-ons** (30–31 July 2026). Reserved primes for
  the derivatives internal to an ODE extension and wrote Δ-derivatives as
  subscripts throughout, the two having been used interchangeably (`0e34162`);
  showed r in the base ring of the first two ansatz diagrams, which had omitted
  it (`47e9360`); **rewrote the Organization paragraph** (`a59b8ae`), stale in
  three ways, the serious one being that it promised a general algorithm built
  on Rosenfeld–Gröbner that the paper does not contain — the one surviving live
  mention of Rosenfeld–Gröbner names it as computationally infeasible and
  therefore not taken; moved the constants-as-a-block rationale from §1 to sit
  beside the cells it explains, where it is load-bearing rather than
  unmotivated (`0b3c482`, `e953bc8`); and **reordered the bibliography by first
  citation** per Elsevier style (`125b8ae`), with the order computed from the
  *typeset* text only — `\begin{comment}` blocks and `%`-comments stripped,
  since citations invisible in the PDF must not influence numbering — and
  verified entry by entry for all 25 cited entries. Flagged and left to the
  author: six bibliography entries are uncited, which Elsevier will not accept.
- **Errors of mine during this arc, and how they were caught.** Twice I
  over-claimed and had to withdraw. First, that moving to the radical made the
  per-component test exact — the author's question about the Thomas
  decomposition not using saturated ideals exposed that the inequations need
  saturating too. Second, that the change from union to intersection was global
  — the author pushed back ("those ideals are ideals in the space of constants —
  we still want union"), and the resolution is union-across-cells with
  intersection-within, settled only once the cells were understood to overlap
  rather than partition. A duplicated `Input:` line was introduced and removed
  during an edit of the algorithm. On the other side, one confusion of the
  author's was clarified: the ∃-projection *gains* points rather than imposing
  an additional condition, since the open-set condition is the ∀ step over
  independent variables, which both variants perform.
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

**A new worked example, and the generalization to systems** (2–3 August 2026;
model Claude Opus 5).

- **Wrote the Navier–Stokes example section** (`06a67dc`). Figure 8 of the paper
  had drawn a nonlinear ansatz whose four top-level polynomials were declared to
  be the dependent variables of a system aimed at Navier–Stokes, and the paper
  promised the resulting system twice without ever writing it down. It was
  written down, reduced in closed form by hand, and the closed form turned out
  to be informative: **the ansatz reaches exactly the parallel shear flows and
  nothing else**, for two structurally distinct reasons — incompressibility
  annihilates the convective nonlinearity, the only nonlinear jet monomial's
  coefficient carrying the same factor as the continuity residual; and the
  residual is degree 1 in the leader while the ansatz's ODE is degree 2, so
  pseudo-division returns it unchanged and the projection forces μ·a·σ = 0
  outright, leaving no viscous solution with spatial structure. Written up as a
  new section presenting a **negative** result, which is what makes it worth
  including. The algebra was verified by an independently written sympy script
  committed alongside; the paper's Figure 8 notation collision (v used both for
  a velocity component and for the ODE's independent variable) was caught and
  renamed in the same commit.
- **Generalized the whole development from a single PDE to a system**
  (`93050df`, plus `0e5e5d7`, `0d77891`, `a10e350`, `3cc000f`). The
  restructuring is the author's, and his assessment of it was right and is worth
  recording as the reason it was safe: **theoretically this is a conjunction over
  j and nothing more**, because membership in an ideal is a condition on one
  differential polynomial at a time. The work was therefore mostly notational,
  but it touched abstract, introduction, algorithm, completeness theorem and
  corollary, since "the PDE" was woven through the prose. Four places needed
  more than substitution, and they are the substance: the *existence* loci do
  not decompose the way the membership locus does; power products must be
  collected within each remainder and never across them; the proof's conjunction
  over j works only because E:q^∞ is independent of j; and the coarse-splitting
  variant acquires a second defect that appears only once t > 1. Algorithm 3 was
  compacted so that every line is a one-line formula, which required naming the
  coefficient set — `Coeffs`, with a second argument naming the coefficient
  ring, renamed from `Coeff` and given its precedents in the literature after
  a short survey (`~/project/reports/newmethod-coefficient-set-notation-survey.md`).
- **Diagrammed the normalized ansatz and repaired the figure key** (`9500166`,
  `2583a62`). Introduced a diagram convention for normalized ansätze — a slot
  pinned by a normalization is written out literally, and the degree bound is
  written only when every generator of the ring may appear to that degree — and
  documented it in the key rather than leaving the new figure silently
  inconsistent with Figure 1. Also `a27858c`, a two-line alignment fix to
  Proposition 1.

The bullets that follow record a second connected arc, 4–5 August 2026 (model
Claude Opus 5), carried out at the author's direction: admitting inequations
into the target system, and following them through the loci, the completeness
theorem and the algorithm.

- **Traced the "involutive ⟹ coherent" citation to its end, and proved a usable
  substitute** (`44d26e9`, `747002c`). The author left a placeholder, "cite
  Li–Wang 1999, Theorem 3", against the claim that an involutive system is
  coherent. The companion note had already established that Li–Wang Thm. 3 is
  about *Wu*-passivity; checking the alternative showed that Bächler et al. §1
  asserts the implication in its introduction with no proof and no forward
  reference (three occurrences of "coheren" in the entire paper, two prose and
  one a bibliography entry), and that Robertz's book appendix, both CIRM
  versions, and Gerdt (2018) say nothing about coherence in their bodies at all.
  Gerdt's own DARCA 2008 slides list the involutive/coherent relation as a
  relation *to be found*, pointing at an unpublished "Minzlaff'06" that does not
  resolve to anything citable. Rather than cite an assertion, supplied a proof of
  the weaker statement a Thomas-based argument actually needs — a differentially
  simple system is coherent — from **Prop. 2.2.50** (a Δ-polynomial lies in E,
  hence pseudo-reduces to zero) and **BLOP Lemma 20** (a critical pair whose
  Δ-polynomial reduces to zero is solved), with (A1) supplying Lemma 20's
  distinct-leaders hypothesis; checked for circularity, Robertz's route to
  2.2.50 not passing through Rosenfeld's Lemma. Recorded at length as a new §5 of
  the companion note, which also flags a terminological trap the note's own
  implication lattice had hidden: Li–Wang's d-sim is *defined* as a coherent
  d-tri system, so "simple ⟹ coherent" is a tautology there, whereas Bächler
  et al.'s Def. 3.3 does not assume coherence and needs the theorem. The author
  kept a one-line remark in the paper and cut the proof (`155e716`).
- **Made the decomposition identity a numbered proposition, and corrected it**
  (`d0f6263`), at the author's request. Two defects surfaced in making the
  statement self-contained: the q on the left of √(E:q^∞) = ⋂(E^(i):(q^(i))^∞)
  is the product of the **inequations** of the undecomposed system, not of the
  initials and separants, and the running-text version had inherited E and q
  from the preceding paragraph — where q *is* the initials-and-separants product
  — so the displayed identity saturated by the wrong thing; and the components
  were indexed S_1…S_s while s already denoted the number of equations, leaving
  E^(s) ambiguous (now S_1…S_r, as in Robertz).
- **Put inequations into the target system and removed the homogeneous existence
  locus** (`aa3015f`). The restructuring is the author's proposal; supplied the
  analysis and the text. 𝒫 becomes a system (𝒫^=, 𝒫^≠) with Q the product of the
  inequations; V_{∃∖Ψ} disappears, being the ordinary existence locus with Ψ
  among the inequations, and three loci become two. Established that homogeneity
  was never part of the old definition but a statement of when the locus is worth
  computing — the formula makes sense without it — so it demotes to motivation;
  and that this gate is what had been blocking Navier–Stokes, whose convective
  term is quadratic in the unknown and so fails the paper's own homogeneity test,
  leaving it with no homogeneous existence locus at all. Inequations enter the
  membership side through a separate non-degeneracy locus N_Q, chosen over
  folding them into V_∀ so that the completeness theorem and its corollary need
  not be touched; the containment V_∀ ∩ N_Q ⊆ V_∃ is then unconditional, N_Q
  implying the ansatz-consistency hypothesis the old proof carried separately.
  Two points recorded because they are easy to get wrong: Q must be the *product*
  rather than the inequations taken separately, since √[A(c\*)] is radical but in
  general not prime, so the product demands a single solution satisfying all of
  them at once; and an inequation on a constant works through the same formula,
  via I:0^∞ = R.
- **Reformulated the membership locus on the saturated radical ideal**
  (`af0670c`). Proposed by the author; supplied the analysis, the enabling lemma
  and the downstream changes. V_∀ becomes {c\* : 𝒫^= ⊆ √[A(c\*)]:Q^∞}, which by a
  closure argument says exactly that each equation vanishes at every ansatz
  solution f with Q(f) ≠ 0 — the *guard* reading of the inequations. This repairs
  a real defect: the unsaturated condition rejects a value of the constants
  whenever **any** ansatz solution violates the system, including one on a branch
  the inequations were written to exclude. The saturation is free, by a two-line
  lemma (**Lemma 3**): for radical I, P ∈ I:Q^∞ iff QP ∈ I, since
  (QP)^k = (Q^k P)P^{k−1} and I is radical. A single power of Q suffices, so no
  saturation exponent is searched for, every reduction becomes a reduction of
  Q·P_j, and the Thomas decomposition is still computed from the ansatz alone —
  the paper's cost claim survives intact. Theorem 7, its proof, Corollary 8 and
  Algorithm 5 follow mechanically. Two cautions written into the text: the
  saturation is vacuous unless some irreducible component of the ansatz's
  solution set lies inside V(Q) — for a homogeneous linear ansatz with Q = Ψ the
  zero function is a limit of non-zero solutions rather than a component, so
  nothing moves, and Example 4 is now worked to show exactly that — and V_∀ is
  satisfied **vacuously** off N_Q, so the locus to report is V_∀ ∩ N_Q and
  omitting N_Q now over-reports where before it under-reported.
- **Errors of mine during this arc.** Three, each caught by the author's
  follow-up questions rather than by me. I claimed V_∀ was Zariski-closed and
  that admitting inequations into it would push set differences through the
  downstream apparatus; the algorithm's own Output specification already promises
  a constructible set presented as (𝔭, h) pairs, and the assembly already
  computes complements, so the objection was empty — found on re-reading the
  specification that the author's request for elaboration sent me back to. I
  framed "Q ≠ 0" as admitting exactly two readings, a universal assertion and an
  existential one, and recommended the existential; the author's saturated-radical
  proposal supplied a third, the guard reading, better than both and now the one
  the paper uses. And I recommended folding the non-degeneracy condition into
  V_∀'s definition before reading the completeness theorem, which is stated
  per-component over 𝒫^= and would have carried the change into its proof;
  retracted on reading it. A fourth, smaller: my first proof of Lemma 3 went
  through a decomposition into primes, where radicality and two lines suffice.
- **Flagged, not applied.** §5.2 describes the Thomas decomposition as
  partitioning the constant space into "quasi-projective varieties"; the cells
  are locally closed and in general reducible, so the term is wrong — they are
  constructible, or locally closed, sets. Raised twice and left to the author.
  On the author's question of where Chevalley's theorem on constructible sets is
  to be found: the original is Cartan–Chevalley, *Géométrie Algébrique*,
  Séminaire Cartan–Chevalley 1955/56, exposé 7 (Hartshorne's own attribution,
  verified in his bibliography), with **Hartshorne Ex. II.3.19** the citable
  statement and proof, **Matsumura Ch. 2 §6** the commutative-algebra treatment,
  and **EGA IV₁ Thm. 1.8.4** the general form.
- **A reference the author decided not to use.** For the parenthetical that in
  characteristic zero the radical of a differential ideal is again a differential
  ideal, located **Kaplansky, *An Introduction to Differential Algebra*
  (Hermann, 1957), Ch. I §4, Lemma 1.8**, with his Lemma 1.7
  (a^n ∈ I ⟹ (a′)^{2n−1} ∈ I, the division by n being where characteristic zero
  enters) and his characteristic-2 counterexample showing the hypothesis is
  needed. Kolchin Ch. I §2 Lemma 2 supplies the same machinery in his generality
  but leaves the radical corollary to an exercise. Drafted and then reverted at
  the author's request.

- **Renamed the existence locus to the consistency locus** (`d3dd938`, 5 August
  2026). The rename is the author's proposal, on the grounds that the paper has
  no analytic existence proof and the algebraic term is that the ideal is
  consistent; supplied the supporting check and the downstream changes. The
  check strengthened his case: **the paper never mentions a Nullstellensatz
  anywhere** — zero occurrences of "Nullstellensatz" or "Raudenbush" — so
  "existence locus" was not a name that slightly overreached a proved result but
  one for which the paper offered no warrant at all, while "consistent" was
  already a defined term in it. Three things travelled with the rename. The
  symbol V_∃ was **kept**, the subscript being load-bearing for the quantifier
  table and for V_{∃∀}, with a sentence noting that it records the quantifier
  pattern rather than a proved theorem. The quantifier table, not the name, is
  where the existence claim actually lived — its last column reads "vanishes at
  some one solution", whose use as a *characterization* needs consistency ⟹ a
  solution exists, the direction requiring a Nullstellensatz; the column is now
  marked as a reading rather than a definition, at no cost, since the paper's
  proofs only ever use the free direction. And the defining sentence was made
  explicit that the condition is on the ansatz *together with* the system, since
  "consistent" is separately used in the paper for the ansatz alone. Also
  recommended **against** the alternative of keeping the name and buying it with
  a citation — Robertz states the Ritt–Raudenbush analytic Nullstellensatz
  beside a proposition the paper already cites — because even then "existence"
  would mean an analytic solution on some domain rather than in a prescribed
  function space, and the hydrogen section turns on exactly that gap. Flagged in
  passing, and confirmed deliberate by the author: 92 lines of the *Physical
  Interpretation* subsection sit inside a `\begin{comment}` block and are not
  typeset.
- **Drafted the colloquium slide deck** (`289af84`, 5 August 2026), at the
  author's request, from the paper: `NewMethod-talk.tex`, 33 slides in Beamer,
  for a 50–60 minute general audience, with the optional frames marked so a
  ~25-minute version can be cut. Two choices were made on the author's behalf
  and flagged as easily reversed: making the completeness certificate the spine
  of the talk, and stating the non-physicality of the J₀ solution outright
  rather than eliding it. The deck is not part of the submitted paper.

The bullets that follow record a third connected arc, 17–18 August 2026
(model Claude Sonnet 5), carried out while the author was working directly
in the completeness subsection: cleaning up a notation leftover, formalizing
a paragraph he supplied into a corollary, correcting a saturation bug in the
membership-locus algorithm that he found, and adding a simpler
consistency-locus algorithm ahead of it at his request.

- **Cleaned up the A/S^= notation left over from the previous session's
  theorem rewrite, and formalized a paragraph into a corollary**
  (`b981e50`). Found unprompted, while reviewing the section at the
  author's general request to look at where he'd last been working:
  Theorem 11 had been reworded from "$A$" to "$S^=$" but its proof still
  said "$A$" throughout ($h(c^*)$, $[A(c^*)]$, etc.), a leftover from
  before the rename. At the author's direction — "let's just use S for
  the equations of a differentially simple system, since we don't use the
  inequations for anything in the theorem" — replaced both the leftover
  $A$ and the theorem statement's $S^=$ with a single $S$, matching the
  paper's own existing convention for a generic system $S$ (§2.1); fixed
  three "constuctable"/"constructable" typos in the same pass. Then, given
  a paragraph the author had written ("Theorem 10 shows how a variety...
  give a constructible set") and asked to have "worked up into a
  Corollary," formalized it as the new **Corollary 12**: equations
  $\mathcal{E}$ and inequations $\mathcal{N}$, each reduced against the
  same component, assembled by intersection-minus-union into one
  constructible set, with a proof from the definition of constructible
  plus the per-polynomial theorem. Corollary 13 (the assembly step) was
  rewritten to cite Corollary 12 explicitly — taking
  $\mathcal{E}=\mathcal{P}^=$, $\mathcal{N}=\emptyset$ — instead of
  invoking the single-polynomial theorem directly, closing a gap where it
  had asserted a per-component constructible set without ever justifying
  where it came from.
- **Corrected a saturation bug in Algorithm MembershipLocus that the
  author found, getting it wrong once before he supplied the fix**
  (`f99f416`). The author read line 11 and objected: $W_i$ reused $h_i$ —
  the cell's own constant-only inequation product — as its saturator,
  when $W_i$ needs to know whether a component degenerates the *target's*
  declared inequations, a different fact entirely. **My first attempt was
  wrong, twice over, and the author corrected both parts in one message**:
  I had proposed reducing the target's combined inequation $Q$ directly
  inside the per-component loop and multiplying its coefficients into
  $h_i$; the author rejected this — "I don't think we want Q on line 11 -
  we want to reduce the system's $S_i^{\ne}$, and on line 12 I don't think
  we need $h_i$" — without yet saying how $Q$ should enter at all. Once
  redirected, reducing the (now correctly scoped) $S_i^{\ne}$ as a single
  product resolved the second part on its own: an already-constant
  inequation reduces to itself under `FullReduce` and passes through the
  product unchanged, so $h_i$'s contribution is subsumed automatically and
  no separate multiplication is needed, matching the author's second
  objection exactly. That left $Q$'s own role unaccounted for, which I
  flagged as an open question; the author supplied the actual resolution
  in the next message — feed $Q$ into the `DifferentialDecompose` call on
  line 1 as the system's inequation, $(A(c),\{Q\})$, so the Thomas
  decomposition's own case-splitting absorbs $Q$ into every component's
  $S_i^{\ne}$ automatically, which the already-fixed reduction step then
  picks up with no further change. Updated the surrounding prose and the
  $W_i=\mathcal{V}_i\cap C_i$ identity (now $W_i=\mathcal{V}_i\cap
  C_i\cap N_i$) to match, and declared
  $\mathcal{P}^{\ne}=\{Q_1,\ldots,Q_w\}$ in the algorithm's Input, which
  the pre-existing Output line had referenced without ever declaring.
- **Added Algorithm ConsistencyLocus** (`62a9f8e`), at the author's
  request and to his specification: "it just runs DifferentialDecompose on
  the whole system and then we have to pull out the constants." Supplied
  the formal write-up — decompose $(A(c)\cup\mathcal{P}^=,\{Q\})$ in one
  call, then read $E_i,h_i$ off each component with no per-equation
  reduction loop, no $g_i$ refinement and no SMPD assembly — and the
  justification for why the simpler version is sound: the target's
  equations and inequations are decomposed as *part of* $S_i^=$ and
  $S_i^{\ne}$ rather than tested afterward, so cell membership alone, via
  the same \cite{ThomasDecomp} Remark 2.3 extension argument the
  completeness section already uses, guarantees a genuine solution —
  unlike MembershipLocus, which needs the finer per-equation
  ideal-membership test because it decomposes the ansatz alone. Inserted
  as a new §2.4 ahead of MembershipLocus (which shifts to §2.5 without any
  label changes), both for completeness and, per the author's own stated
  reasoning, "to ease the reader in with a simpler algorithm" before the
  more elaborate one. Softened the §2.3 closing paragraph, which
  previously described the direct approach only to declare it infeasible
  and drop it, and added a one-sentence transition into MembershipLocus
  tying the two algorithms together.
- **Typed ideals as ideals throughout both algorithms, after a research
  dispatch and one self-caught wrong fix** (`4491431`). The author pushed
  on `K_i`/`H_i`'s precision across several exchanges — first asking for
  the exact `(T,h)`-regular-system representation the CTD paper's
  Difference algorithm builds, then, on being told `H_i` (a product of
  several possibly-non-principal primes) couldn't honestly be forced into
  a single polynomial without either an approximation or a combinatorial
  blow-up, saying plainly "I don't want an approximation... let's just
  write it that way in the paper" ($V(\mathfrak p)\setminus V(\mathfrak
  q)$, literally) — and then noticed the deeper issue himself: the
  pseudocode's own convention (bare symbols are finite generator lists,
  per `minAss`'s stated signature) made `𝔮⊆𝔭` and `∏𝔮` ill-typed as
  written, and asked whether treating ideals as a primitive type is
  standard practice in a professional journal article. Rather than settle
  that by assertion, dispatched a **Claude Opus** subagent with the two
  algorithms' full text and instructions to read the actual cited sources
  (the CTD paper, its 2008 ConstructibleSetTools companion, and — found
  independently by the subagent — Kurata–Nabeshima 2024, a directly
  analogous comprehensive-primary-decomposition algorithm that is fully
  ideal-typed with no bracketing at all) and come back with a concrete,
  worked recommendation rather than a principle. Verdict: type the ideals
  as ideals. Applied in full: `minAss` restated as `ideal → set of
  ideals` (`minAss((1))=∅`), which alone makes `K_i`, `H_i` and the
  containment test correct exactly as already written; ideals of
  $\mathbb Q[c]$ set in fraktur; a new $\mathrm Z(\cdot)$ operator (citing
  CTD's own §2 notation) separating a pair's representation from the
  point set it denotes, closing several places (`C_i`, `W_i`, `D`) where
  the same symbol had been doing both jobs at once; `H_i` renamed to the
  ideal $\mathfrak h_i$, letting the ad hoc "product of finite sets"
  definition be deleted outright; the now-redundant `𝔭≠(1)` test dropped.
  The subagent also caught two defects unprompted: a mathematical slip in
  my own prose from the `K_i`/`H_i` pass ("$\langle q\rangle$ is prime" —
  false whenever $q$ is reducible; corrected to go through $q$'s
  irreducible factors), and that `SMPD`'s citation was wrong for what it
  was actually being fed — CTD's own Algorithm 3 works only on *regular
  systems* (a chain plus a single inequation), not the general
  `(ideal, ideal)` pairs this algorithm produces. **My first pass at that
  second fix was itself wrong**: I re-cited the operation to the 2008
  ConstructibleSetTools paper (which does handle general constructible
  sets) but left the name `SMPD` — CTD's own name for its own narrower
  algorithm — unchanged, so the name and citation no longer agreed. The
  author caught it directly ("didn't the report say that SMPD doesn't
  work on pairs of ideals and suggests using some other algorithms from
  the Chen paper?"). Corrected by renaming to `RefiningPartition`
  (ConstructibleSetTools's actual name) and using its real return shape —
  pieces tagged with which family members they came from — to simplify
  the assembly step's containment test into a direct provenance lookup,
  which is what the subagent's report had recommended and my first pass
  had dropped.

**Computational infrastructure (no paper text).** The following is tooling and
measurement on the directory's computation scripts, recorded here for
completeness rather than because it produced prose; 5 August 2026 (model
Claude Opus 5) and 17–18 August 2026 (model Claude Sonnet 5).

- **Added `--gtz-subprocess` to `thomas-ansatz-solve.sage`** (`7f3e1a6`). The
  occasion was a lost run: the `navier-stokes` / ansatz-25 `--generic-cell`
  computation ran 32 hours inside `minimal_associated_primes` on c200-1 and was
  killed by a mains outage, leaving no trace past the line "entering GTZ". The
  underlying finding is a negative one, and worth recording because it is not
  obvious: the in-process call **cannot be made to report progress**.
  primdec.lib's `minAssGTZ` → `minAssGTZ_i` chain carries no `dbprint`
  instrumentation, so `printlevel` adds nothing, and `option(prot)` does not
  reach stdout through Sage's libsingular wrapper — verified both ways,
  `opt['prot'] = True` and `opt_ctx(prot=True)`, against `minAssGTZ` and
  against a plain cyclic-6 `std`, with no output either way, while the
  identical computation in a standalone Singular prints the full protocol
  stream. The new mode therefore runs the same `minAssGTZ` in a standalone
  Singular under `option(prot)` and `option(mem)`, writing the primes to their
  own file so the protocol stream never has to be untangled from the result. It
  buys a tailable log, a child process that `--gtz-timeout` can bound and kill,
  and an on-disk cache keyed by the ideal's content hash so an interrupted run
  resumes rather than restarts. Checked against the in-process path on
  hydrogen/5 `--generic-cell`: identical five primes and two GENUINE varieties,
  a cache hit on re-run, and the timeout firing with no orphaned children. Two
  defects were found and fixed while testing — a generated identifier `res`
  collides with one primdec.lib already owns, and **Singular exits 0 even after
  a hard error**, so the exit status alone is not a success test and the log is
  now scanned for error lines.
- **Two diagnostic findings from the same run.** The stray line `ZERO` that had
  been noted as unexplained output from `minAssGTZ` is a leftover debug print in
  Singular's own `primdec.lib:9401`, inside `primdec`'s independent-set fast
  path (`if (dim(j)!=d) {"ZERO"; return(0);}`) — it records only that the fast
  path declined, and is not a result. Separately, the coherence checks for
  `navier-stokes` ansatzes 25, 25.2 and 25.3, previously recorded as never
  launched, in fact completed on 4 August: all three are **passive** (264, 192
  and 180 integrability conditions respectively, every one reducing to zero).
  With 25.34, that is four generic cells measured coherent, so for these the
  pseudo-remainder against the hand-built generic cell is a genuine normal form
  and the corresponding loci are not merely lower bounds. This is a statement
  about four measured cells, not a theorem about generic cells in general.
- **Conformed `thomas-ansatz-solve.sage` to the paper's updated Algorithm
  MembershipLocus** (`51d748f`, 17–18 August 2026), at the author's
  direction ("the script should conform to the paper's algorithm").
  Dropped the script's older "Lemma 3" mechanism — reduce $Q\cdot P_j$
  for each target equation against the ansatz-alone decomposition,
  relying on the cell's ideal being radical to fold the saturation in —
  since the paper no longer works that way: $Q$ is now fed into
  `DifferentialDecompose` itself as the system's inequation, so every
  component's own inequations already carry $Q$'s non-degeneracy, and
  each $P_j$ is reduced plain. Removed the script's separate `Jq`/
  `cell_off_nq`/"off N\_Q" bookkeeping and `--split-inequations` flag,
  since Q-derived inequations are now ordinary members of a cell's
  $S^{\ne}$ and the script's pre-existing `ineq_coeffs` per-inequation
  pruning already catches them for free — net $-69$ lines. Net change
  84 insertions / 153 deletions. Validated on hydrogen/5: `16` cells
  (down from the stale `29`-cell baseline the header comment had
  asserted, expected once $Q$ enters the decomposition), full run
  26m38s, 6 distinct genuine solution varieties recovered including
  recognizable hydrogen eigenvalues ($E=-1/8$, $E=-1/2$), no errors.
  Two operational mishaps along the way, neither affecting the result:
  a background fork stalled with zero tool calls on its first attempt
  and had to be explicitly resumed to actually do the work, and the
  resumed fork accidentally launched the same 26-minute validation run
  twice in parallel (caught and the duplicate process killed). The
  commit itself landed under the wrong message — an unrelated
  concurrent task-runner auto-commit sweep picked up the uncommitted
  script change and stamped it with an unrelated task's commit message
  ("task-runner: conf-audio-autojoin-start-limit") — caught by the
  author, and amended to describe what actually changed.

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

**The computation scripts.** This directory also holds the Sage and Python
scripts that run the method — `joca.sage` and its variants, `joca-thomas.sage`,
`joca-rg.sage`, `joca-rg-combined.sage`, `thomas-ansatz-solve.sage`,
`ansatz-library.sage`, `rg_basefield.py`, the coherence and comparison probes,
and the directory's `README.md` documenting them. **Treat all of it as
AI-written**, from the author's specifications and mathematics and reviewed by
him. Some sixty further commits, not itemized here, do nothing but develop these
scripts. They are given a blanket attribution rather than a per-commit one for
two reasons: essentially the whole of `thomas-ansatz-solve.sage` and
`ansatz-library.sage` was written by the AI, so there is nothing to
disambiguate; and the scripts are separate artifacts from the paper, to which
the journal's disclosure requirement attaches. Where a computation changed what
the paper *says*, that is recorded above as a paper contribution, not here.

**Artifacts.** The **Graphical Abstract** is AI-generated (`6db3f5b`). The
colloquium slide deck `NewMethod-talk.tex` (`289af84`), the twelve companion
notes listed above, and this record itself are AI-produced.

**Note for the formal declaration — an inventory, not a decision.** Elsevier's
generative-AI policy for authors is framed around assistance with *language and
readability*, and its declaration template is worded for that scope; AI use in
the *research* process and AI-generated *images* fall under separate provisions.
Drafting the submitted statement is a later step. What follows is the inventory
that drafting will need, so the classification does not have to be redone from
scratch; the journal's current requirements must be confirmed at that time,
since these policies have been revised repeatedly. Items recorded above that
fall outside "language and readability":

- **Mathematical drafting.** The completeness theorem and its proof, drafted by
  the AI and inserted by the author (`843f9ea`, inlined and edited `874a9c7`),
  and the regularity counterexamples that entered with it and were later moved
  out of the paper into a companion note (`933ef24`).
- **Structural and expository work on the author's ideas.** The generalization
  from a single PDE to a system (`93050df` and follow-ons), the formalization of
  the algorithm together with the 30–31 July citation audit that corrected it
  (`8de5053` and siblings), and the rewritten Organization paragraph
  (`a59b8ae`).
- **A section written outright.** The incompressible Navier–Stokes example
  (`06a67dc`) — a negative result: the ansatz reaches exactly the parallel shear
  flows.
- **Literature search and prior-art assessment,** including the AI's negative
  verdict on the novelty of a pipeline it had itself proposed (`5d56eb6`,
  `bc7a8c4`), and citation suggestions from a second tool, OpenAI's GPT-5
  (`55487de`).
- **Results the paper states that came from AI-written computation.** The
  scripts are separate artifacts (see above), but where a run changed what the
  paper *says*, it is the research-process provision rather than the writing
  provision that is the relevant one.
- **An AI-generated image.** The Graphical Abstract (`6db3f5b`). Of the three
  tracks this is the most restrictive, and it may bear on whether the artifact
  can be submitted at all, not merely on how it is described.

Items outside the submission itself — the twelve companion notes, the
colloquium deck (`289af84`), and this record — are noted here for completeness;
whether they require disclosure at all is part of the deferred decision.

---

## Appendix — commit index

Every commit to the paper or to a companion note that the AI wrote or that
landed AI-drafted material, in date order. The narrative above is selective;
this is not. **Scope and its limits:** the table covers `NewMethod.tex`,
`NewMethod-talk.tex` and the companion `.tex`/`.md` notes. It deliberately
excludes three classes — the computation scripts and their `README.md` (blanket
attribution above; roughly sixty further commits), commits that only rebuild the
PDF or only maintain this record (about twenty), and the author's own prose
commits, of which there are some two hundred and which are his work. So the
table is complete for what it covers and says nothing about what it does not.

*Whose idea* distinguishes **author** (his proposal, question or objection; the
AI supplied the execution), **AI** (the AI's own finding or initiative,
reviewed by the author), and **joint** (developed in discussion, neither party's
alone). It is taken from the session transcripts, not from the commit's author
field. The seven rows in *italics* are the author's own commits that carried AI
material, listed here because the author field alone would hide them.

To update this record: run `git log --author='Claude Code'` over `NewMethod/`,
diff it against the Commit column, and only the new rows need thought.

| Commit | Date | Artifact | What changed | Whose idea |
|---|---|---|---|---|
| *`226a2a5`* | 2025-12-31 | NewMethod.tex | AI's fixes for the overfull-hbox warnings | AI |
| *`55487de`* | 2026-02-27 | NewMethod.tex | citation suggestions from **GPT-5** | GPT-5 |
| *`6db3f5b`* | 2026-04-22 | GraphicalAbstract | AI-generated graphical abstract | AI |
| *`843f9ea`* | 2026-05-28 | NewMethod.tex | inserted the AI-drafted completeness theorem and proof | joint |
| `02b57b8` | 2026-05-28 | NewMethod.tex | bibitem for the squarefree-regular-chain definition | AI |
| `874a9c7` | 2026-05-29 | NewMethod.tex | inlined the theorem, proof and regularity counterexamples | joint |
| `57ce904` | 2026-05-29 | NewMethod.tex | specialization result promoted to Lemma 1; amsthm numbering | AI |
| `67aeb0f` | 2026-05-29 | NewMethod.tex | hypothesis (3) stated generically on the parametric ansatz | AI |
| `cc6ee77` | 2026-05-29 | NewMethod.tex | all four hypotheses verified in the ODE-ansatz example | AI |
| `7c18e49` | 2026-05-29 | NewMethod.tex | hypothesis check moved into §3 as a capstone | AI |
| `52664d7` | 2026-05-29 | NewMethod.tex | stale motivation of the [A]∩R₀=(A) conjecture fixed | AI |
| `6e1ee75` | 2026-05-29 | NewMethod.tex | stopped re-defining leader / reduction in the completeness §| AI |
| `aba1b6e` | 2026-05-29 | NewMethod.tex | Projection subsection tightened, reconciled with Theorem 1 | AI |
| `0aee721` | 2026-05-29 | NewMethod.tex | system of equations displayed rather than run in-line | AI |
| `2b71ff5` | 2026-05-29 | NewMethod.tex | bibliography reformatted to Elsevier numbered style | author |
| `e2a7449` | 2026-06-04 | NewMethod.tex | restructured into core + general algorithm | author |
| `c05c04a` | 2026-06-04 | NewMethod.tex | §2 folded into §1; "Why not Rosenfeld–Gröbner alone?" added | author |
| `ff291ad` | 2026-06-04 | NewMethod.tex | comprehensive step recast as cited prior art (Fakouri; Thomas) | AI |
| `102fd6a` | 2026-06-04 | NewMethod.tex | positioned against the symmetry-reduction tradition | AI |
| `806b526` | 2026-06-04 | NewMethod.tex | **corrected** the separability claim (parabolic / Runge–Lenz) | AI |
| `3c0818c` | 2026-06-04 | NewMethod.tex | one shared counter for all theorem-like environments | AI |
| `315ab7c` | 2026-06-05 | NewMethod.tex | **corrected** "RG returns the ansatz unchanged" | AI |
| `cc09d29` | 2026-06-05 | NewMethod.tex | Clarkson–Mansfield 1994 cited beside the direct method | AI |
| `59ccab9` | 2026-06-05 | rg-saturation note | **new note**: the dropped-saturation bug and the bad locus as a pole | AI |
| `b659ce8` | 2026-06-05 | rg-saturation note | layout fix (P2/P3 bled off the page) | AI |
| `243b96d` | 2026-06-08 | ck-direct-method-examples | **new note**: worked CK examples (heat kernel, hydrogen J₀) | AI |
| `3d39b0a` | 2026-06-11 | method-comparison.md | **new note**: 525-line comparison against CTD/CGS/GTZ, P-RG, Thomas | AI |
| `18ddc11` | 2026-06-12 | NewMethod.tex | **corrected** "the bad locus is empty" — B is {a₀=a₁=0} | AI |
| `d965dbb` | 2026-06-12 | NewMethod.tex | **corrected** the general hypothesis-(4) linear-leader claim | AI |
| `65fe337` | 2026-06-12 | method-comparison.md | worked Δ(Ψₓ,Ψ_y) split as an atlas-vs-image example | AI |
| `901ce6a` | 2026-06-18 | partial-strata note | **new note**: names the partial-solution strata | joint |
| `5d56eb6` | 2026-07-02 | bounded-prolongation | algorithm exposition + prior-art survey (negative novelty verdict) | joint |
| `bc7a8c4` | 2026-07-05 | bounded-prolongation | **corrected**: elimination must be (J:Ψ^∞)∩K[c], not J∩K[c] | joint |
| `3bd3828` | 2026-07-07 | gcrd-closure note | **new note**: GCRD construction recovering partial strata | joint |
| `7d083b0` | 2026-07-07 | gcrd-closure note | wrote out "equation (2)" in algorithm step 1 | AI |
| `6f2d20a` | 2026-07-07 | riquier note | **new note**: Cauchy–Kovalevskaya, Riquier, ansatz-5 data | AI |
| `c9f682a` | 2026-07-07 | riquier note | converted from markdown to LaTeX | author |
| `72fe9f9` | 2026-07-09 | partial-strata note | PDE counterexample: autoreducedness does not preclude strata | AI |
| `b8e2cec` | 2026-07-10 | prolongation-projection | **new note**: the author's algorithm developed; quantifier discipline | author |
| `95d687c` | 2026-07-13 | prolongation-projection | the staged route's loss located; fixed-chain prolongation | AI |
| `97f718c` | 2026-07-18 | ansatz-method-provenance | **new note**: four-layer provenance of the "ansatz method" term | AI |
| `4063c0c` | 2026-07-20 | hierarchy note | **new note**: coherent/regular/simple/passive pinned to sources | AI |
| `fd8f64d` | 2026-07-20 | hierarchy note | duality caveat; parametric RG as a union of intersections | AI |
| *`933ef24`* | 2026-07-23 | regularity-counterexamples | AI subsections cut from the paper, preserved in a new note | author |
| `f3fda04` | 2026-07-23 | equivalence note | **new note**: membership ↔ homogeneous existence | AI |
| `50cf452` | 2026-07-23 | equivalence note | "hypotheses and their necessity" remark | AI |
| `f3ba3f5` | 2026-07-23 | partial-strata note | order-matching counterexample; ansatz-5 has no genuine strata | AI |
| *`53c4eba`* | 2026-07-27 | NewMethod.tex | author's rewrite incorporating the AI's repaired containment proofs | joint |
| `7532ee5` | 2026-07-28 | NewMethod.tex | preliminaries rewritten onto simple systems | author |
| `68af580` | 2026-07-28 | NewMethod.tex | Rosenfeld's Lemma stated inline, fixing a bug the AI had just introduced | AI |
| `41a7102` | 2026-07-28 | hierarchy note | **corrected**: Wu-passive separated from Janet-involutive | AI |
| `d015a61` | 2026-07-28 | NewMethod.tex | completeness rebuilt on Robertz Prop. 2.2.50 (−54 lines) | AI |
| `3d92072` | 2026-07-29 | NewMethod.tex | denominator hypothesis dropped; algebraic simplicity as (A1)–(A3) | AI |
| `93a031b` | 2026-07-29 | NewMethod.tex | **corrected** the Robertz citation and proposition numbers | AI |
| `6617bdd` | 2026-07-29 | NewMethod.tex | Robertz proposition numbered into the paper's sequence | author |
| `d9c4c8d` | 2026-07-30 | NewMethod.tex | membership locus on the radical; theorem becomes an equivalence | author |
| `83a4465` | 2026-07-30 | NewMethod.tex | how a component's cell is obtained | AI |
| `26ff25e` | 2026-07-30 | NewMethod.tex | assembly by intersection, not union (after author's pushback) | joint |
| `1380cf6` | 2026-07-30 | NewMethod.tex | second, cheaper assembly branch; the intermediate locus | author |
| `232cec8` | 2026-07-30 | NewMethod.tex | prime decomposition hoisted into the algorithm's loop | author |
| `07ab265` | 2026-07-30 | NewMethod.tex | a computational recipe for step 7a | AI |
| `ab1166c` | 2026-07-30 | NewMethod.tex | hand-rolled recipe **replaced** by the published disjoint refinement | AI |
| `4afff4d` | 2026-07-30 | NewMethod.tex | Algorithm 3 reformatted to the cited papers' house style | AI |
| `b585d77` | 2026-07-30 | NewMethod.tex | ansatz consistency becomes a computed by-product | author |
| `f4996bb` | 2026-07-30 | NewMethod.tex | per-branch refinement; **corrected** the AI's own Output spec | AI |
| `3670235` | 2026-07-30 | NewMethod.tex | algorithm carries finite generating sets, not ideals | author |
| `8de5053` | 2026-07-30 | NewMethod.tex | **corrected**: SMPD, replacing a subroutine name the AI invented | author |
| `bf45291` | 2026-07-30 | NewMethod.tex | table of the algorithm's four subroutines | AI |
| `3efd8bd` | 2026-07-30 | NewMethod.tex | ring notation migrated to Ritt/Kolchin, 17 sites | AI |
| `dd996e2` | 2026-07-30 | NewMethod.tex | Example 2 notation; numbered subroutine citations; lines 4–6 | AI |
| `2553068` | 2026-07-30 | NewMethod.tex | **corrected**: line 3 is Ritt full reduction | author |
| `631e6e9` | 2026-07-30 | NewMethod.tex | FullReduce cited to Ritt §6 | author |
| `260d919` | 2026-07-30 | NewMethod.tex | **corrected** the Reduce row — Alg. 3.3, and what it returns | author |
| `d35e9cd` | 2026-07-30 | NewMethod.tex | minAss cited to GTZ Cor. 3.2(v), §9 | AI |
| `47e9360` | 2026-07-30 | NewMethod.tex | r shown in the base ring of the first two ansatz diagrams | AI |
| `0e34162` | 2026-07-30 | NewMethod.tex | subscripts for Δ-derivatives; primes only inside ODE extensions | AI |
| `a59b8ae` | 2026-07-31 | NewMethod.tex | organization paragraph rewritten (it promised a missing algorithm) | AI |
| `0b3c482` | 2026-07-31 | NewMethod.tex | forward reference from §1 dropped | AI |
| `e953bc8` | 2026-07-31 | NewMethod.tex | constants-as-a-block rationale restated beside the cells | AI |
| `125b8ae` | 2026-07-31 | NewMethod.tex | bibliography ordered by first citation, per Elsevier | author |
| `06a67dc` | 2026-08-02 | NewMethod.tex | **new section**: incompressible Navier–Stokes, a negative result | AI |
| `93050df` | 2026-08-03 | NewMethod.tex | whole development generalized from a PDE to a system | author |
| `0e5e5d7` | 2026-08-03 | NewMethod.tex | Algorithm 3 lines 8–9 as subring intersections | AI |
| `0d77891` | 2026-08-03 | NewMethod.tex | coefficient set named, collapsing line 6 | AI |
| `a10e350` | 2026-08-03 | NewMethod.tex | Coeff given a second argument naming the coefficient ring | author |
| `3cc000f` | 2026-08-03 | NewMethod.tex | renamed Coeffs; precedents cited after a notation survey | AI |
| `9500166` | 2026-08-03 | NewMethod.tex | normalized ansatz 25.34 diagrammed at the head of §4 | AI |
| `2583a62` | 2026-08-03 | NewMethod.tex | Figure 9 degree labelling fixed; diagram key extended | AI |
| `a27858c` | 2026-08-03 | NewMethod.tex | two text lines of Proposition 1 left-aligned | AI |
| `44d26e9` | 2026-08-04 | NewMethod.tex | coherence of simple systems proved, replacing an uncitable claim | author |
| `747002c` | 2026-08-04 | hierarchy note | new §5: BGLR-simple ⟹ coherent, with the terminological trap | AI |
| `d0f6263` | 2026-08-04 | NewMethod.tex | decomposition identity made a Proposition; two defects **corrected** | author |
| *`155e716`* | 2026-08-04 | NewMethod.tex | author kept a one-line remark and cut the AI's coherence proof | author |
| `aa3015f` | 2026-08-05 | NewMethod.tex | inequations into the target system; V_{∃∖Ψ} removed | author |
| `af0670c` | 2026-08-05 | NewMethod.tex | membership locus on the saturated radical ideal; Lemma 3 | author |
| `289af84` | 2026-08-05 | NewMethod-talk.tex | **new artifact**: 33-slide colloquium deck | author |
| `d3dd938` | 2026-08-05 | NewMethod.tex | existence locus renamed the consistency locus | author |
| `b981e50` | 2026-08-17 | NewMethod.tex | A/S^= notation **corrected**; typos fixed; Corollary 12 formalized | joint |
| `f99f416` | 2026-08-17 | NewMethod.tex | **corrected** Algorithm 7's saturation bug, after two wrong AI attempts | joint |
| `62a9f8e` | 2026-08-18 | NewMethod.tex | **new algorithm**: ConsistencyLocus, ahead of MembershipLocus | author |
| `4491431` | 2026-08-18 | NewMethod.tex | ideals typed as ideals throughout; **corrected** RefiningPartition naming | joint |

*/Claude Opus 4.8 (this record drafted by the AI it documents); updated
27 July 2026, 30 July 2026 and 5 August 2026 by /Claude Opus 5, and audited
against the full git history and the session transcripts on 5 August 2026 by
/Claude Opus 5; updated 18 August 2026 by /Claude Sonnet 5.*
