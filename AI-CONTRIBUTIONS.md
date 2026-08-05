# A.I. use — record of AI contributions to the *NewMethod* paper

*Working record, kept so the required journal statement on generative-AI use
can be written from fact rather than memory. Reconstructed from this
repository's git history; covers contributions through 5 August 2026 and should
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
27 July 2026, 30 July 2026 and 5 August 2026 by /Claude Opus 5.*
