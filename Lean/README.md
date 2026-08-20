# `NewMethod/Lean` — the two algorithms in Lean 4 / Mathlib

A Lean 4 / Mathlib formalization of the **two algorithms** of
`~/Papers/NewMethod/NewMethod.tex`:

| Algorithm | `NewMethod.tex` | Lean |
|---|---|---|
| 1. `ConsistencyLocus` (`alg:consistency`) | line 1594 | `NewMethod.AlgorithmInput.consistencyLocus` |
| 2. `MembershipLocus` (`alg:core`) | line 1657 | `NewMethod.AlgorithmInput.membershipLocus` |

together with every supporting notion their Input/Output clauses mention.

The goal is **fidelity of statement, not proof-completion**.  A reader should
be able to put `NewMethodAlgorithms.lean` beside the LaTeX and check, line by
line, that the Lean says what the paper says.  Every `axiom` names its source,
the source's section/theorem/definition number, and the local file the
statement was read from.  Every departure from a source's own phrasing carries
a `NOTE (departure from source):` line; every suspected error in the paper or a
cited source carries a `NOTE (suspected error in source):` line and is **not**
silently repaired — the Lean follows the paper.

## Layout

There is a **single self-contained file**, `NewMethodAlgorithms.lean`
(1487 lines), with `import Mathlib` at the top.  This is the constraint imposed
by the checker (see *Checking* below): AXLE compiles one file.  No modular
mirror is provided — a stale mirror would be worse than none — so the file is
organized with `/-! # §N … -/` module-doc headers instead:

| § | Contents |
|---|---|
| §1 | differential ideals, saturation, radical |
| §2 | the input data of both algorithms (the Input clauses) |
| §3 | systems, solutions, simple systems (ThomasDecomp §§2–3) |
| §4 | constructible sets, `V(·)`, `Z(·)` |
| §5 | cells and the projection fact |
| §6 | the four subroutines |
| §7 | the three constant loci |
| §8 | the two algorithms and their output correctness |
| §9 | bridge to `~/project/Completeness.lean` |
| §10 | Propositions `simple-membership` and `decomposition-membership` |

## Checking

```
~/axiommath.ai/check.py ~/Papers/NewMethod/Lean/NewMethodAlgorithms.lean
```

Last run (environment `lean-4.29.0`), verbatim:

```
okay: True
errors: 0, warnings: 8, infos: 0
--- lean warnings (8) ---
-:1083:8-1083:32: warning: declaration uses `sorry`
-:1272:8-1272:36: warning: declaration uses `sorry`
-:1297:8-1297:38: warning: declaration uses `sorry`
-:1310:8-1310:16: warning: declaration uses `sorry`
-:1373:4-1373:20: warning: declaration uses `sorry`
-:1377:4-1377:16: warning: declaration uses `sorry`
-:1381:4-1381:11: warning: declaration uses `sorry`
-:1425:6-1425:16: warning: declaration uses `sorry`
```

All eight warnings are the declared `sorry` debts tabulated below.  There are
no errors and no other warnings.

## Section-by-section map: `NewMethod.tex` → Lean

| `NewMethod.tex` | what it says | Lean declaration |
|---|---|---|
| 1082–1085 | differential ideal, `[g₁,…,g_n]` | `IsDiffIdeal`, `diffIdealBy` |
| 1087 | saturation `I : S^∞` | `satBy` (single element), `satByIdeal` (ideal) |
| 1089 | consistency `1 ∉ I` | `Consistent` |
| 1097–1108 | ranking, compatibility conditions | `AlgorithmInput.rank`, `.rank_deriv_mono`, `.rank_lt_deriv` |
| 1120–1126 | leader, initial, separant, `H_A = I_A ∪ S_A` | `.leader`, `.initial`, `.separant`, `initSep` |
| 1148–1152 | Ritt's full reduction contract | `FullReduce`, `fullReduce_spec` |
| 1211–1213 | a system `(S^=, S^≠)`; `Sol` | `DiffSystem`, `Sol` |
| 1215–1226 | Janet division, non-multiplicative prolongation | (folded into `Involutive`, stub) |
| 1229–1231 | involutivity | `Involutive` (stub) |
| 1236–1250 | algebraic simplicity (A1)–(A3) | `AlgebraicallySimple` |
| 1253–1263 | differential simplicity (S1)–(S4) | `DifferentiallySimple` |
| 1265–1266 | differential Thomas decomposition | `IsThomasDecomposition` |
| 1281–1297 | Prop. `simple-membership` | `simple_membership` |
| 1300–1324 | Prop. `decomposition-membership` | `decomposition_membership` |
| 1336–1345 | `𝒫 = (𝒫^=, 𝒫^≠)`, `Q = Q₁⋯Q_w` | `.Peq`, `.Pne`, `Q` |
| 1360–1362 | non-degeneracy locus `N_Q` | `Nq` |
| 1371–1373 | consistency locus `V_∃` | `Vexists` |
| 1387–1389 | membership locus `V_∀` | `Vforall` |
| 1404–1432 | Thm. `locus-containment`, `V_∀ ⊆ V_∃` | `Vforall_subset_Vexists` (**proved**) |
| 1537–1545 | CTD Def. 6, constructible set | `IsConstructible` |
| 1547–1552 | the cell of a differentially simple system | `cell` |
| 1554–1561 | the Remark 2.3 projection fact | `projection_eq_cell` |
| 1571–1592 | `⟨F⟩`, `Z(𝔞,𝔟)`, `Z(𝒲)`, `V(∏) = ⋃V`, `V((1)) = ∅` | `CPair`, `Zp`, `Zw`, `Vz_mul`, `Vz_listProd`, `Vz_one` (**proved**) |
| **1594–1626** | **Algorithm 1** | **`consistencyLocus`** |
| 1618–1619 | its Output clause `Z(𝒲) = V_∃` | `consistencyLocus_correct` (`sorry`) |
| 1628–1641 | `C_i = V(E_i) \ V(h_i)` | `cell_eq_Zp` (**proved**) |
| **1657–1712** | **Algorithm 2** | **`membershipLocus`** |
| 1678–1679 | the flag `μ ∈ {∀, ∃∀}` | `Assembly` |
| 1685–1689 (lines 3–7) | `J_i` | `Ji` |
| 1690–1691 (lines 8–9) | `E_i`, `h_i` | `Ei`, `hi` |
| 1692 (line 10) | `C_i` | `Ci` |
| 1693 (line 11) | `𝔥_i` | `hfrak` |
| 1694 (line 12) | `W_i` | `Wi` |
| 1704 (line 15) | the `RefiningPartition` family | `refiningFamily` |
| 1705 (line 16) | the `C_i ∈ ι ⟹ W_i ∈ ι` filter | `keepBlock` |
| 1706 (line 17) | `minAss(𝔞 : 𝔟^∞)` | `minAss`, `satByIdeal` |
| 1738–1760 | `K_i`, the discarded-prime rule, `𝔥_i ⊇ h_i` | see finding (F1) below |
| 1786–1800 | `Coeffs(r, A)` | `.coeffs` |
| 1815–1822 | minimal primes ↔ irreducible components | `minAss_spec` (→ `Ideal.minimalPrimes`) |
| 1845 | `V_∀ ⊆ V_{∃∀} ⊆ V_∃` | `VforallAssembled_subset_VexistsForall` (**proved**), `VexistsForall_subset_Vexists` (`sorry`) |
| 1855–1865 | `Z(W_i) = 𝒱_i ∩ Z(C_i) ∩ N_i` | `Zp_Wi_eq` (`sorry`) |
| 1866–1872 | `V_∀ = ⋂(Z(W_i) ∪ ∁Z(C_i))`, `V_{∃∀} = ⋃ Z(W_i)` | `VforallAssembled`, `VexistsForall`, `mem_VforallAssembled_iff` (**proved**) |
| 1876–1911 | the four-subroutine table | `DifferentialDecompose`, `FullReduce`, `minAss`, `RefiningPartition` |
| 1919–1943 | Lemma `lem:specialization` | `specialization_simple` |
| 1946–1994 | Thm. `thm:completeness` | `ritt_remainder_zero` (§9, from `Completeness.lean`) |
| 2093–2110 | Cor. `cor:assembly` | `mem_VforallAssembled_iff` (**proved**) |

## Axiom inventory

Eighteen axioms.  "Read?" says whether the source's own text was read directly
during this formalization.

| Lean axiom | What it quotes | Source, with the source's numbering | Read? |
|---|---|---|---|
| `radical_isDiffIdeal` | radical of a differential ideal is differential (char 0) | Ritt 1950, Ch. I; Kolchin, *DAAG* §I.2 | Ritt Ch. I read via `djvutxt`; statement is standard |
| `projection_eq_cell` | projection of `Sol(S)` onto `ℂⁿ` is the cell | ThomasDecomp **Remark 2.3** (p. 3) | **yes**, verbatim |
| `specialization_simple` | initials/separants survive specialization at a cell point | NewMethod Lemma `lem:specialization` (tex 1919) | yes (the paper) |
| `DifferentialDecompose` | the subroutine | ThomasDecomp **Alg. 3.6** (p. 22) | **yes**, verbatim |
| `differentialDecompose_spec` | its Output clause | ThomasDecomp Alg. 3.6 + **Def. 2.4** | **yes** |
| `FullReduce` | the subroutine | Ritt 1950 **§I.6** ("REDUCTION", p. 5) | **yes**, verbatim |
| `fullReduce_spec` | `h·P − r ∈ [S]`, `r` reduced | Ritt 1950 §I.6, p. 5 | **yes**, verbatim |
| `minAss` | the subroutine | GTZ 1988 **§9** ("Applications to Computing Radicals and Associated Primes") | **yes** — see finding (F3) |
| `minAss_spec` | `minAss 𝔞 = 𝔞.minimalPrimes` | GTZ §9; ties to Mathlib `Ideal.minimalPrimes` | **yes** |
| `RefiningPartition` | the subroutine | CSTools 2008 **§2.4, Example 8** | **yes** — see finding (F4) |
| `refiningPartition_disjoint` | blocks pairwise disjoint | CSTools §2.4 Ex. 8 | **yes** |
| `refiningPartition_cover` | each family member is a union of blocks | CSTools §2.4 Ex. 8 | **yes** |
| `refiningPartition_index` | `ι` records containment | NewMethod tex 1897–1902 (stronger than CSTools) | yes — departure flagged |
| `simple_membership` | `E : q^∞` radical; membership ⟺ pseudo-remainder 0 | Robertz 2018 **Prop. 3.31** = [Rob14] **Prop. 2.2.50** | **yes**, verbatim |
| `decomposition_membership` | `√(E:q^∞) = ⋂ᵢ (E^{(i)}:(q^{(i)})^∞)` | Robertz 2018 **Prop. 3.32** = [Rob14] **Prop. 2.2.72** | **yes**, verbatim |
| `rittReduce` | Ritt reduction (Completeness.lean form) | Ritt 1950 §I.6 | **yes** |
| `rosenfeld` | Rosenfeld's Lemma | BLOP 2009 **Theorem 23** | via `~/project/Completeness.lean` |
| `ritt_remainder_zero` | Thm. `thm:completeness`, kernel form | NewMethod tex 1946; `~/project/Completeness.lean` | yes |

Local paths of the sources read:

* ThomasDecomp — `~/project/papers/bachler-gerdt-lange-hegermann-robertz-2012-differential-thomas-decomposition.pdf`
* Ritt 1950 — `~baccala/Books/Differential Algebra/Ritt J.F. Differential algebra (AMS, 1950)(T)(189s).djvu` (via `~/miniforge3/envs/djvu/bin/djvutxt`; the `.pdf` copy is a scan with no text layer)
* GTZ 1988 — `~baccala/Books/Algebra/Primary Decomposition/GianniTragerZacharias.pdf`
* CSTools 2008 — `~/project/papers/chen-lemaire-li-moreno-maza-pan-xie-2008-constructiblesettools-regularchains.pdf`
* CTD 2007 — `~/project/papers/chen-golubitsky-lemaire-moreno-maza-pan-2007-comprehensive-triangular-decomposition.pdf`
* Robertz — `~/project/papers/robertz-2018-formal-methods-systems-pdes-cours-cirm-published.pdf`
* BLOP 2009 — `~/project/papers/boulier-lazard-ollivier-petitot-2009-radical-differential-ideal.pdf`
* Kolchin — `~baccala/Books/Differential Algebra/Kolchin - DAAG.pdf`

## `sorry` inventory

Eight, classified as **(a)** definitional stub, **(b)** cited theorem (true,
proved in the literature or in the paper, not mechanized here), **(c)** genuine
gap (not known to be true, or known false as stated).

| Line | Declaration | Kind | What would have to be proved |
|---|---|---|---|
| 1083 | `consistencyLocus_correct` | (b) | `Z(consistencyLocus I) = V_∃`.  Needs `projection_eq_cell` plus a differential Nullstellensatz linking `1 ∉ [·] : Q^∞` to existence of a solution — the link the paper deliberately declines to assert (tex 1381–1384). |
| 1272 | `VexistsForall_subset_Vexists` | (b) | `V_{∃∀} ⊆ V_∃`.  Asserted at tex 1845, explained at 1826–1836; no separate proof in the paper. |
| 1297 | `membershipLocus_forAll_correct` | (b) | `Z(membershipLocus I ∀) = V_∀`.  Cor. `cor:assembly` + the `RefiningPartition` contract + `Z(𝔞,𝔟) = ⋃_{𝔭 ∈ minAss(𝔞:𝔟^∞)} Z(𝔭,𝔟)`; the last needs `V(𝔞) = ⋃ V(𝔭)` over `ℚ[c] → ℂⁿ`. |
| 1310 | `Zp_Wi_eq` | (b) | `Z(W_i) = 𝒱_i ∩ Z(C_i) ∩ N_i` (tex 1855–1865). |
| 1373 | `TriangularSet.PartiallyReduced` | (a) | "no derivative of a leader of `A` occurs in `r`" — needs a jet-variable / degree apparatus Mathlib lacks. Carried over from `Completeness.lean`. |
| 1377 | `TriangularSet.FullyReduced` | (a) | as above, plus the degree condition. |
| 1381 | `TriangularSet.Regular` | (a) | BLOP Definition 22 (differentially triangular + separants partially reduced + coherent critical pairs). |
| 1425 | `lazard_gap` | **(c)** | *The paper's own known gap.*  The published argument claims fully-reduced monomials form a `K(N)`-basis of each `Fᵢ = Frac(R/𝔭ᵢ)`; they span but need not be independent.  Counterexample `A = {u³ − u²}`, `r = u − 1` (fails regularity; whether regularity restores the conclusion is open).  See `~/project/docs/completeness-theorem-split.md`. |

Three further notions in §3 — `Involutive` (S2), `MinimalEqs` (S3),
`IneqsIrreducible` (S4) and `SquareFreeAt` (A3) — are stated as `True` rather
than `sorry`, so they raise no warning, but they are **the same kind of debt as
(a)** and are labelled `NOTE (definitional stub)` in the file.  They are
`Prop`-valued placeholders: `DifferentiallySimple` therefore currently asserts
only (S1) with real content.  Stating them properly needs Janet division's
multiplicative split and ThomasDecomp's `Reduce` (Alg. 3.3).

`sorry2lemma.py` was run against the file.  It found exactly the five *tactic*
sorries (the four §8 theorems plus `lazard_gap`) and lifted each into a
`…​.sorried` top-level lemma — but since those five are already top-level named
obligations with fully explicit hypotheses, the lift was a no-op in substance,
and the tool's transformed output does not compile (it re-declares the universe
name `uu` in §9).  The transformation was therefore **not** applied.

## What Mathlib provided vs. what was built here

| Notion | Mathlib | Built here |
|---|---|---|
| differential rings | nothing usable (`Derivation` is single, module-valued, and does not give commuting families over a polynomial ring) | `IsDerivations`, and `DifferentialRing` carried over from `Completeness.lean` |
| differential ideals `[S]` | nothing | `IsDiffIdeal`, `diffIdealBy` (**a real definition**, replacing the `sorry` stub in `Completeness.lean`), with `diffIdealBy_le` / `subset_diffIdealBy` / `diffIdealBy_isDiffIdeal` proved |
| saturation `I : q^∞` | `Submodule.colon` only (single-step) | `satBy`, `satByIdeal`, `colonIdeal`; `satBy_satBy`, `satBy_radical_satBy`, `satBy_isDiffIdeal`, `deriv_pow_mul` all **proved** |
| radical | `Ideal.radical`, `Ideal.le_radical` — used directly | — |
| minimal primes | `Ideal.minimalPrimes` — used directly; `minimalPrimes_top` proved here | — |
| `V(·)` and the ideal–variety correspondence | `MvPolynomial.aeval` and `RingHom.ker`; **no** `V(·)`/`Z(·)` for `ℚ`-ideals evaluated over `ℂ` | `Vz`, `Zp`, `Zw`, `evalHom`; `Vz_span`, `Vz_top`, `Vz_one`, `Vz_bot`, `Vz_antitone`, `Vz_mul`, `Vz_listProd`, `Vz_sup` all **proved** |
| constructible sets | nothing (Mathlib has `IsConstructible` for topological spaces, but not the CTD ideal-pair representation) | `IsConstructible` after CTD Def. 6 |
| simple systems, Thomas decomposition, Ritt reduction | nothing | axiomatized from the sources |

## Findings — discrepancies noticed while translating

**(F1) `𝔥_i`: the pseudocode and the prose are different computations.**
Algorithm 2 line 11 (tex 1693) sets
`𝔥_i ← ∏_{q ∈ S_i^≠} ⟨Coeffs(q, ℚ[c])⟩`.  The prose at tex 1738–1750 instead
defines `K_i` — the union over `q ∈ S_i^≠` of the *minimal associated primes*
of `⟨Coeffs(q, ℚ[c])⟩` — sets `𝔥_i = ∏_{𝔮 ∈ K_i} 𝔮`, and adds a rule
discarding any `𝔭 ∈ minAss` that contains some `𝔮 ∈ K_i`.  Neither `K_i` nor
the discarding rule appears in the pseudocode.  `Z(W_i)` is unaffected
(`V(𝔞) = V(√𝔞) = ⋃ V(𝔭)`), but the two are not the same computation, and the
`K_i`-based discarding step is genuinely absent from the algorithm.  The Lean
follows the pseudocode (`hfrak`, with the discrepancy recorded in its
docstring).

**(F2) the paper's definition of saturation is malformed.**  tex 1087 reads
`I : S^∞ = {p : ∃ q ∈ R, q p ∈ I}`.  The multiplier ranges over all of `R`
rather than over powers of `S`, and `S` does not occur on the right-hand side;
taken literally the set is all of `R`.  Every use in the paper is the standard
saturation, which is what `satBy` formalizes.  (The same passage calls a
differential ideal a "subring", which it is not.)

**(F3) `minAss` is cited to GTZ Cor. 3.2(v), which is about saturation.**
tex 1893–1896 attributes `minAss` to "\cite{GTZ}, Cor. 3.2(v), §9".  GTZ
Corollary 3.2(v) states that `I·R[x]_f ∩ R[x]` is computable for a
nonzerodivisor `f` — i.e. it is the *saturation* primitive `I : f^∞`.  The
minimal-primes content is §9.  Cor. 3.2(v) is exactly the right citation for
the `𝔞 : 𝔟^∞` of Algorithm 2 line 17, so the two citations look to have been
attached to the wrong subroutine.

**(F4) `RefiningPartition` is cited to CSTools §3; it is in §2.4.**  §3 of that
paper is the `ParametricSystemTools` module; `RefiningPartition` is documented
in §2.4 (`ConstructibleSetTools`), Example 8.  The paper's tabulated
description is also *stronger* than CSTools' own: CSTools claims a pairwise
disjoint refinement with provenance indices, while the paper additionally
requires `ι` to be exactly the set of family members containing `B`.  Algorithm 2
line 16 needs the stronger reading, so it is what `refiningPartition_index`
formalizes (departure flagged).

**(F5) `V_∀ ⊆ V_{∃∀}` needs the covering hypothesis, and cannot be stated
without it.**  The paper attaches the hypothesis "whenever the ansatz is
consistent throughout `ℂⁿ`" (tex 1675–1677) but does not remark that it is
*necessary*.  Formalizing made this unavoidable: a point lying in no cell
belongs to `⋂ᵢ (Z(W_i) ∪ ∁Z(C_i))` vacuously but to no `Z(W_i)`, so the
containment fails outright without it.  See
`VforallAssembled_subset_VexistsForall`, where the hypothesis is an explicit
argument.

## Relationship to `~/project/Completeness.lean`

§9 reproduces that file's `DifferentialRing`, `TriangularSet`,
`sepMonoid`/`algIdeal`/`diffIdealOf`, `satByMonoid` (there `Ideal.satBy`),
`rittReduce`, `rosenfeld`, `lazard_gap`, and `ritt_remainder_zero`, keeping the
names so the two developments can be read together.  Two deliberate changes:

* `diffIdeal` is no longer a `sorry` — §1's `diffIdealBy` gives it a real
  definition, and §9's `diffIdeal` is the same `sInf` construction phrased for
  the `DifferentialRing` class.
* `Completeness.lean`'s `lazard` axiom is **not** carried over.  It quotes BLOP
  Theorem 3 in terms of `nonLeaderRing` and `nonLeaderDim`, two further
  definitional stubs that nothing in *this* file uses; `lazard_gap` — the actual
  obligation — is carried over.

`ritt_remainder_zero` is what makes `membershipLocus`'s inner loop (lines 5–6)
correct: for `c*` in the cell, `c* ∈ V(Coeffs(r, ℚ[c]))` iff
`P ∈ [S(c*)] : q(c*)^∞`, i.e. `J_i` cuts out exactly the constants at which
`P_j` is a member.  It is restated, not re-proved; see `~/project/Completeness.lean`
and `~/project/docs/completeness-theorem-split.md`.
