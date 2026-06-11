# NewMethod and the Parametric Decomposition Landscape: a Comparison

**Date:** 2026-06-11
**Context:** Written after a working session (2026-06-10/11) that walked the
full chain: GTZ ↔ CTD ↔ comprehensive Gröbner systems on the commutative
side, then Rosenfeld–Gröbner, parametric RG, and differential Thomas
decomposition on the differential side, ending with where the NewMethod
algorithm (Algorithm 2 / Theorem 3 of `NewMethod.tex`, work in progress)
sits among them. Companion to the litsearch report
`~/project/reports/newmethod-parametric-differential-elimination-litsearch.md`
(2026-06-04), which carries the full prior-art citations and verdicts; this
report is the *conceptual* comparison, with the worked examples from the
session.

---

## TL;DR

All the methods surveyed here describe the same solution set; they differ in
**which question about it they answer** and **what triggers a case split**.

- **GTZ / minimal primes** decompose the solution set itself — intrinsic,
  canonical, projection-free. No parametric information beyond what creates
  whole components.
- **CTD, comprehensive Gröbner systems, Thomas decomposition** (algebraic or
  differential) answer the *fibered* question: stratify parameter space so
  that above each cell the solutions have one uniform description. They
  split wherever the **description** degenerates.
- **NewMethod** answers the *projection* question: compute the set of
  constants `c*` at which the ansatz solves the PDE — the **image** of the
  solving set in constants space, returned as canonical minimal primes. It
  splits only where **Theorem 3's hypotheses fail** (the bad locus), not
  wherever the description changes. In the good case (coherent ansatz,
  empty bad locus — hydrogen is one) it performs *no* differential case
  analysis at all: one Ritt reduction over `ℚ(c)`, then one commutative
  minimal-prime computation.

Capability-wise, the differential Thomas decomposition can reach the same
answer (litsearch verdict: position-and-cite). Route-wise, nothing else
computes it this way, and the route is where the efficiency lives: the
constants stay in the coefficient field, out of the differential engine.

---

## 1. The commutative baseline

### 1.1 GTZ versus the fibered question

Given `F ⊂ K[u, x]` (parameters `u`, unknowns `x`), GTZ-style minimal-prime
decomposition treats all variables on an equal footing and returns the
irreducible components of the total variety `V(F) ⊆ K̄^{d+m}` — canonical,
intrinsic, no projection chosen.

The fibered question — *for each parameter point `u₀`, what does the fiber
`V(F(u₀, x))` look like?* — is **not** answered by that decomposition,
because an irreducible variety can fiber non-uniformly. The minimal example
is `⟨ax − b⟩ ⊂ ℚ[a,b,x]`: it is **prime** (one component, nothing to
decompose), yet its fibers over the `(a,b)`-plane have three regimes
(`a ≠ 0`: one point; `a = 0, b ≠ 0`: empty; `a = b = 0`: a line). The fiber
dimension jumps *inside* a single irreducible component, invisibly to any
decomposition of the total space. GTZ sees parameter-special behaviour only
when it creates a separate component (`⟨ux⟩ = ⟨u⟩ ∩ ⟨x⟩` — the vertical
sheet `u = 0` appears as its own prime).

Recovering the fiber structure from the minimal primes requires, per
component: projection closure (elimination), generic fiber dimension, the
dimension-jump/collision loci via iterated resultants of initials and
discriminants, and recursion on those loci — i.e., re-deriving exactly the
discriminant-variety machinery that CTD packages. In
algebraic-geometry terms, CTD is roughly a computable **flattening
stratification** of the family `V(F) → K^d`; the prime decomposition is a
description of the total space with the map forgotten.

Conversely, the fibered methods hand each cell a per-parameter *solver*
(specialize the chain, back-substitute); the minimal primes give no such
thing.

### 1.2 Comprehensive Triangular Decomposition (CTD)

Chen–Golubitsky–Lemaire–Moreno Maza–Pan, CASC 2007 (LNCS 4770, 73–101);
Hierarchical CTD: Chen–Tang–Xia 2014. Both archived in
`~/project/papers/`.

Output: a finite partition of `Π_U(V(F))` into constructible cells, each
carrying a family of regular chains that **specializes well** (initials
nonvanishing, squarefreeness surviving) at *every* point of the cell, the
specialized quasi-components uniting to the exact fiber. Construction: the
generic decomposition is valid off a **discriminant variety** (cut out by
iterated resultants of initials/discriminants, eliminated to the
parameters — Lazard–Rouillier); on that lower-dimensional bad locus the
construction **recurses**, terminating by strict parameter-dimension
descent.

The paper's own Example 1 — `F = {vxy + ux² + x, uy² + x²}`, parameters
`u > v`, unknowns `x > y` — shows the behaviour on a **reducible** variety:
a global pool of four chains, each cell selecting a subset:

| Cell | Chains | Fiber |
|---|---|---|
| `C₁: u(u³+v²) ≠ 0` | `{T₃, T₄}` | origin + two points |
| `C₂: u = 0` | `{T₂, T₃}` | origin + the line `x = 0` (dimension jump) |
| `C₃: u³+v² = 0, v ≠ 0` | `{T₁, T₃}` | origin + one point (roots merged) |

with `T₁ = {vxy+x−u²y², 2vy+1, u³+v²}`, `T₂ = {x, u}`, `T₃ = {x, y}`,
`T₄ = {vxy+x−u²y², u³y²+v²y²+2vy+1}`. Components persist across cells
(`T₃`), collide on discriminant strata (`T₄ → T₁`), inflate in dimension on
degenerate strata (`T₂`), and die where their fibers go empty. The cell
partition is the common refinement, over the chain pool, of the
specialize-well, nonemptiness, and collision loci.

Caveats: the chains are quasi-components — equidimensional, not
irreducible (the single chain `{x² − uy²}` specializes well for `u ≠ 0` yet
every such fiber is two lines); and reducible inputs are where the case
count multiplies, the motivation for the Hierarchical variant.

### 1.3 Comprehensive Gröbner systems and the Gröbner Cover

Same fibered question, Gröbner-school answer. Distinguish **CGB**
(Weispfenning 1992: one parametric set that is a Gröbner basis under every
specialization — no triangular analogue exists) from **CGS** (pairs
`(Cᵢ, Gᵢ)`, a basis valid throughout each cell — the true CTD mirror). The
**Gröbner Cover** (Montes–Wibmer 2010) is the canonical CGS: locally closed
segments with distinct leading-power-product structure, reduced bases —
uniquely determined by ideal and term order, unlike any CTD cell structure.

Both schools split where their representation's *leading structure* breaks:
leading coefficients (GB) versus initials+squarefreeness (chains). Per-cell
deliverables differ: uniform Hilbert function / membership test (GB) versus
back-substitution solving and equidimensional fiber decomposition (chains).
A CGS is the most scheme-faithful fiberwise (it carries the honest
specialized ideal, multiplicities included); CTD is set-theoretic with
squarefreeness enforced; minimal primes describe only the radical.

Singular's `grobcov.lib` computes the Gröbner Cover and is the one
comprehensive tool runnable locally. Session demos:

- CTD Example 1 → **5 segments**: the three CTD cells with `v = 0` split
  off the generic cell and the origin split off `u = 0` (canonicity forces
  splits at every lpp change, finer than specialize-well requires).
- `ax − b`, split `(a,b | x)` → 3 segments (basis `{ax−b}` / `{1}` / `{0}`).
- `ax − b`, split `(a,x | b)` → **one segment covering everything**
  (basis `{b − ax}`): the same prime ideal, fibered along a different
  projection, has no degeneration at all. Degeneracy is a property of the
  projection, not the variety; the parameter/unknown split is mathematically
  free (any subset works; the well-posed choice is a transcendence basis of
  dimension `d = dim V(F)`), and the total variety reassembles identically
  from any split's cells. GTZ is what you compute when you refuse to choose.

### 1.4 Algebraic Thomas decomposition

A third commutative school (J.M. Thomas 1937; Wang ~1998; the modern
algorithm in Bächler–Gerdt–Lange-Hegermann–Robertz 2012, §2). Decomposes
any constructible set into **disjoint** *simple systems* — triangular,
nonvanishing initials, squarefree, with inequations as first-class
citizens. Disjointness (no other triangular method has it) is bought with
extra splitting and pays for exact solution counting (Plesken's counting
polynomial) and automatic generic/special separation:

| | GTZ / minimal primes | Regular chains (Kalkbrener/Lazard) | Algebraic Thomas |
|---|---|---|---|
| Pieces | irreducible components | quasi-components (may overlap) | simple systems (disjoint) |
| Canonical? | yes | no | no |
| Inequations native? | no | bolted on | yes |
| Counting per piece | no | partial | yes, by inspection |

---

## 2. The differential side: what breaks

### 2.1 No differential GTZ

The existence theorem survives: Ritt–Raudenbush gives every radical
differential ideal `{F}` a unique finite intersection of minimal prime
differential ideals. Coverage survives too: Rosenfeld–Gröbner decomposes
`{F}` into characterizable components whose solutions union to the *entire*
solution set, singular solutions included; Ritt's classical refinement
(with factorization) even makes the components prime.

What is missing is **minimality**: deciding containment between components
is the **Ritt problem**, open since the 1930s. Ritt's Low Power Theorem
settles special shapes (`(y′)² − 4y`: the singular component `y = 0` is
essential), but no general algorithm prunes redundant components. Further
degradations: differential ideals are not finitely generated as ideals
(only radical ones, as radical ideals), there is no terminating
differential Gröbner basis in general, membership is decidable only at the
radical — so primary-decomposition-style scheme structure has no
differential counterpart at all. And the implementations (BLAD /
`DifferentialAlgebra`, Maple) are factorization-free for efficiency: their
components are characterizable (radical, possibly reducible), not prime —
Hubert's *Notes II* (LNCS 2630) is the standard treatment of the
characterizable-vs-prime distinction.

Net: on the commutative side both schools are effective (GTZ and
triangular/comprehensive); on the differential side **only the triangular
school survives**. This asymmetry is why NewMethod's strategy — project to
the commutative side and run honest minimal-prime machinery there — does
real work.

### 2.2 Constants are the only fiberable directions

A parametric specialization `u ↦ c` must commute with the derivations to
mean anything differentially. For a differential indeterminate (`δu ≠ 0`)
it cannot: the differential ideal `[u − c]` contains `δ(u−c) = u′, u″, …`,
so the "fiber" `u = c` is not a slice of solution space but the locus *`u`
is identically constant* — the slice collapses rather than restricts.
Geometrically, the total-derivative flow in jet space runs *through* the
slices `{u = c}`; only constant directions are transverse the way a
parameter axis needs to be. For a constant indeterminate (`δa = 0`),
`[a − c] = (a − c)` and fibers behave exactly as in the commutative case.

Consequences: the commutative "any split works" symmetry is broken — the
legitimate parameters are exactly (a subset of) the **constants**, which
conveniently form a field; the residual freedom is only which constants are
parameters versus unknowns (`E` versus the `cᵢ`). The escape hatch —
specializing a non-constant to a *function* — is a valid differential
specialization, but the base becomes a function space and the
constructible-cell technology has no home there. This is also why plain
Rosenfeld–Gröbner has no parameter space to stratify: its splits (initials,
separants) are regularity conditions on the *jet* variables, partitioning
the solution set, not a parameter space; only after designating constants
does the comprehensive question become well-posed.

---

## 3. The differential comprehensive methods

### 3.1 Parametric Rosenfeld–Gröbner (Fakouri–Rahmany–Basiri 2018)

Cogent Math. & Stat. 5(1):1507131. Mechanism: run RG with a comprehensive
Gröbner system maintained on the parameter conditions — output "regular
representations for all possible states of the parameters": a case split of
constant-parameter space, regular differential chains per cell. So it
stratifies the constants the way a *single* comprehensive pass does. Per
the litsearch assessment it does **not** confirmably provide CTD's two
deeper guarantees: squarefree specialization at every special value in a
cell, and termination via the iterated discriminant-variety recursion with
dimension descent — the ingredients that make a cover *comprehensive*
rather than just case-split, and exactly what NewMethod's recursion adds.

Project empirics (port: `~/parametric-rg`, task 384, report
`~/project/reports/parametric-rg-port.md`): all three of the paper's
examples reproduce; on the hydrogen ansatz the **parameter branching is
cheap and clean** (worklist stays at 3–4 vertices) while the **differential
reduction blows up** (per-vertex coefficient swell, the same BLAD
factorwise-reduction wall) — the run is time-boxed, not branch-exploded.
The cross-check that matters: the parametric-RG cells on which the
Schrödinger PDE becomes redundant modulo the ansatz coincide **exactly**
with the five minimal associated primes computed by `joca.sage` — the two
routes agree where both can reach.

### 3.2 Differential Thomas decomposition (Bächler–Gerdt–Lange-Hegermann–Robertz 2012)

J. Symbolic Comput. 47(10):1233–1266 (archived). The algebraic Thomas
machinery of §1.4 with two subroutines made differential. A *differential
simple system* is algebraically simple in its jet variables **and
Janet-involutive** (passive): all non-multiplicative prolongations reduce
to zero, i.e. all integrability conditions are accounted for.

Two ingredients matter for the comparison:

- **Squarefreeness is load-bearing**: the initial of every nontrivial
  prolongation is the *separant*, and a squarefree polynomial shares no
  roots with its separant — so the algebraic simplicity conditions are
  precisely what licenses reduction modulo all prolongations. The separant
  loci that Ritt's theory treats as lost singular strata ride along as
  explicit inequations.
- **Janet division** assigns each leader a private cone of multiplicative
  derivations (scan derivations in order; `∂_l` is multiplicative for `w`
  iff `w` has maximal `l`-th exponent among generators agreeing in the
  first `l−1` exponents). Cones are disjoint, every principal derivative
  has a unique reductor (deterministic reduction — no strategy choices),
  the finitely many non-multiplicative prolongations are the integrability
  conditions (replacing RG's pairwise Δ-polynomials), and the disjoint
  cones form a Stanley decomposition from which solution-counting (the
  differential counting polynomial, Lange-Hegermann 2014) reads off
  directly. The complementary free derivatives are Riquier's initial data.

Parameters: model each constant as `ċᵢ = 0` (legitimate per §2.2); the
squarefreeness/initial inequations *in the constants* are then
parameter-degeneracy conditions, and disjointness forces automatic
generic/special separation (demonstrated in the parametric-control
companion paper, Lange-Hegermann–Robertz 2013). This makes differential
Thomas **the one published method with a CTD-grade output guarantee
differentially** — but reached *lazily*, by local inequation splits: no
discriminant variety is ever named as a parameter-space object, and
termination is a Noetherian ranking argument, not parameter-dimension
descent. RG-school tools are stronger on ODE systems in practice; Thomas
pays its disjointness tax in finer splitting and buys counting and clean
specialization. All implementations (Wang's Epsilon; Aachen's
`AlgebraicThomas`/`DifferentialThomas`) are Maple-hosted.

---

## 4. NewMethod

### 4.1 What it is

**Core loop** (Algorithm 1, `NewMethod.tex` §"Solution Algorithm"):
Ritt-reduce the PDE `P` modulo the ansatz `A` over `ℚ(c)` — constants in
the **coefficient field**, never ring variables — obtaining
`h·P − r ∈ [A]`, `r` fully reduced, `h` a power product of initials and
separants. Collect `r = Σ m_α p_α(c)` by power products in the
non-constants; the candidate solution locus is `𝒱 = V(p_1, …, p_N) ⊆ ℂⁿ`,
decomposed into **minimal primes by commutative GTZ** (the entire
helium/Singular pipeline applies).

**Theorem 3** (completeness): if `c*` solves (`P ∈ [A(c*)]`), the
denominator does not vanish (`h(c*) ≠ 0`), `A` is a parametric regular
differential system, and `A(c*)` is a squarefree regular chain, then
`c* ∈ 𝒱`. Its specialization lemma is explicitly the differential analogue
of CTD's "specializes well" (the paper cites CTD Def. 11/Prop. 4); the
squarefree half of the hypotheses need not specialize — that is what the
bad locus tracks.

**General algorithm** (Algorithm 2): RG over `ℚ(c)` → squarefree regular
chains generically; core loop on each; compute the **bad locus** `B` (where
specialization to a squarefree regular chain fails — the differential
discriminant variety); recurse over the function fields `ℂ(Bʲ)` of its
irreducible components; terminate by strict dimension descent. Corollary:
completeness for an arbitrary ansatz, unconditional but for the denominator
locus.

### 4.2 The two structural differences from everything in §§1–3

1. **Image, not atlas.** Every comprehensive method outputs a
   stratification of parameter space with per-cell solving descriptions —
   the structure of the map `V → K^d`. NewMethod outputs the **projection
   of the solving set to constants space**: where solutions *exist*, as a
   union of irreducible varieties. A smaller object, and — as minimal
   primes — a **canonical** one, which no comprehensive method's cell
   structure is. (The natural cost: no per-cell solution descriptions; the
   eigenfunctions are recovered by re-instantiating the ansatz at `c*`. A
   deliberate non-goal, worth stating explicitly in the paper.)
2. **Splits indexed by theorem failure, not description change.** Thomas
   and CTD split wherever the description's leading structure degenerates,
   whether or not that changes anything about existence; parametric RG
   branches at every mid-elimination parameter ambiguity. NewMethod's case
   tree is the failure set of Theorem 3's hypotheses and nothing else. In
   the good case — coherent ansatz, `B = ∅`; the hydrogen ansatz is exactly
   this — there are **zero splits**: one Ritt reduction, one GTZ call.

### 4.3 Head-to-head

| | Parametric RG (Fakouri) | Differential Thomas | NewMethod |
|---|---|---|---|
| Output | regular representations per parameter cell | disjoint simple systems, all solutions | solution locus in constants space (minimal primes) |
| Splits when… | parameter status ambiguous mid-elimination | description non-uniform (squarefree/initial inequations) | Theorem-3 hypotheses fail (bad locus only) |
| Constants live… | branched on explicitly | indeterminates with `ċ = 0` | coefficient field `ℚ(c)` |
| Squarefree at special values | unconfirmed | yes (definitional) | yes (forced by recursion) |
| Termination proof | (unconfirmed) | Noetherian ranking on systems | parameter-dimension descent |
| Generic-nice-case cost | full branch tree | full disjoint cover | 1 reduction + 1 commutative GTZ |
| Canonical output? | no | no | yes (minimal primes) |
| To get the existence locus | post-filter cells where `P` reduces to 0 | post-filter | it *is* the output |

The practical lever underneath the table: **keeping the constants out of
the differential engine**. Direct RG with constants as ring variables does
not terminate on hydrogen (BLAD coefficient swell — measured, task 384);
Thomas carries the constants as `ċ = 0` indeterminates through every
involutive reduction; Fakouri branches on them throughout. NewMethod
touches them differentially not at all in the generic pass.

Relationship to the others (litsearch verdict, 2026-06-04):
position-and-cite. The *capability* — a comprehensive cover valid at all
parameter values — is achievable via differential Thomas; the *assembly* —
Hubert/RG characteristic decomposition over `ℚ(c)` + an explicit
differential discriminant variety + function-field recursion with dimension
descent — is published only algebraically (CTD 2007, Hierarchical CTD
2014). NewMethod is the differential CTD in the
"specializes-well-off-a-discriminant-variety" form, with the comprehensive
engines (Fakouri, Thomas) cited as the fallback skeleton its Algorithm 2
steps 1/3/4 instantiate.

### 4.4 Open flanks (work-in-progress notes)

- **The denominator locus** `h(c*) = 0` is genuinely outside Theorem 3.
  The paper notes it could be folded into the comprehensive recursion and
  declines (no solutions of interest there for the ansätze considered); a
  referee may push — having the one-paragraph answer ready is cheap.
- **Theorem 3's hypothesis (4)** (squarefree regular chain at `c*`) is the
  load-bearing nontrivial hypothesis; the 8.4-style completeness gap noted
  in the Hubert errata work is adjacent territory — keep the statements
  aligned.
- The claimed output is a *variety*; the existence locus for the weaker
  reading ("some solution of `A(c*)` solves `P`") is in general only
  constructible (see §5.4). The strong reading (`P ∈ [A(c*)]`, every
  solution of the ansatz solves the PDE) is the one Theorem 3 uses;
  keeping the two readings visibly distinct in the paper forestalls a
  natural referee confusion.

---

## 5. Roads not taken (and why)

### 5.1 Direct elimination ranking ("RG alone")

Adjoin `P` to `A`, rank all non-constants above all constants, run RG; the
constants-only equations among the output are the solution-locus
conditions. The classical attack, and the paper's §"Why not
Rosenfeld–Gröbner alone?" addresses it; the project has now also *measured*
it: with the constants as ring indeterminates, BLAD's factorwise reduction
hits coefficient swell and does not terminate on the hydrogen ansatz
(task 384; per-vertex cost 8 s → 48 s over vertices 40 → 50 with the
worklist flat — pure swell, not branching). Wu-style differential
characteristic-set elimination is the same family with weaker guarantees
(no radical membership test, zero decomposition only).

### 5.2 Differential resultants

The "in principle" single-object answer to differential projection: a
resultant whose vanishing characterizes when the system has a common
solution — evaluated on `(P, A)`'s coefficients it would directly cut out a
hypersurface containing the existence locus.

State of the theory: Carrà Ferro gave determinantal differential resultants
for two ODEs ([AAECC 8:539–560, 1997](https://link.springer.com/article/10.1007/s002000050090);
[AAECC 1997, LNCS 1255](https://link.springer.com/chapter/10.1007/3-540-63163-1_5));
Gao–Li–Yuan (2010) gave the first rigorous general definition for `n+1`
differential polynomials in `n` indeterminates; Li–Yuan–Gao defined the
**sparse differential resultant** for Laurent differentially essential
systems with a single-exponential computation algorithm
([Found. Comput. Math. 15(2):451–517, 2015](https://link.springer.com/article/10.1007/s10208-015-9249-9);
[ISSAC 2011](http://mmrc.iss.ac.cn/~weili/papers/ISSAC2011-PAPER.pdf)).

Why it is not a competitor here: the theory covers **generic** (essential)
systems — an ansatz with structured, non-generic coefficients in `ℚ(c)` is
exactly the regime where the classical resultant can vanish identically or
fail to exist; it produces one necessary condition (a hypersurface), not
the full ideal of the locus, and nothing stratified or comprehensive; and
the single-exponential algorithm is not practical at this scale. Verdict: a
related-work sentence ("the resultant-theoretic route to differential
projection exists for generic systems [CarràFerro97, LiYuanGao15] but not
at this generality"), not an alternative.

### 5.3 Comprehensive involutive / characteristic decompositions

Gerdt–Hashemi, *Comprehensive Involutive Systems*
([CASC 2012, LNCS 7442:123–138](https://arxiv.org/abs/1206.0181)) and
*Comprehensive Characteristic Decomposition of Parametric Polynomial
Systems* ([ISSAC 2021, DOI 10.1145/3452143.3465536](https://doi.org/10.1145/3452143.3465536)):
the "comprehensive X" pattern applied to involutive bases and to
Ritt–Wu-style characteristic decomposition — both **algebraic only**. They
confirm the naming convention is standard and well-populated on the
polynomial side and that the differential lift is the open slot; cite for
anchoring, nothing more.

### 5.4 Model theory: why the target object is tame

Differential elimination theory (Seidenberg; model-theoretically,
quantifier elimination for differentially closed fields of characteristic
0, Blum) guarantees that the set of constants `c` satisfying any
first-order condition mixing the ansatz and the PDE — e.g. "every solution
of `A(c)` satisfies `P`" — is a **definable** subset of the constants,
hence (by stable embeddedness of the constant field) a **constructible**
subset of `ℂⁿ` in the ordinary Zariski sense. This is the abstract
guarantee that NewMethod's deliverable is a finite, algebraic object at
all; the algorithm's content is computing it *efficiently and canonically*
(closed form, minimal primes), which QE-style procedures do not provide.

### 5.5 Differential Galois / Kovacic-type algorithms

Kovacic's algorithm (J. Symbolic Comput. 2(1):3–43, 1986) and its
descendants decide Liouvillian solvability for low-order *linear* ODEs —
complete algorithms over a **frozen** solution class. Philosophically the
opposite trade from NewMethod: total completeness for a fixed ansatz
family versus a parameterized-ansatz framework with completeness relative
to the chosen ansatz (Theorem 3). For second-order linear reductions,
Kovacic is also a useful *oracle* for spot-checking NewMethod outputs. A
positioning sentence in related work suffices.

### 5.6 Honourable mentions

- **Thomas decomposition as post-filtered competitor**: run differential
  Thomas with `ċ = 0`, then keep the simple systems in which `P` reduces to
  zero and project their constant-conditions. Capability-equivalent,
  description-driven splitting throughout, Maple-only. This is the honest
  "someone else could compute the same set" baseline the paper should keep
  acknowledging.
- **Numerical parameter continuation** (Bertini/paramotopy-style) — wrong
  category (numeric, no completeness certificate) for the paper's claims,
  though potentially useful for sanity checks on specific ansätze.

---

## 6. Implementations and local infrastructure

| Tool | What | Where | Status here |
|---|---|---|---|
| Singular `grobcov.lib` | Gröbner Cover (canonical CGS) | `~/Singular/install/bin/Singular` | working; session demos above |
| Maple `RegularChains` | regular chains, Triangularize, **CTD** | regularchains.org | **GPL v3 source, release 2548** — downloaded to `~/project/regularchains-upstream/`; site's `download.php` is broken, direct URL works (`…/RegularChains_2548_src.tgz`) |
| RegularChains→Sage port | CTD core + full library, 6 phases | tasks 395–400 (`regularchains-port*`) | created, pending |
| BLAD / `DifferentialAlgebra` | Rosenfeld–Gröbner | sage env + `~/diffalg-venv` | working |
| Fakouri parametric RG | parametric RG port | `~/parametric-rg` (task 384) | working; hydrogen cross-check passed |
| `joca.sage` | NewMethod core loop on hydrogen | `~/Papers/NewMethod/joca.sage` | working; 5 minimal primes |
| Epsilon, `AlgebraicThomas`/`DifferentialThomas` | (algebraic/differential) Thomas | Maple-only | not runnable locally |

---

## References

(Beyond those embedded above.) CTD: Chen, Golubitsky, Lemaire, Moreno Maza,
Pan, CASC 2007, LNCS 4770:73–101. Hierarchical CTD: Chen, Tang, Xia, 2014.
Gröbner Cover: Montes, Wibmer, J. Symbolic Comput. 45 (2010). CGB:
Weispfenning, J. Symbolic Comput. 14:1–29 (1992). Thomas decomposition:
Bächler, Gerdt, Lange-Hegermann, Robertz, J. Symbolic Comput.
47(10):1233–1266 (2012); Lange-Hegermann, Robertz, parametric control
(arXiv:1212.0377); Robertz, LNM 2121 (2014); Lange-Hegermann, differential
counting polynomial (arXiv:1407.5838). Parametric RG: Fakouri, Rahmany,
Basiri, Cogent Math. & Stat. 5(1):1507131 (2018). RG: Boulier, Lazard,
Ollivier, Petitot, AAECC 20 (2009). Characteristic decompositions: Hubert,
*Notes on Triangular Sets and Triangulation-Decomposition Algorithms II*,
LNCS 2630 (2003). Kovacic, J. Symbolic Comput. 2(1):3–43 (1986).
Differential resultants: Carrà Ferro (1997); Gao, Li, Yuan (2010); Li,
Yuan, Gao, Found. Comput. Math. 15(2):451–517 (2015). Comprehensive
involutive/characteristic: Gerdt, Hashemi, CASC 2012; ISSAC 2021,
10.1145/3452143.3465536. QE for DCF₀: Blum (1977); Seidenberg's elimination
theory (1956). Local artifacts: litsearch report (2026-06-04),
`parametric-rg-port` report, task prompts `regularchains-port*`.

---

*Researched and written by Claude on behalf of Brent Baccala, from a
working session of 2026-06-10/11; worked examples (grobcov runs, source
inspections) were executed live during the session.*

/claude-fable-5
