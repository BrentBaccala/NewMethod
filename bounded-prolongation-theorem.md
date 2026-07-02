# Bounded-Prolongation Algebraization — a proof attempt

*Goal.* Replace the differential Thomas decomposition of a finite-type ansatz
by the **algebraic** Thomas decomposition of a finite, bounded-order
prolongation — so the whole computation is pure polynomial algebra (Singular /
RegularChains), with no differential pseudo-remainder. Unlike Theorem 1 of the
paper, this is **saturation-free**: it keeps every cell, including the
bad-locus strata where initials/separants vanish (the ~29 hydrogen cells).

Status: attempt. The skeleton is complete and rests on three classical inputs
(formal integrability of involutive systems; termination of the differential
Thomas algorithm; passive ⟺ involutive). The one quantity left to pin
explicitly is the involutive order `N₀` for this ansatz class; §4 bounds it and
flags it for computational confirmation.

---

## 1. Setup and notation

Work in the differential polynomial ring `R = K{u₁,…,u_m}` with commuting
derivations `Δ = {δ₁,…,δ_d}` (the independent variables). `Θ` is the free
commutative monoid on `Δ`; for `θ ∈ Θ`, `|θ|` is its order. The **jet
variables** are the symbols `θuⱼ`. Fix a ranking `<` on the jets satisfying the
two compatibility conditions of §"Differential Reduction" of the paper
(`u < δu`; `u < v ⇒ δu < δv`).

For `N ∈ ℕ` let
```
   J_N  =  K[ θuⱼ : |θ| ≤ N ]           (a polynomial ring, finitely many vars)
```
be the **order-N jet ring** (adjoin the base variables and the parameters `c`
as needed). The **prolongation of A to order N** is the finite polynomial set
```
   A^[N]  =  { θa : a ∈ A,  θ ∈ Θ,  |θ| + ord(a) ≤ N }  ⊂  J_N,
```
i.e. every derivative of every equation, up to total order `N`. (Throughout,
"order" = differentiation order, **not** polynomial degree — the distinction
matters in §5.)

`DTD(A)` denotes the differential Thomas decomposition of `A`: a finite set of
differential **simple systems** (cells) `(E_i = 0, NE_i ≠ 0)`, with `E_i`
differentially triangular and **passive** (every Δ-polynomial reduces to zero),
`NE_i` the attached initials/separants, the cells' solution sets disjoint and
covering `Sol(A)` (Bächler–Gerdt–Lange-Hegermann–Robertz 2012; Robertz LNM 2121).

`ATD(S)`, for a finite `S ⊂ J_N`, denotes the **algebraic** Thomas decomposition
of `S` in `J_N` (same notion, no derivations), under the ranking obtained by
restricting `<` to jets of order `≤ N`.

Define the **involutive order** `N₀(A)`: the least `N` such that, on each cell of
the running decomposition, the equation set prolonged to order `N` is
**Janet-involutive / passive** — every non-multiplicative prolongation reduces to
zero (Robertz LNM 2121, Def. 2.35; Bächler et al. 2012, Def. 3.5) — equivalently,
the multiple-closed set `S` of principal derivatives has reached its Janet
completion and no integrability condition is born above order `N`. (This is the
Riquier–Janet form of Goldschmidt's involutivity; the two agree under a Riquier
ranking. The ascending chain of `S`'s terminates — Robertz LNM 2121, Lemma 2.9 —
so `N₀` is finite.) Write `R_N = V(A^[N]) ⊆ J_N`.

---

## 2. The theorem

> **Theorem (bounded-prolongation algebraization).**
> Let `A` be a finite-type differential system (the ODE-tower-plus-algebraic
> ansatz), parameters `c`, over `K`. Let `N ≥ N₀(A)`. Then `ATD(A^[N])` and
> `DTD(A)` correspond cell-for-cell: for each differential cell
> `(E_i, NE_i)` of `DTD(A)`, truncating `E_i` and `NE_i` to order `≤ N` yields
> an algebraic cell of `ATD(A^[N])`, and conversely every algebraic cell of
> `ATD(A^[N])` arises this way; the higher-order part of each `E_i` is the
> prolongation of its order-`≤ N` part and introduces no further split.
> No saturation is performed; bad-locus strata appear as their own cells.

> **Corollary (solving the PDE).** Applying the theorem to `A ∪ {P}` (the
> ansatz together with the PDE, each prolonged to order `N`), the solution
> strata of the PDE within the ansatz — including those on the bad locus — are
> exactly the constant-space projections of the cells of `ATD((A∪{P})^[N])`,
> computed entirely by polynomial algorithms.

The corollary is the payoff: it is Theorem 1's conclusion **without**
hypothesis (2) (`h(c*) ≠ 0`), because nothing is divided out.

---

## 2bis. Algorithm-agnostic form (preferred)

The §2 statement is phrased as "differential Thomas decomposition = algebraic
Thomas decomposition," which over-commits: it ties the reduction to one
algorithm. The content actually lives one level lower — at the **ideal/variety**,
not the decomposition. State it there and the choice of algebraic algorithm
(GTZ primary decomposition, triangular / regular-chain decomposition, algebraic
Thomas, or plain Gröbner elimination) becomes a free downstream choice.

> **Theorem (differential→algebraic reduction; algorithm-agnostic).**
> Let `A` (ansatz) and `P` (PDE) be finite type over `K` (the field of the
> constants, or `K[c]` comprehensively), and let `N ≥ N₀(A∪{P})`. Treating every
> jet of order `≤ N` as an independent indeterminate, form the polynomial ideal
> ```
>     J  =  ( (A ∪ {P})^[N] )  ⊆  K[ jets_{≤N} ],     (no saturation, no division)
> ```
> i.e. the equations and all their derivatives up to order `N`. Then:
>
> **(a) Variety = truncated solutions.** `V(J)` is exactly the set of order-`≤N`
> truncations of formal solutions of `A∪{P}`. (`⊆` is trivial; `⊇` is the
> involutivity input — every point extends, §3 Lemma 3 / Lange-Hegermann 2014
> Thm 4.1 + Bächler et al. 2012 §3; elementary 2nd-order-ODE existence for this
> ansatz.)
>
> **(b) Projection to the constants.** With `J_c = J ∩ K[c]` (eliminate the
> jets), `V(J_c)` is the Zariski closure of `{ c : A(c) solves P }`, and the
> solving constants are recovered as the projection `π(V(J))`.
>
> **(c) No stratum lost.** Because `J` is **not saturated** by the initials /
> separants `H_A`, the strata on the bad locus `H_A = 0` survive: there the
> ansatz polynomials degenerate to the correct lower-order constraints (the ODE
> drops to first order), whose zeros stay in `V(J)`. These are exactly the
> components Theorem 1's projection deletes through hypothesis (2).
>
> **(d) Algorithm freedom.** Any algebraic decomposition of `J` (or `J_c`)
> computes the solution structure of the *same* fixed variety `V(J)`: primary
> decomposition (GTZ) → minimal primes = irreducible solution components;
> regular-chain / triangular decomposition (Lazard–Kalkbrener, Wang) →
> quasi-components; algebraic Thomas → disjoint constructible cells; Gröbner →
> membership and elimination. The differential content is fully discharged by the
> prolongation to order `N`; the decomposition algorithm is downstream and free.

**Why this is the right framing.** The single load-bearing differential fact is
(a) — that prolonging to the involutive order `N` makes the *algebraic* variety
`V(J)` coincide with the truncated solution set. Everything after is commutative
algebra on a fixed ideal. In particular the §2 Thomas-vs-Thomas statement is just
the (d)-instance "decompose `J` with algebraic Thomas"; GTZ or regular chains are
equally valid and answer the same question.

**No a-priori bound needed (prolong-until-stable).** One need not compute `N₀` in
advance. Prolong incrementally, `N = 0, 1, 2, …`, and watch the elimination ideal
`J_c^{(N)} = (A∪{P})^[N] ∩ K[c]`. The chain `J_c^{(0)} ⊆ J_c^{(1)} ⊆ …` ascends in
the Noetherian ring `K[c]` and **stabilizes** (Lange-Hegermann 2014, §5.3, the
ascending-chain argument); stop when `J_c^{(N+1)} = J_c^{(N)}`. For finite type
this halts, and the stable ideal is the constant-space answer — computed with
nothing but Gröbner elimination, no involutivity analysis, no Thomas. (`N₀ = 2`
of §4 then just says where it stabilizes for hydrogen.)

**Honest caveat — granularity, not correctness.** The listed algorithms present
`V(J)` at *different granularities*. A minimal-prime decomposition (GTZ) gives the
irreducible components of the solution **set**, and *absorbs* a lower-dimensional
stratum into the closure of a higher one when it is a limit of it; a
Thomas / comprehensive decomposition keeps every stratum as its own disjoint cell
with defining in/equations, singling out the bad-locus pieces. Both represent the
same `V(J)` — they differ only in whether a limit-stratum is displayed separately
or absorbed. So: for **enumerating the solution set**, GTZ on the un-saturated `J`
is the simplest complete choice; for **stratifying** it (which `c` give which
ODE-order), use a constructible (Thomas/comprehensive) decomposition. The
reduction (a)–(c) is identical either way; only the report differs.

---

## 2ter. The algorithm (executable form)

The §2/§2bis theorems say *what is true*; this is *what to run*. The whole point
is that after one bounded prolongation the differential engine is never touched
again — every remaining step is a call to a standard commutative-algebra routine.

> **Algorithm (bounded-prolongation, constructible-split, algebraic-finish).**
> **Input:** a finite-type differential system `A` (the ODE-tower-plus-algebraic
> ansatz), optionally together with a PDE `P`; parameters `c`; base field `K`
> (or `K[c]` for the comprehensive/parametric run).
> **Output:** the solution strata — for `A ∪ {P}`, the solving constants `c`,
> stratified, with no bad-locus stratum lost.
>
> 1. **Prolong to the involutive bound.** Set `N = max(N₀(A∪{P}), ord P)` and
>    form `(A∪{P})^[N]` — every equation and all its derivatives up to total
>    differentiation order `N`. This is the *only* differentiation performed;
>    it discharges all differential content (Lemmas 2–3).
>    *If `N₀` is not known a priori,* use **prolong-until-stable** (§2bis): raise
>    `N = 0, 1, 2, …` and stop when the elimination ideal
>    `J_c^{(N)} = (A∪{P})^[N] ∩ K[c]` stops growing. No involutivity analysis
>    is needed to run this — it is Noetherian termination, not a bound.
>
> 2. **Algebraize.** Treat each jet `θuⱼ` of order `≤ N` as an independent
>    indeterminate. The result is a finite polynomial system in
>    `K[ jets_{≤N}, c ]` — no derivations, no differential pseudo-remainder.
>
> 3. **Constructible split** (optional; = an algebraic Thomas decomposition).
>    Case-split on the distinguished polynomials `= 0` vs `≠ 0`:
>    - **initials** and **separants** of the equations linear in their leader;
>    - **discriminants** (squarefree/gcd conditions) of any equation *nonlinear*
>      in its leader — here the radius relation `r² − x²−y²−z²`.
>    Each branch is a cell `(E = 0, NE ≠ 0)`, `E` triangular, `NE` the retained
>    inequations. **The bad-locus strata are exactly the `separant = 0` (resp.
>    `initial = 0`) branches** — kept as their own cells rather than divided out,
>    which is why no stratum is lost and no saturation is performed.
>
> 4. **Algebraic finish** (§2bis(d), algorithm-free). On each leaf, complete the
>    computation with any sound commutative-algebra routine — they all decompose
>    the *same* fixed variety and differ only in granularity:
>    - **GTZ / primary decomposition** → minimal primes = irreducible solution
>      components (absorbs limit-strata into closures);
>    - **regular chains / triangular decomposition** → quasi-components;
>    - **Gröbner elimination** `J ∩ K[c]` → the solving constants directly;
>    - **algebraic Thomas** → disjoint constructible cells (keeps every stratum).
>    For a GTZ leaf, saturate that *individual* leaf by its own inequations `NE`
>    (an ideal quotient) — this is local to the cell and is **not** the global
>    saturation of Theorem 1 that would delete the bad locus.
>
> **Proof obligation (the one hypothesis that makes step 1's "stop" legitimate).**
> `N₀` must be the involutive order **uniformly over every branch of step 3**,
> since the theorem takes `N₀ = max over cells`. The danger is a degenerate
> branch (e.g. `a₀=a₁=0`) whose altered leader structure has a critical pair with
> least common derivative at order `N₀+1`: stopping at the generic `N₀` would then
> under-constrain *that leaf only*, admitting spurious solutions there. For the
> hydrogen ansatz the claim is that degeneration only *lowers* the order (the ODE
> drops 2nd→1st on the bad locus), so the generic `N₀ = 2` dominates — but this is
> **§5 gap 1**, argued and not yet verified. Two safe ways to honour it: run
> **prolong-until-stable** per branch (step 1's fallback needs no bound), or take
> `N` a hair above the conjectured `N₀` (harmless — extra generators leave `V(J)`
> and hence the decomposition unchanged, §5 item 5).

**Two routes, same theorem.** Step 3 is *optional*. Because every split is on a
distinguished polynomial of order `≤ N₀` (Lemma 4), those polynomials already sit
inside `(A∪{P})^[N]`, so an algebraic algorithm applied to the **single
un-split ideal** rediscovers exactly the same cases. Thus:

- **Explicit-split route** (steps 1–2, then 3, then 4 per leaf) — do the
  branching yourself when you want to *control* the stratification and force the
  bad-locus cells to appear separately.
- **Single-ideal route** (steps 1–2, then straight to 4 on `J = (A∪{P})^[N]`) —
  hand the whole ideal to one decomposition and let GTZ hand back components
  (absorbing limit-strata) or algebraic Thomas hand back cells. Less machinery;
  the right choice when you only need to *enumerate* the solving constants.

Either way the differential work is one prolongation to `N₀`; everything after is
commutative algebra on a fixed polynomial ideal.

### 2ter′. Step-by-step summary (operational)

The formal box above, walked through the way you actually run it for the hydrogen
ansatz, with the choice at each step made explicit.

**Step 0 — pin the involutive bound `N₀` (offline, once).** For this ansatz,
prove it structurally: `N₀ = 2` (§4 — all integrability obligations are the
mixed-partial commutations `Ψ_xy = Ψ_yx`, which close at order 2). If you cannot
pin it in advance, there is no Step 0 — Step 1(b) discovers it.

**Step 1 — reach a bounded involutive prolongation.** Two modes:
- **(a) bound known (the hydrogen case).** Differentiate every equation to total
  order `N₀ = 2` in one shot. **No completion loop, no stabilization test** — the
  proof that `N₀ = 2` is what licenses stopping. This is the whole payoff of
  Step 0.
- **(b) bound unknown (general case).** Run **one** completion mechanism —
  **Janet** involutive completion *or* **Rosenfeld–Gröbner** coherence (pick one;
  they certify the same thing) — as a `prolong → reduce → adjoin-any-nonzero-
  residue` loop. Each nonzero residue is a new integrability condition. Stop when
  a full round produces nothing new (involutive / coherent / passive); the order
  reached is `N₀`.

Either mode, the differential engine is now finished and never touched again.

**Step 2 — algebraize.** Read every jet of order `≤ N₀` as an independent
indeterminate: the system becomes a polynomial ideal `J` in `K[ jets_{≤N₀}, c ]`,
no derivations. **Do not saturate by the initials/separants `H_A`** — that
saturation is exactly what deletes the bad locus (Theorem 1's hypothesis (2)).

**Step 3 — (optional) constructible split.** Case-split on initials / separants /
discriminants (`= 0` vs `≠ 0`) into disjoint cells `(E=0, NE≠0)`. What it buys:
- **skip it** → Step 4 returns the **Zariski closure** of the solving-constant
  set — cheap (one elimination), but may carry spurious lower-dimensional boundary
  constants;
- **do it** → Step 4 returns the **exact constructible** existence set, and the
  **bad-locus strata** (`separant = 0`, e.g. `a₀ = a₁ = 0`) appear as their own
  labeled cells — the concrete evidence that this method keeps what Theorem 1
  drops (the cells carrying `E = −1/2`, `E = −1/8`).

**Step 4 — project onto the constants (existence).** Eliminate all jets:
`J_c = J ∩ K[c]` (per cell if you split in Step 3, else on the whole ideal). Then
`V(J_c)` is the (closure of the) set of constants admitting a solution, and its
**minimal associated primes are the irreducible components** of that set — the
hand-off to `minprimes.sage`. Because Step 2 did not saturate, the bad locus is
included. By Lemma 3 (`N ≥ N₀`), a nonempty fibre over `c` *is* a genuine
truncated solution, so `π(V(J))` is exactly "a solution exists for `c`."

**Decision cheat-sheet.**
- *Just need the region of solving constants* → Steps 1(a), 2, 4; skip 3; accept
  the closure.
- *Need the exact existence set, or want to exhibit the bad-locus solutions* → add
  Step 3.
- *Don't fully trust `N₀ = 2` on a degenerate branch* → run Step 1(b)'s completion
  loop on that branch alone and confirm it stabilizes at order `≤ 2` (§5 gap 1).

**Remark (skipping the split is a superset, never a subset).** Omitting Step 3
costs nothing in completeness. Working saturation-free at `N ≥ N₀`, the
elimination ideal `J_c = J ∩ K[c]` is the **Zariski closure** of the true
solving-constant set `π(V(J))`, hence a **superset**: `π(V(J)) ⊆ V(J_c)`. No
solving constant is ever omitted — the discrepancy `V(J_c) \ π(V(J))` is a set of
*spurious* constants (the elimination ideal vanishes but the fibre is empty, so no
solution exists there), never a missing one. Moreover `π(V(J))` is dense in
`V(J_c)`, so it meets every irreducible component in a dense subset: **every
minimal associated prime of `J_c` is a genuine solution-component** (solving on a
dense open subset), and the spurious constants are confined to a
lower-dimensional sub-locus *within* those components — precisely what Step 3
carves off. Both hypotheses are load-bearing: *no saturation* is what preserves
the bad locus, and `N ≥ N₀` is what makes every point of `V(J)` a real truncated
solution (Lemma 3). Under-prolonging (`N < N₀`) and GTZ-absorption of a
limit-stratum both err the *same* safe way — a superset that may carry spurious
constants but never hides a true solution.

---

## 3. Proof

The argument is four lemmas. Lemmas 1, 2, 4 are elementary; Lemma 3 carries the
analytic content and cites the classical formal-integrability theorem.

### Lemma 1 (reduction coincides at bounded order)

*If `p ∈ J_N` and all reductors lie in `A^[N]`, then differential
pseudo-reduction of `p` modulo `A^[N]` (never prolonging past order `N`) equals
algebraic pseudo-reduction of `p` modulo `A^[N]` in `J_N`, under the matching
ranking.*

**Proof.** The only operation in differential reduction absent from algebraic
reduction is **prolongation** — replacing a reductor `a` by a derivative `θa` to
match a derivative of the leader appearing in `p`. By hypothesis every such `θa`
with `|θ|+ord(a) ≤ N` is already an element of `A^[N]`, so the differential
reducer at each step is literally an element of the algebraic set `A^[N]`.
Modulo that substitution, both procedures are the same sequence of
pseudo-remainder steps in the polynomial ring `J_N`. ∎

### Lemma 2 (integrability conditions are captured)

*For `N ≥ N₀(A)`, every Δ-polynomial of every critical pair of `A` (and of the
cells arising during the differential algorithm), reduced or not, lies in the
ideal `(A^[N]) ⊆ J_N`.*

**Proof.** For a critical pair `{p₁,p₂}` with least common derivative `θ₁₂u`,
the Δ-polynomial `Δ(p₁,p₂) = s₁·(θ₁₂/θ₂)p₂ − s₂·(θ₁₂/θ₁)p₁` is, by inspection,
a `J_N`-linear combination of the derivatives `(θ₁₂/θ_i)p_i`, each of order
`|θ₁₂u| ≤ N` once `N ≥ N₀` (definition of `N₀`: no critical pair has its least
common derivative above order `N₀`). Hence `Δ(p₁,p₂) ∈ (A^[N])`. Reduction only
subtracts further ideal elements, so the reduced form stays in `(A^[N])`. ∎

Lemma 2 is what makes a **single-shot** `ATD(A^[N])` equal to the *iterative*
prolong-and-decompose differential algorithm: the integrability constraints the
differential algorithm discovers one prolongation at a time are already present,
as ideal members of `(A^[N])`, in the input to one algebraic Thomas
decomposition. (The zero set `V(A^[N])` is therefore the same whether or not the
reduced Δ-polynomials are appended, so `ATD` sees them.)

### Lemma 3 (no spurious cells)

*For `N ≥ N₀(A)`, every point of `R_N = V(A^[N])` lying at the generic point of
an algebraic cell of `ATD(A^[N])` is the order-`≤ N` truncation of a formal
power-series solution of the full differential system `A`. Hence each algebraic
cell is the truncation of a genuine differential cell, and no algebraic cell is
an artifact of truncation.*

Three ingredients, the first two now *native* to the Thomas framework rather
than borrowed, the third the only classical existence input.

**(i) The cells are passive by definition, not by a borrowed equivalence.**
A cell of `DTD(A)` is a *differential simple system* (Bächler–Gerdt–Lange-
Hegermann–Robertz 2012, Def. 3.5): algebraically simple **and**
*Janet-involutive* — every non-reductive (non-multiplicative) prolongation of
its equation set reduces to zero modulo that set. Janet-involutivity is exactly
Robertz's *passivity* (LNM 2121, Def. 2.35: `NF(v·bᵢ, T, >) = 0` for every
non-multiplicative `v`), which makes the equation set a **Janet basis** with a
*unique Janet normal form* (LNM 2121, Thm. 2.38(d)). So the
"passive ⟺ involutive" step my first draft borrowed from Goldschmidt is here a
definitional property of `DTD`'s output — nothing to import.

**(ii) Parametric/principal split fixes the data.** A Janet basis partitions the
jet monomials into the multiple-closed set `S` of **principal derivatives**
(leaders of the equations and their multiplicative prolongations) and its
complement `C` of **parametric derivatives** (Janet 1929; Robertz LNM 2121,
§"parametric vs. principal"). The Taylor coefficients of a formal solution at the
parametric derivatives may be assigned **freely**; those at the principal
derivatives are then **uniquely determined** — by the cell's solved-form
equations (orthonomic: each leader expressed through lower data, its initial a
unit by the cell's inequations). For a **finite-type** system `C` is finite, and
the Janet completion of `S` stabilizes at a finite order — this *is* `N₀`.

**(iii) Existence: the assignment extends.** That a consistent assignment of the
parametric Taylor coefficients extends to an actual formal power-series solution
is the **Riquier–Thomas existence theorem** for orthonomic passive systems in a
Riquier ranking (Riquier 1910; Thomas 1937; the nonlinear statement is Bächler
et al. 2012's "every (differentially) simple system has a solution in `E`", `E`
the formal-power-series ring). Convergence to an analytic solution is the
separate Cauchy–Kovalevskaya/Riquier analyticity statement and is not needed
for the cell correspondence.

**Assembling.** Fix `N ≥ N₀`. A point of `R_N` assigns every Taylor coefficient
of order `≤ N`. Because `S` has stabilized by order `N₀ ≤ N`, the constraints on
these coefficients are exactly the cell's solved-form relations expressing the
principal coefficients through the parametric ones — **no constraint among the
order-`≤ N` parametric coefficients is first imposed at order `> N`** (that is
precisely what `N ≥ N₀` buys). Hence the point is a consistent
parametric/principal assignment, and (iii) extends it to a formal solution whose
order-`≤ N` truncation it is. The converse inclusion is immediate: a differential
cell's equations and inequations, truncated to order `≤ N`, are a triangular
algebraic system with the same leaders, i.e. an algebraic cell. So the cells
biject and their solution sets agree. ∎

**Published form (this is essentially a known theorem).** The semantic core of
Lemma 3 — and indeed of the whole §2 theorem — is **Lange-Hegermann (2014),
"The Differential Counting Polynomial," Theorem 4.1**: for an algebraically
restricted system `S` of differential equations and any order `ℓ`,
```
   Sol_ζ(S)_{≤ℓ}  =  ⊔_{S̃ ∈ C}  Sol_ζ(S̃)_{≤ℓ},
```
where `C` is a set of **simple *algebraic* systems** in the order-`≤ℓ`
power-series-coefficient ring `C[G]_{≤ℓ}`, and `Sol_ζ(S)_{≤ℓ}` is the order-`ℓ`
truncation of the genuine formal-power-series solution set. The identification
"derivative ↔ independent algebraic indeterminate" is his explicit bijection
`ρ : {U}_Δ → G(U,Δ)`. His proof (§5.3) is exactly the construction here: apply
`ρ` to the equations *and all their iterated derivatives*, intersect with order
`≤ ℓ`, and run the algebraic Thomas decomposition. Crucially, the step that lets
the **reductive** prolongations to order `ℓ` (i.e. `A^[ℓ]`) stand in for *all*
derivatives —
> "the non-reductive prolongations are redundant, as `S` is involutive (cf.
> [Bächler-Gerdt-Lange-Hegermann-Robertz 2012, §3])"
— is the published form of my Lemmas 2–3, and it is cited there to **Bächler
et al. 2012, §3**. So:

* gap 2′ (existence ⇒ bijection, nonlinear) is closed by Theorem 4.1, which is
  an *equality of solution sets*, not mere existence, valid for nonlinear,
  inequation-bearing systems;
* the involutive-redundancy of non-reductive prolongations (my Lemma 2 + the "no
  new constraint above `N₀`" of Lemma 3) is Bächler et al. 2012 §3, not
  something to reprove.

What is **not** supplied by Theorem 4.1, and is where the work here actually
lives: (a) Lange-Hegermann allows the decomposition `C` to be *countable* (the
inequation tails); finiteness for our ansatz is inherited from Bächler's
*differential* Thomas decomposition terminating with finitely many simple
systems (29 for hydrogen, observed). (b) The **explicit involutive order**
`N₀ = 2` for the ODE-tower ansatz (§4) — Theorem 4.1 is stated per-`ℓ` and does
not by itself say which `ℓ` suffices; the finite-type/single-`v` structure does.
(c) The **saturation-free / bad-locus** reading and the **solving** corollary
(§2) — that running this on `A ∪ {P}` recovers exactly the strata Theorem 1's
projection drops on `h = 0`. These three are the contribution; the core
equivalence is Lange-Hegermann's.

**Specialization to this ansatz (removing the Riquier hypothesis).** The
hydrogen ansatz is finite type with its entire differential structure factoring
through the single variable `v`: every `x,y,z`-derivative of `Ψ` is, by the
chain-rule relations, a function of the `v`-derivatives of `Ψ`, and these satisfy
one 2nd-order linear ODE. So on each cell the parametric derivatives are just
`{Ψ, Ψ'}` (generic cell) — `C` has ≤ 2 elements — and (iii) reduces from the
general Riquier theorem to the **elementary existence theorem for a 2nd-order
ODE solved for `Ψ''`** (the recursion `aₖ ↦ aₖ₊₁` on Taylor coefficients, valid
where the initial `a₀+a₁v` is a unit — guaranteed by the cell's inequation).
Composition with `v = v₁x+v₂y+v₃z+v₄r` yields the solution as a function of the
spatial variables; the integrability `Ψ_xy = Ψ_yx` is automatic because `Ψ`
depends on `x,y,z` only through `v` and `v`'s mixed partials commute. In this
specialization the **δ-regularity / Riquier-ranking hypothesis disappears**
(one effective independent variable), so Lemma 3 holds for the hydrogen ansatz
with no genericity caveat. The general statement above retains it.

### Lemma 4 (splits and inequations are order-bounded)

*Every initial and separant on which `DTD(A)` splits is attached to a leader of
order `≤ N`, for `N ≥ N₀`.*

**Proof.** A split occurs when an initial or separant of a current equation is a
zero-divisor. New equations are produced only by prolongation and Δ-reduction,
all of order `≤ N` by Lemma 2; no equation, hence no leader, hence no
initial/separant carrying a split, has order `> N`. So the inequation set
`NE_i ∩ {order ≤ N}` of each cell already contains all of its splits. ∎

### Assembling

Fix `N ≥ N₀`. By Lemmas 1 and 2 the differential reductions performed by the
differential Thomas algorithm, restricted to order `≤ N`, are algebraic
reductions in `J_N` of elements of `(A^[N])`; by Lemma 4 every split is on an
order-`≤ N` initial/separant; so the *cell-splitting tree* of `DTD(A)` truncated
to order `≤ N` is exactly the splitting tree that `ATD(A^[N])` builds under the
matching ranking. The two decompositions therefore produce the same triangular
equation sets and the same inequation sets at order `≤ N`. Lemma 3 supplies the
two inclusions of solution sets (no spurious algebraic cell; every differential
cell truncates to an algebraic one), upgrading the syntactic match to a
bijection of cells whose solution sets agree. The higher-order part of each
differential cell is the (split-free) prolongation of its order-`≤ N` part,
again by Lemma 2 and `N ≥ N₀`. The corollary follows by running the same
argument on `A ∪ {P}` (`P` of order `≤ N` after the bound is enlarged to
`max(N₀, ord P)`). ∎

---

## 4. The bound `N₀` for the hydrogen ansatz

The ansatz (paper §"Completeness for this ansatz"):

- ODE element `(a₀+a₁v)Ψ'' + (b₀+b₁v)Ψ' + (c₀+c₁v)Ψ` — order **2** in `Ψ`
  along `v`;
- chain-rule relations `Ψ_x = Ψ'·v_x`, … — order **1**, expressing the
  `x,y,z`-derivatives through the single variable `v`;
- variable relation `v = v₁x+v₂y+v₃z+v₄r` — order 0 (algebraic);
- radius relation `r² − x²−y²−z²` — order 0 (algebraic), the only element
  nonlinear in its leader.

PDE `P`: `−½(Ψ_xx+Ψ_yy+Ψ_zz) − (1/r)Ψ = EΨ` — order **2**.

**Claim: `N₀ ≤ 2`, and `N := max(N₀, ord P) = 2` suffices.**

*Reasoning.* All differentiation runs through the single variable `v`: a mixed
spatial derivative reduces by the chain rule to `v`-derivatives of `Ψ`, e.g.
`Ψ_xy = Ψ''·v_x v_y + Ψ'·v_xy`. The only integrability obligations are the
commutations `Ψ_xy = Ψ_yx`, `Ψ_xz = Ψ_zx`, `Ψ_yz = Ψ_zy`, which hold **iff**
the mixed partials of `v` commute — and they do, since `v` is a fixed function
of `x,y,z` (with `r` algebraic, `r_xy = r_yx` from `r²=x²+y²+z²`). So every
critical pair closes at order 2, no integrability condition is born above order
2, and `R_2` is involutive. The PDE adds order-2 jets already present. Hence
prolonging `A∪{P}` to order 2 — a finite polynomial system in the order-`≤2`
jets of `Ψ` plus `x,y,z,r,v` and the 11 constants — carries the full content.

This is the same fact the paper invokes informally ("the second derivatives of
`v` are symmetric", §2379, lines 2422–2424); here it is the statement that the
involutive order is 2.

*Note (order vs degree).* `N₀ = 2` bounds the differentiation **order**. It says
nothing about polynomial **degree** — and the round-1 swell of the differential
port (task 444: reducing `r²⁰` modulo `r²−x²−y²−z²` makes `(x²+y²+z²)¹⁰`) is a
*degree* phenomenon in the algebraic reduction, not an order phenomenon. The
theorem moves that degree blow-up out of the differential engine and into a
**pure-algebra** computation, where a `degrevlex`/lazy ordering keeps
`r²−x²−y²−z²` as a triangular generator instead of substituting it to high
power. Controlling the degree is then a ranking/strategy choice handed to
Singular, not a swell baked into the reduction.

---

## 5. Gaps and what must be checked

1. **Exact `N₀`.** §4 argues `N₀ ≤ 2` from the single-`v` structure; this should
   be confirmed by an involutive (Janet) completion of the prolonged ansatz, and
   re-checked on the **degenerate cells** (where the ODE drops to first order on
   `a₀=a₁=0`): there the order can only drop, so `N=2` remains an upper bound,
   but the claim "`N` uniform over all cells" should be verified, since the
   theorem takes `N` = max over cells.

2. **Passive ⟺ involutive — RESOLVED (no longer borrowed).** The first draft
   imported this from Goldschmidt; the tightened Lemma 3 removes the import. A
   `DTD` cell is *by definition* (Bächler et al. 2012, Def. 3.5) algebraically
   simple **and** Janet-involutive, and Janet-involutivity is verbatim Robertz's
   passivity (LNM 2121, Def. 2.35). So passivity is a property of the output, not
   a theorem to be matched against Goldschmidt's symbol-regularity hypotheses.
   What remains genuinely external is only the **Riquier–Thomas existence
   theorem** (parametric assignment ⇒ formal solution) — a classical result, and
   for the hydrogen ansatz it degrades to the elementary 2nd-order-ODE existence
   theorem (Lemma 3, specialization), which carries no hypotheses to check.

2′. **Existence vs. bijection in the nonlinear case — RESOLVED.** Closed by
   **Lange-Hegermann 2014, Theorem 4.1**: for any algebraically restricted
   (hence nonlinear, inequation-bearing) differential system,
   `Sol_ζ(S)_{≤ℓ}` *equals* the truncated solution set of an algebraic Thomas
   decomposition in `C[G]_{≤ℓ}`, via the explicit bijection
   `ρ : {U}_Δ → G(U,Δ)`. This is an equality of solution sets, so it gives the
   cell bijection directly, no "existence-only" gap. The redundancy of
   non-reductive prolongations for an involutive system (so that `A^[N]` — the
   reductive prolongations to order `N` — suffices) is cited there to
   **Bächler et al. 2012, §3**. Net: the §2 theorem is essentially
   Lange-Hegermann's Theorem 4.1 specialized; only the explicit bound `N₀`
   (gap 1), the finiteness inheritance from Bächler's terminating differential
   Thomas decomposition, and the saturation-free solving corollary are added
   here.

3. **Comprehensive / parametric uniformity.** The cells split on the constants
   `c` too (the bad locus). The argument should be run with `c` as parameters
   (comprehensive ATD over `ℂ[c]`, à la Chen–Moreno Maza / the algebraic
   counterpart of parametric RG), and `N₀` shown independent of `c` — true if
   the system is of finite type uniformly in `c`, which holds here (the order
   structure is `c`-independent; only degrees/leaders drop on the bad locus,
   into separate cells).

4. **The two inclusions of Lemma 3 in the *parametric* setting.** Formal
   integrability is usually stated over a field; over `ℂ[c]` it must hold
   generically on each parameter cell, which is exactly what the
   cell-with-inequations structure provides. Worth making explicit.

5. **No spurious *components* vs no spurious *cells*.** Lemma 3 rules out
   truncation artifacts for an *involutive* `R_N`. If `N` is taken larger than
   `N₀` "to be safe", `A^[N]` may contain redundant generators but the zero set
   and hence `ATD` are unchanged; this is harmless but should be noted so the
   bound is not mistaken for a sharp invariant.

6. **Verification.** The cleanest external check: run an **algebraic** Thomas /
   comprehensive-triangular decomposition of `(A∪{P})^[2]` (Singular, or
   RegularChains as reference) and confirm it reproduces the 29 cells and their
   surviving solution varieties from the open-maple differential run
   (`~/Papers/NewMethod/joca-thomas-openmaple.out`). Agreement is direct
   evidence for the theorem and pins `N₀`; disagreement localizes the wrong
   lemma.

---

## 6. Related work and prior art

A literature check (June 2026) places this note against four established bodies
of work. The short version: **each of the three ingredients is published
separately, no single named method fuses all three, and the general form of the
one quantity we pin — the up-front prolongation bound `N₀` — is an explicitly
open problem.** So the contribution here is not the reduction idea but its
specialization to a finite-type ansatz where `N₀` becomes an explicit constant.

### 6.1 The reduction is classical — Rosenfeld's Lemma

"Differential solvability/membership collapses to *algebraic* solvability once a
coherent, regular system is reached" is **Rosenfeld's Lemma** (Rosenfeld 1959),
made algorithmic by **Boulier, Lazard, Ollivier & Petitot**, *Computing
representations for radicals of finitely generated differential ideals* (AAECC
20(1):73–121, 2009; algorithm orig. 1995) — the **Rosenfeld–Gröbner** algorithm.
There, Theorem 23 (Rosenfeld's Lemma) and Cor. 24 give: in a regular differential
system, a partially-reduced differential polynomial lies in `[A]:S^∞` iff it lies
in the *algebraic* ideal `(A):S^∞`; Theorem 25 lifts Lazard's Lemma to a
bijection of minimal differential primes with minimal algebraic primes. **What it
still does differentially:** it *computes the coherent/regular decomposition
itself* — autoreduction, Δ-polynomial/critical-pair coherence, prolongation as
needed — before the algebraic test applies. So Rosenfeld–Gröbner is the ancestor
of "defer to an algebraic algorithm," but it interleaves differential reasoning
throughout; it does not isolate the differential work into one bounded
prolongation.

### 6.2 The closest single result — Lange-Hegermann 2014

The precise construction of §2 — *treat every derivative of order `≤ N` as an
independent algebraic indeterminate, prolong only reductively to the involutive
order, then decompose algebraically* — is **essentially proved** in **Markus
Lange-Hegermann, "The Differential Counting Polynomial"** (arXiv:1407.5838, 2014;
and his 2014 RWTH-Aachen thesis). His §5.5 (proof of Thm 4.6) defines `S̄` by
applying the bijection
```
   ρ : {U}_Δ → G(U,Δ),   u_μ ↦ g_μ   (a fresh indeterminate per derivative),
```
extending the differential ranking to an algebraic ranking on `C[G]`, to the
equations **and only their reductive (bounded) prolongations**, and proves
```
   Sol_ζ(S)_{≤ℓ}  =  Sol_ζ(S̄_{≤ℓ})_{≤ℓ},
```
with the key line: *"the non-reductive prolongations are redundant, as `S` is
involutive."* That is our Lemmas 2–3, and our "derivative = algebraic
indeterminate + reductive prolongation only" (§2, §2ter step 2) verbatim.

**Sharp caveat (corrects this note's earlier looser citation).** Lange-Hegermann's
theorem **presupposes that `S` is already involutive** — a differential *simple
system*, the output of a differential Thomas decomposition. It therefore supplies
the algebraic-**finish** step *given involutivity*; it does **not** by itself
supply the a priori numeric bound `N₀` at which involutivity is *reached*. So we
may cite LH 2014 for the reduction (§2bis (a), Lemma 3) but **not** for §4's
`N₀ = 2`, which is the separate, added quantity. His Theorem 4.1 is also stated
per-order `ℓ` and allows a countable decomposition; our finiteness is inherited
from the terminating differential Thomas decomposition, not from his statement.

### 6.3 The algorithm-agnostic, strata-keeping finish — differential Thomas

The "keep every bad-locus stratum as its own cell, split on initials/separants/
discriminants" behaviour is a built-in guarantee of the **differential and
algebraic Thomas decomposition**: **Bächler, Gerdt, Lange-Hegermann & Robertz**,
*Algorithmic Thomas decomposition of algebraic and differential systems* (JSC
47(10):1233–1266, 2012; arXiv:1008.3767, 1108.0817). It returns a **disjoint**
partition into *simple* systems of equations **and inequations**; algebraic
simplicity = triangularity + squarefreeness + non-vanishing initials (so it
case-splits on initials, and via squarefreeness on **discriminants**); the
differential version additionally certifies **involutivity/passivity**. This is
the §2ter step-3/step-4 machinery, already published — our use of it is standard.
For the *algebraic* finish itself in the **parametric** case, the natural tool is
a **comprehensive triangular decomposition** (Chen, Golubitsky, Lemaire,
Moreno Maza & Pan, 2007) or comprehensive/algebraic Thomas.

**Why Thomas, specifically — it fuses both completion traditions.** The two
"prolong-until-nothing-new" completion tests that certify `N₀` come
from different lineages: **Janet / involutive-basis completion** (Janet–Riquier;
Gerdt–Blinkov involutive bases; Robertz LNM 2121) enumerates *non-multiplicative
prolongations* under an involutive division, reduces them *involutively* (unique
normal form), and yields the **Stanley/cone decomposition** — i.e. the
**parametric/principal derivative split** and the counting function; the
**Rosenfeld–Gröbner coherence test** (Ritt–Kolchin; Boulier–Lazard–Ollivier–
Petitot) instead forms *Δ-polynomials* for leader-overlapping pairs, reduces by
*differential pseudo-reduction modulo `:H_A^∞`*, and yields **regular differential
chains** (the membership reduction of Rosenfeld's lemma). By Seiler's unification
(*Involution*, 2010) these certify the *same* formal-integrability condition
(coherence ⟺ involution under a δ-regular ranking) — they are the combinatorial
and the membership-theoretic faces of one `N₀`. The **differential Thomas
decomposition is the hybrid that runs both at once**: Janet-division involutive
completion for the involutivity certificate, plus RG-style initial/separant
*inequation-splitting* for the disjoint, saturation-free cells. That is exactly
why this note leans on it — Lemma 3(ii)'s parametric/principal split is the Janet
side, and the "differential = algebraic at bounded order" membership claim is the
coherence side, and Thomas supplies both in one disjoint decomposition.

### 6.4 The bound `N₀` in general is an open problem — Golubitsky et al. 2008

The theory that *guarantees* a finite involutive order is **Cartan–Kuranishi /
Seiler** (W. Seiler, *Involution: The Formal Theory of Differential Equations*,
Springer 2010): any regular system reaches involution after finitely many
prolongations and projections, and "the prolongation of an involutive symbol is
again involutive" (whence *no new integrability condition is born above `N₀`* —
exactly our stopping property). **But the Cartan–Kuranishi geometric statement is
a termination/stability result, not an effective numeric bound.**

*(Open-access: the book is expensive, but its algebraic/combinatorial core — the
part this note uses — is Seiler's two arXiv papers "A Combinatorial Approach to
Involution and δ-Regularity" I (arXiv:math/0208247, AAECC 20:207–259) and II
(arXiv:math/0208250, AAECC 20:261–338), both 2009; local copies in
`~/project/papers/`. These cover involutive/Pommaret bases, δ-regularity, and the
completion mechanics. Part II carries the one result closest to an **effective**
involutive bound: the **degree of a Pommaret basis equals the Castelnuovo–Mumford
regularity** of the ideal (in δ-regular coordinates) — i.e. the algebraic
completion degree is `reg(I)`, a computable invariant. That bounds the
completion order for the **algebraic** (Pommaret/involutive-basis) side; it does
not by itself pin the **differential** prolongation order `N₀`, but it is the
right handle if we ever want an unconditional bound for the algebraized system
`A^[N]` rather than the structural `N₀ = 2` argument of §4.)*

A *computable* up-front bound is known only in the **ordinary** (single-derivation)
case: **Golubitsky, Kondratieva, Moreno Maza & Ovchinnikov, "A bound for the
Rosenfeld-Gröbner algorithm"** (JSC 43(8):582–610, 2008; arXiv:math/0702470)
prove `M(F) ≤ (n−1)!·M(F₀)` on the summed derivative orders of a modified
Rosenfeld–Gröbner run whose algorithm "performs all differentiations in the
beginning and then uses a purely algebraic decomposition algorithm." In the same
paper they state the general question as **open** and conjecture its difficulty:

> "It would be good to have a bound that would tell us how many times we need to
> differentiate the original system in the beginning… so that the rest of the
> computation can be performed by a purely algebraic decomposition algorithm… we
> do not provide such a bound and, moreover, **conjecture that it would have
> solved the Ritt problem**."

**Calibration (do not overstate).** This is a *conjectural, one-directional*
link — (a general computable up-front differentiation bound) ⇒ (the Ritt problem
is solvable) — for the **Rosenfeld–Gröbner** differentiation count, **not** a
proven *equivalence*, and **not** specifically about the Janet-involutive
completion order. The companion paper **Golubitsky–Kondratieva–Ovchinnikov,
*On the generalised Ritt problem as a computational problem*** (J. Math. Sci.
163(5):515–522, 2009; arXiv:0809.1128 — local copy in
`~/project/papers/`) lists **six** equivalent formulations of the (generalised)
Ritt problem (Thm 1: char-set→generators; inclusion of two primes; non-redundant
prime decomposition; primality of a radical ideal by generators; prime
decomposition with generators; and `f` a zero-divisor mod `I`). **"Compute the
involutive / prolongation bound" is NOT among them, and no cited result makes it
an equivalent formulation.** The reason is instructive and *strengthens* the
positioning below rather than weakening it.

**Why the bound is not a reformulation but the crux.** The construction of §2/§2ter
is *literally* the one inside GKO-2009's proofs of (4⇒1) and (6⇒1): from a
characteristic set `C` they form the **algebraic** ideals
`Jᵢ = (C^(i)) : H_C^∞`, where `C^(i) = {θf : ord θ ≤ i, f ∈ C}` is the
**prolongation to order `i`** with derivatives read as ordinary indeterminates,
take a Gröbner basis `Fᵢ`, and note *"by the basis theorem, there exists an index
`i` such that `{Fᵢ} = P`."* That existential `i` **is** our `N₀`. Its existence is
free — but by the **Ritt–Raudenbush basis theorem**, not by Hilbert: the `{Fᵢ}`
are radical *differential* ideals forming an ascending chain
`{F₀} ⊆ {F₁} ⊆ ⋯ ⊆ P`, and Ritt–Raudenbush is exactly the ACC on *radical*
differential ideals (the differential analogue of Hilbert's theorem — and it
holds *only* for radical/perfect differential ideals; a differential polynomial
ring is not Noetherian as a differential ring). GKO name it explicitly on p. 7
("by the Ritt–Raudenbush theorem, a strictly increasing chain of radical
differential ideals terminates"); the (4⇒1) proof on p. 5 shorthands it as "the
basis theorem." Our own **prolong-until-stable** (§2bis) needs only the *weaker*
ordinary ACC: it watches the elimination ideal `J_c ∩ K[c]` inside the fixed,
finitely-generated commutative ring `K[c]`, whose chain stabilizes by Hilbert's
basis theorem alone — no Ritt–Raudenbush required. Either way, what no one can do
in general is compute that `i` *a priori* without a primality/zero-divisor oracle
— which is why the bound is not a separate equivalent problem but the
**effectivity gap** the whole Ritt problem turns on. So the earlier `N₀`-vs-Ritt link is at most: *an a-priori bound would
remove the oracle from the (4⇒1)/(6⇒1) search* — the conjectural one-directional
implication of Golubitsky et al. 2008, nothing stronger.

**Two distinctions this sharpens (both citable against GKO-2009).**
1. **Saturation.** GKO form `(C^(i)) : H_C^∞` — the `: H_C^∞` divides out the
   initials/separants, *deleting the bad locus*, because they chase a single prime
   `P`. This note is deliberately **saturation-free** (§2bis (c)): the bad-locus
   strata are kept as their own cells. We keep precisely what their construction
   discards.
2. **Effectivity.** GKO locate `i` by an increasing search with a primality
   oracle; §4 computes `i = N₀ = 2` **a priori** from the single-`v` structure —
   an effective bound for a finite-type class where the general effective bound is
   open. That is the one place the note answers, for this ansatz family, the "there
   exists an `i`" that their proof leaves non-effective.

**Honesty note.** The prolong-to-order-`i` / algebraic-ideal / Gröbner
construction is therefore **not novel** — it is textbook by GKO-2009 (and older).
The note should *cite* it and claim only (1) the saturation-free bad-locus reading
and (2) the explicit `N₀ = 2`, not the pipeline itself.

### 6.5 Parametric / comprehensive analogues

The parametric case-split on parameters (our bad-locus `c`-splitting) mirrors
**comprehensive Gröbner systems** (Weispfenning 1992; Suzuki–Sato 2006;
Kapur–Sun–Wang, ISSAC 2010; Montes' Gröbner Cover), and on the involutive side
**Gerdt & Hashemi, "Comprehensive Involutive Systems"** (CASC 2012,
arXiv:1206.0181) — a Janet-basis analogue that partitions parameter space into
cells each carrying a uniform involutive basis. A **parametric Rosenfeld–Gröbner**
(all regular representations across all parameter states) is given by Hashemi et
al. (*Cogent Mathematics*, 2018). These are the closest parametric relatives; the
verified sources do **not**, however, show any of them fused with a *single
bounded involutive prolongation* and an algorithm-agnostic algebraic finish.

### 6.6 Net positioning of this note

- **Not new:** the derivative↔indeterminate reduction and "reductive prolongations
  suffice for an involutive system" (Lange-Hegermann 2014); the
  prolong-to-order-`i` → algebraic-ideal `(C^(i)):H_C^∞` → Gröbner →
  increasing-search construction and its Noetherian termination
  (Golubitsky–Kondratieva–Ovchinnikov 2009, proofs of Thm 1 (4⇒1)/(6⇒1)); the
  disjoint, involutive, strata-keeping algebraic finish (Thomas decomposition
  2012); comprehensiveness over parameters (Weispfenning; Gerdt–Hashemi 2012).
- **The contribution:** (i) the **explicit `N₀ = 2`** for the finite-type
  ODE-tower ansatz — an *effective* value where the general bound is open
  (Golubitsky et al. 2008) — proved from the single-`v` structure (§4);
  (ii) the **saturation-free / bad-locus** reading, i.e. running the reduction on
  `A ∪ {P}` recovers exactly the strata Theorem 1's projection deletes on
  `h = 0` (§2bis (c), Corollary); (iii) the explicit **algorithm-agnostic**
  packaging (§2bis (d), §2ter) that lets GTZ / regular chains / Gröbner /
  algebraic Thomas each finish the *same* fixed ideal.

Honest framing for the paper: **specialize a known equivalence (LH 2014) to a
finite-type ansatz where the otherwise-open prolongation bound (Golubitsky et al.
2008) becomes an explicit `N₀ = 2`, and read the un-saturated form for the bad
locus.** That is well-cited and still a genuine contribution; it is *not* a new
general reduction theorem.

---

## 7. Relationship to Theorem 1 of the paper

This is not prior art — it is the paper's own companion result, and it can look
as if the two contradict: Theorem 1 restricts its guarantee to `h(c*) ≠ 0`, while
the saturation-free reading claims to find *all* solutions (bad locus included),
possibly with some spurious. **There is no contradiction.** They are guarantees
about *two different algorithms*, and they agree on every axis that matters.

**What Theorem 1 (`thm:completeness`) actually claims — and disclaims.** Under
hypotheses (1)–(4), if `c*` solves the PDE *and* `h(c*) ≠ 0`, then `c* ∈ 𝒱`. The
paper states both directions explicitly:
- **completeness** on `{h ≠ 0}` — "every valid choice of constants will be found";
- **not sufficiency** — "there may exist points in `𝒱` that do not actually solve
  the PDE, for example because the denominator vanishes there";
- and it may **miss the bad locus** — where `h = 0`, "the reduction degenerates
  and the algorithm may miss solutions."

So Theorem 1's output `𝒱` is *already* a **superset of the generic (`h ≠ 0`)
solutions that may contain spurious points** — the "complete but not sufficient"
character is Theorem 1's, not something the bounded-prolongation reading
introduces.

**Same character, saturation removed.** Line up the two algorithms:

| | Theorem 1 (generic, **saturating**) | bounded-prolongation (**saturation-free**) |
|---|---|---|
| completeness | all solutions with `h ≠ 0`; **may miss** `h = 0` | **all** solutions, `h = 0` included |
| sufficiency | **not** guaranteed — `𝒱` may hold spurious | **not** guaranteed — closure may hold spurious |
| mechanism | invert `H_A` (divide by `h`), then project | prolong to `N₀`, **don't** divide, then project |

On **sufficiency** the two *agree* — both are supersets that may carry spurious
constants (for the saturation-free method these are the closure-minus-image locus
of §2ter′'s superset remark). On **completeness** the saturation-free method is
strictly *stronger* — it recovers the bad-locus solutions Theorem 1 may drop.
Stronger on completeness, equal on sufficiency, is an **improvement, not a
contradiction**; and it cannot contradict Theorem 1, which never claimed the bad
locus (hypothesis (2) fences it off).

**The single point of divergence is the saturation.** Theorem 1's proof passes to
`sat(A(c*))` via Rosenfeld's Lemma, which is exactly what forces `h(c*) ≠ 0` and
exactly what deletes the components inside `{h = 0}`. The bounded-prolongation
reading replaces "invert `H_A`" with "prolong to `N₀`," never saturates, and so
keeps those components. **Removing the saturation is precisely what lifts
hypothesis (2).**

Net: Theorem 1's `𝒱` is a superset of the *generic* solutions that may contain
spurious points; the bounded-prolongation reading is the *same* complete-but-not-
sufficient superset with the saturation removed — extending completeness to *all*
constants at no cost on the sufficiency side.
