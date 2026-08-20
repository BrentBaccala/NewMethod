/-
# The two algorithms of *NewMethod* in Lean 4 / Mathlib

This file states, in Lean 4 against Mathlib, the two algorithms of

  B. Baccala, *A New Method for Solving Systems of Differential Equations
  with an Ansatz*, `~/Papers/NewMethod/NewMethod.tex`

together with every supporting notion their Input/Output clauses mention.

* **Algorithm 1**, `ConsistencyLocus` — `NewMethod.tex:1594` (`alg:consistency`),
  in `\subsection{Computing the Consistency Locus}` (`sec:consistency-locus`,
  `NewMethod.tex:1530`).
* **Algorithm 2**, `MembershipLocus` — `NewMethod.tex:1657` (`alg:core`),
  in `\subsection{Computing the Membership Locus}` (`sec:membership-locus`,
  `NewMethod.tex:1649`).

## What this file is for

Fidelity, not proof-completion.  A reader should be able to put this file
beside the LaTeX and check line by line that the Lean says what the paper
says.  Every `axiom` names its source, the source's own numbering, and the
local file the statement was read from.  Departures from a source's own
phrasing are flagged `NOTE (departure from source):`.  Suspected errors in
the paper are flagged `NOTE (suspected error in source):` and are *not*
silently repaired — the Lean follows the paper.

## Provenance of reused material

§9 (`DifferentialRing`, `TriangularSet`, `rittReduce`, `rosenfeld`,
`lazard`, `lazard_gap`, `ritt_remainder_zero`) comes from
`~/project/Completeness.lean`, the earlier Lean formalization of this same
paper's projection-completeness theorem (`sec:completeness`,
`NewMethod.tex:1911`).  Names are preserved so the two developments can be
read together.  `diffIdeal`, a `sorry` stub there, is given a real
definition here (`diffIdealBy`, §1).

## Contents

* §1  differential ideals, saturation
* §2  the input data of both algorithms (the Input clauses)
* §3  systems, solutions, simple systems (ThomasDecomp §§2–3)
* §4  constructible sets, `V(·)`, `Z(·)` (CTD Def. 6; NewMethod.tex:1568-1592)
* §5  cells (NewMethod.tex:1550-1567)
* §6  the four subroutines
* §7  the three constant loci (NewMethod.tex:1327-1500)
* §8  the two algorithms and their output correctness
* §9  bridge to `~/project/Completeness.lean`

## Checking

    ~/axiommath.ai/check.py ~/Papers/NewMethod/Lean/NewMethodAlgorithms.lean

See `README.md` in this directory for the last-run output and the axiom /
`sorry` inventories.
-/

import Mathlib

noncomputable section

namespace NewMethod

/-! ############################################################
# §1.  Differential ideals and saturation

Sources
* `NewMethod.tex:1082-1089` (`sec:preliminaries`) — differential ideal
  `[g₁,…,g_n]`, radical `√I`, saturation `I : S^∞`, consistency `1 ∉ I`.
* Kolchin, *Differential Algebra and Algebraic Groups*, §I.1–I.2
  (`~baccala/Books/Differential Algebra/Kolchin - DAAG.pdf`).
############################################################ -/

section DiffIdeal

variable {A : Type} [CommRing A] {m : ℕ}

/-- `d` is a family of `m` pairwise commuting derivations.

Source: `NewMethod.tex:1082`, and the Input clause of both algorithms
(`NewMethod.tex:1596-1600`, `NewMethod.tex:1659-1663`), which fixes
`Δ = {δ₁,…,δ_m}`. -/
structure IsDerivations (d : Fin m → A → A) : Prop where
  add  : ∀ i x y, d i (x + y) = d i x + d i y
  leib : ∀ i x y, d i (x * y) = d i x * y + x * d i y
  comm : ∀ i j x, d i (d j x) = d j (d i x)

/-- An ideal closed under every derivation — a *differential ideal*.

Source: `NewMethod.tex:1082-1084`: "A differential ideal `I` is a subring of
a differential ring `R`, closed under arbitrary derivation and
multiplication by arbitrary ring elements."

NOTE (suspected error in source): the paper says *subring*; an ideal is not
a subring (it need not contain `1` — indeed a *consistent* one must not).
We take the standard reading — an ideal closed under the derivations —
which is what every later use in the paper requires. -/
def IsDiffIdeal (d : Fin m → A → A) (J : Ideal A) : Prop :=
  ∀ i x, x ∈ J → d i x ∈ J

/-- `[S]` — the smallest differential ideal containing `S`.

Source: `NewMethod.tex:1084-1085`: "Given a set of generators `g₁,…,g_n`,
the differential ideal they generate is denoted `[g₁,…,g_n]`."

NOTE: in `~/project/Completeness.lean` this was the `sorry` stub
`def diffIdeal (S : Set R) : Ideal R := sorry`.  Here it is a real
definition, the infimum of the differential ideals containing `S`. -/
def diffIdealBy (d : Fin m → A → A) (S : Set A) : Ideal A :=
  sInf {J : Ideal A | S ⊆ (J : Set A) ∧ IsDiffIdeal d J}

lemma diffIdealBy_le {d : Fin m → A → A} {S : Set A} {J : Ideal A}
    (hS : S ⊆ (J : Set A)) (hd : IsDiffIdeal d J) : diffIdealBy d S ≤ J :=
  sInf_le ⟨hS, hd⟩

lemma subset_diffIdealBy {d : Fin m → A → A} {S : Set A} :
    S ⊆ (diffIdealBy d S : Set A) := by
  intro y hy
  simp only [diffIdealBy, SetLike.mem_coe, Submodule.mem_sInf, Set.mem_setOf_eq]
  rintro J ⟨hJ, -⟩
  exact hJ hy

lemma diffIdealBy_isDiffIdeal {d : Fin m → A → A} {S : Set A} :
    IsDiffIdeal d (diffIdealBy d S) := by
  intro i y hy
  simp only [diffIdealBy, Submodule.mem_sInf, Set.mem_setOf_eq] at hy ⊢
  intro J hJ
  exact hJ.2 i y (hy J hJ)

lemma diffIdealBy_mono {d : Fin m → A → A} {S T : Set A} (h : S ⊆ T) :
    diffIdealBy d S ≤ diffIdealBy d T :=
  diffIdealBy_le (h.trans subset_diffIdealBy) diffIdealBy_isDiffIdeal

/-- `I : q^∞` — saturation of the ideal `I` by the powers of `q`.

Source: `NewMethod.tex:1087`.

NOTE (suspected error in source): the paper writes
`I : S^∞ = {p : ∃ q ∈ R, q p ∈ I}`.  As written this is wrong twice over:
the multiplier `q` ranges over all of `R` rather than over the powers of
`S`, and `S` does not occur on the right-hand side at all.  Taken
literally, every `p` qualifies (choose `q = 0`), so `I : S^∞` would be the
whole ring.  Every *use* of the notation in the paper (Prop.
`simple-membership`, the `V_∃`/`V_∀` definitions, Thm. `completeness`) is
the standard saturation, so that is what is formalized here. -/
def satBy (I : Ideal A) (q : A) : Ideal A where
  carrier := {y | ∃ k : ℕ, q ^ k * y ∈ I}
  zero_mem' := ⟨0, by simp⟩
  add_mem' := by
    rintro a b ⟨k, hk⟩ ⟨l, hl⟩
    refine ⟨k + l, ?_⟩
    have h : q ^ (k + l) * (a + b) = q ^ l * (q ^ k * a) + q ^ k * (q ^ l * b) := by
      ring
    rw [h]
    exact I.add_mem (I.mul_mem_left _ hk) (I.mul_mem_left _ hl)
  smul_mem' := by
    rintro t a ⟨k, hk⟩
    refine ⟨k, ?_⟩
    have h : q ^ k * (t • a) = t * (q ^ k * a) := by
      simp only [smul_eq_mul]; ring
    rw [h]
    exact I.mul_mem_left _ hk

lemma le_satBy {I : Ideal A} {q : A} : I ≤ satBy I q :=
  fun _ hx => ⟨0, by simpa using hx⟩

lemma satBy_mono {I J : Ideal A} {q : A} (h : I ≤ J) : satBy I q ≤ satBy J q := by
  rintro y ⟨k, hk⟩; exact ⟨k, h hk⟩

lemma satBy_satBy {I : Ideal A} {q : A} : satBy (satBy I q) q = satBy I q := by
  refine le_antisymm ?_ le_satBy
  rintro y ⟨k, l, hl⟩
  refine ⟨l + k, ?_⟩
  have h : q ^ (l + k) * y = q ^ l * (q ^ k * y) := by ring
  rw [h]; exact hl

/-- The radical of a saturation is itself saturated. -/
lemma satBy_radical_satBy {I : Ideal A} {q : A} :
    satBy (satBy I q).radical q = (satBy I q).radical := by
  refine le_antisymm ?_ le_satBy
  rintro y ⟨k, hk⟩
  obtain ⟨p, hp⟩ := hk
  obtain ⟨l, hl⟩ := hp
  have hy : y ^ p ∈ satBy I q := by
    refine ⟨l + k * p, ?_⟩
    have h : q ^ (l + k * p) * y ^ p = q ^ l * (q ^ k * y) ^ p := by ring
    rw [h]; exact hl
  exact ⟨p, hy⟩

/-- Leibniz for a power times an element:
`q · δ(qᵏ x) = q^{k+1} δx + k (δq)(qᵏ x)`. -/
lemma deriv_pow_mul {d : Fin m → A → A} (hd : IsDerivations d)
    (i : Fin m) (q y : A) (k : ℕ) :
    q * d i (q ^ k * y) = q ^ (k + 1) * d i y + (k : A) * d i q * (q ^ k * y) := by
  induction k with
  | zero => simp
  | succ k ih =>
      have hq : q ^ (k + 1) * y = q * (q ^ k * y) := by ring
      have h1 : d i (q ^ (k + 1) * y)
          = d i q * (q ^ k * y) + q * d i (q ^ k * y) := by
        rw [hq, hd.leib]
      rw [h1, mul_add, ← mul_assoc, mul_comm q (d i q)]
      rw [show q * (q * d i (q ^ k * y)) = q * (q * d i (q ^ k * y)) from rfl]
      have h3 : q * (q * d i (q ^ k * y))
          = q * (q ^ (k + 1) * d i y + (k : A) * d i q * (q ^ k * y)) := by
        rw [ih]
      rw [h3]
      push_cast
      ring

/-- Saturating a differential ideal by a single element keeps it
differential. -/
lemma satBy_isDiffIdeal {d : Fin m → A → A} (hd : IsDerivations d)
    {I : Ideal A} (hI : IsDiffIdeal d I) (q : A) :
    IsDiffIdeal d (satBy I q) := by
  rintro i y ⟨k, hk⟩
  refine ⟨k + 1, ?_⟩
  have hkey := deriv_pow_mul hd i q y k
  have h : q ^ (k + 1) * d i y
      = q * d i (q ^ k * y) - (k : A) * d i q * (q ^ k * y) := by
    rw [hkey]; ring
  rw [h]
  exact I.sub_mem (I.mul_mem_left _ (hI i _ hk)) (I.mul_mem_left _ hk)

/-- **Radical of a differential ideal is a differential ideal.**

Source: Kolchin, *Differential Algebra and Algebraic Groups*, §I.2,
Lemma 1(a) (`~baccala/Books/Differential Algebra/Kolchin - DAAG.pdf`);
also Ritt 1950, Ch. I (`~baccala/Books/Differential Algebra/Ritt J.F.
Differential algebra (AMS, 1950)(T)(189s).djvu`).  Requires characteristic
zero (or a Ritt algebra), which holds here: `R` is a `ℚ`-algebra.

SOURCE CONSULTED: Ritt Ch. I read via `djvutxt`; the statement is the
standard one and is the fact the paper appeals to at `NewMethod.tex:1387`
("`√([A(c*)] : Q^∞)` is a differential ideal"). -/
axiom radical_isDiffIdeal {d : Fin m → A → A} (hd : IsDerivations d)
    {I : Ideal A} (hI : IsDiffIdeal d I) : IsDiffIdeal d I.radical

/-- Consistency: `1 ∉ I`.  Source: `NewMethod.tex:1089`. -/
def Consistent (I : Ideal A) : Prop := (1 : A) ∉ I

lemma consistent_of_le {I J : Ideal A} (h : I ≤ J) (hJ : Consistent J) :
    Consistent I := fun h1 => hJ (h h1)

lemma consistent_radical {I : Ideal A} (h : Consistent I.radical) : Consistent I :=
  consistent_of_le Ideal.le_radical h

end DiffIdeal

/-! ############################################################
# §2.  The input data of both algorithms

The Input clauses of Algorithm 1 (`NewMethod.tex:1596-1613`) and
Algorithm 2 (`NewMethod.tex:1659-1679`) are identical apart from
Algorithm 2's extra flag `μ`.  They are bundled here as one structure.

Type parameters
* `R`  — the differential polynomial ring
  `R = ℚ[x₁,…,x_m]{u₁,…,u_k, c₁,…,c_n}`.
* `Rc` — the ring after specializing the constants at a point of `ℂⁿ`,
  i.e. `ℂ[x₁,…,x_m]{u₁,…,u_k}`.  Needed only to state the loci (§7).
* `F`  — the ring in which a solution takes its values (a ring of
  functions).  Needed only to state `Sol` (§3).
* `Ω`  — the space of solutions.
############################################################ -/

/-- The ring of constants `ℚ[c₁,…,c_n]`.  Source: `NewMethod.tex:1571`. -/
abbrev Cst (n : ℕ) := MvPolynomial (Fin n) ℚ

/-- Constant space `ℂⁿ`.  Source: `NewMethod.tex:1336`. -/
abbrev Pt (n : ℕ) := Fin n → ℂ

/-- Evaluation at a point of `ℂⁿ`, as a ring homomorphism
`ℚ[c₁,…,c_n] → ℂ`. -/
def evalHom {n : ℕ} (c : Pt n) : Cst n →+* ℂ :=
  ((MvPolynomial.aeval (R := ℚ) c : Cst n →ₐ[ℚ] ℂ) : Cst n →+* ℂ)

/-- Evaluation of a constant polynomial at a point of `ℂⁿ`. -/
def evalPt {n : ℕ} (c : Pt n) (p : Cst n) : ℂ := evalHom c p

@[simp] lemma evalPt_one {n : ℕ} (c : Pt n) : evalPt c 1 = 1 :=
  map_one (evalHom c)

@[simp] lemma evalPt_mul {n : ℕ} (c : Pt n) (p q : Cst n) :
    evalPt c (p * q) = evalPt c p * evalPt c q :=
  map_mul (evalHom c) p q

@[simp] lemma evalPt_add {n : ℕ} (c : Pt n) (p q : Cst n) :
    evalPt c (p + q) = evalPt c p + evalPt c q :=
  map_add (evalHom c) p q

@[simp] lemma evalPt_zero {n : ℕ} (c : Pt n) : evalPt c 0 = 0 :=
  map_zero (evalHom c)

/-- A *system* `(S^=, S^≠)`: equations and inequations.

Source: `NewMethod.tex:1211-1213` (`sec:preliminaries`) and
ThomasDecomp §2.1 (`~/project/papers/bachler-gerdt-lange-hegermann-robertz-
2012-differential-thomas-decomposition.pdf`, p. 3): "A finite set of
equations and inequations is called an (algebraic) system over `R`." -/
structure DiffSystem (R : Type) where
  eqs   : List R
  ineqs : List R

/-- The Input clause of Algorithms 1 and 2.

`NewMethod.tex:1596-1613` / `NewMethod.tex:1659-1679`. -/
structure AlgorithmInput (R Rc F Ω : Type) [CommRing R] [CommRing Rc]
    [Algebra ℂ Rc] [CommRing F] [Algebra ℂ F] (m k n : ℕ) where
  /-- The derivations `Δ = {δ₁,…,δ_m}`. -/
  d : Fin m → R → R
  d_isDeriv : IsDerivations d
  /-- The independent variables `x₁,…,x_m`. -/
  x : Fin m → R
  /-- The differential indeterminates `u₁,…,u_k`. -/
  u : Fin k → R
  /-- The constant parameters `c₁,…,c_n`. -/
  c : Fin n → R
  /-- `δᵢ xⱼ = δᵢⱼ` (`NewMethod.tex:1599`). -/
  d_x : ∀ i j, d i (x j) = if i = j then 1 else 0
  /-- `δᵢ c_l = 0` (`NewMethod.tex:1599`). -/
  d_c : ∀ i l, d i (c l) = 0
  /-- The inclusion `ℚ[c₁,…,c_n] ↪ R`. -/
  cst : Cst n →+* R
  cst_X : ∀ l, cst (MvPolynomial.X l) = c l
  cst_inj : Function.Injective cst
  /-- Decides `r ∈ ℚ[c₁,…,c_n]`, returning the preimage.  Used for the
  operation `S ∩ ℚ[c₁,…,c_n]` of lines 3–4 / 8–9. -/
  toCst : R → Option (Cst n)
  toCst_spec : ∀ r p, toCst r = some p ↔ cst p = r
  /-- The ranking `≺` (`NewMethod.tex:1097-1108`). -/
  rank : R → R → Prop
  /-- `u < v → δu < δv` (`NewMethod.tex:1102`). -/
  rank_deriv_mono : ∀ i v w, rank v w → rank (d i v) (d i w)
  /-- `u < δu` (`NewMethod.tex:1103`). -/
  rank_lt_deriv : ∀ i v, rank v (d i v)
  /-- Block ranking `{u} ≫ {c}` (`NewMethod.tex:1610`). -/
  rank_c_lt_u : ∀ (l : Fin n) (j : Fin k), rank (c l) (u j)
  /-- Block ranking `{c} ≫ {x}` (`NewMethod.tex:1610`). -/
  rank_x_lt_c : ∀ (i : Fin m) (l : Fin n), rank (x i) (c l)
  /-- Leader: the highest-ranking derivative occurring in `p`
  (`NewMethod.tex:1120`). -/
  leader : R → R
  /-- Initial: leading coefficient of `p` as a univariate polynomial in its
  leader (`NewMethod.tex:1122`). -/
  initial : R → R
  /-- Separant `s_p = ∂p/∂u_p` (`NewMethod.tex:1123`). -/
  separant : R → R
  /-- `Coeffs(r, ℚ[c₁,…,c_n])` — the coefficients of the unique expansion
  `r = Σ_α m_α p_α` into distinct power products `m_α` in the non-constant
  indeterminates times coefficients `p_α ∈ ℚ[c]`
  (`NewMethod.tex:1786-1800`, eq. `coeff-operator`). -/
  coeffs : R → List (Cst n)
  /-- `ev ω` evaluates a differential polynomial at a solution `ω`. -/
  ev : Ω → R →+* F
  /-- The projection `Ω → ℂⁿ` reading off the values of the constants. -/
  proj : Ω → Pt n
  /-- The projection is compatible with evaluation on `ℚ[c]`. -/
  ev_cst : ∀ ω p, ev ω (cst p) = algebraMap ℂ F (evalPt (proj ω) p)
  /-- Specialization `c ↦ c*`: `R = ℚ[x]{u,c} → ℂ[x]{u} = Rc`. -/
  spec : Pt n → R →+* Rc
  spec_cst : ∀ (cs : Pt n) (p : Cst n),
    spec cs (cst p) = (algebraMap ℂ Rc) (evalPt cs p)
  /-- Derivations on the specialized ring. -/
  dc : Fin m → Rc → Rc
  dc_isDeriv : IsDerivations dc
  /-- The ansatz `A(c)` (`NewMethod.tex:1607`). -/
  ansatz : List R
  /-- The target equations `𝒫^= = {P₁,…,P_t}` (`NewMethod.tex:1602`). -/
  Peq : List R
  /-- The target inequations `𝒫^≠ = {Q₁,…,Q_w}` (`NewMethod.tex:1603`). -/
  Pne : List R

namespace AlgorithmInput

variable {R Rc F Ω : Type} [CommRing R] [CommRing Rc] [Algebra ℂ Rc]
variable [CommRing F] [Algebra ℂ F]
variable {m k n : ℕ}

/-- `Q = Q₁ ⋯ Q_w`, with `Q = 1` when `𝒫^≠` is empty
(`NewMethod.tex:1343-1345`). -/
def Q (I : AlgorithmInput R Rc F Ω m k n) : R := I.Pne.prod

/-- `S ∩ ℚ[c₁,…,c_n]` — the constant-only part of a list of differential
polynomials.  Lines 3–4 of Algorithm 1 (`NewMethod.tex:1618-1619`) and
lines 8–9 of Algorithm 2 (`NewMethod.tex:1690-1691`). -/
def constPart (I : AlgorithmInput R Rc F Ω m k n) (L : List R) : List (Cst n) :=
  L.filterMap I.toCst

/-- The set of initials and separants `H_S = I_S ∪ S_S`
(`NewMethod.tex:1124-1126`). -/
def initSep (I : AlgorithmInput R Rc F Ω m k n) (S : List R) : List R :=
  S.map I.initial ++ S.map I.separant

/-- `q` — the product of the initials and separants of the equations of a
system.  Source: `NewMethod.tex:1281-1282` (the `q` of
Prop. `simple-membership`) and `NewMethod.tex:1978`. -/
def qOf (I : AlgorithmInput R Rc F Ω m k n) (S : List R) : R :=
  (I.initSep S).prod

end AlgorithmInput

/-! ############################################################
# §3.  Systems, solutions, and simple systems

Source: Bächler–Gerdt–Lange-Hegermann–Robertz, *Algorithmic Thomas
Decomposition of Algebraic and Differential Systems*, JSC 47 (2012),
§§2–3, read from
`~/project/papers/bachler-gerdt-lange-hegermann-robertz-2012-differential-
thomas-decomposition.pdf`; restated by the paper at
`NewMethod.tex:1230-1266`.
############################################################ -/

namespace AlgorithmInput

variable {R Rc F Ω : Type} [CommRing R] [CommRing Rc] [Algebra ℂ Rc]
variable [CommRing F] [Algebra ℂ F]
variable {m k n : ℕ}

/-- `Sol(S^=, S^≠)` — the solutions of a system: the equations vanish, the
inequations do not.

Source: ThomasDecomp §2.1, p. 3: "We call `a` a solution of a system `S`
if it is a solution of each element in `S`", where a solution of `p^=` is
`φ_a(p) = 0` and of `p^≠` is `φ_a(p) ≠ 0`.  Restated at
`NewMethod.tex:1211`. -/
def Sol (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : Set Ω :=
  {ω | (∀ p ∈ S.eqs, I.ev ω p = 0) ∧ (∀ q ∈ S.ineqs, I.ev ω q ≠ 0)}

/-- `Sol(S_{<x})` — the solutions of the subsystem of the elements whose
leader ranks strictly below `x`.

Source: ThomasDecomp §2.1, p. 3: "`S_{<x} := {p ∈ S | ld(p) < x}` is a
system over `F[y | y < x]`."

NOTE (departure from source): the source builds the subsystem `S_{<x}` and
then takes its solutions; here the two steps are fused, because building
`S_{<x}` as a `List` would need `I.rank` to be decidable and nothing later
uses the subsystem itself. -/
def SolBelow (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) (v : R) :
    Set Ω :=
  {ω | (∀ p ∈ S.eqs, I.rank (I.leader p) v → I.ev ω p = 0) ∧
       (∀ q ∈ S.ineqs, I.rank (I.leader q) v → I.ev ω q ≠ 0)}

/-- (A3) square-freeness, as a per-equation predicate.

Source: ThomasDecomp Definition 2.2(3): "`S` is square-free if the
univariate polynomial `φ_{<x_i,a}(p) ∈ F[x_i]` is square-free for all
`a ∈ Sol(S_{<x_i})` and `p ∈ S_{x_i}`."  Restated at
`NewMethod.tex:1244-1247`.

NOTE (definitional stub): stating this needs the partial-evaluation
homomorphism `φ_{<x,a}` and a squarefreeness predicate on univariate
polynomials over `F̄`; Mathlib has `Squarefree` but no `φ_{<x,a}`, and
building the jet-space substitution apparatus is out of scope.  The
predicate is therefore a field of the input structure below, and this is a
declared debt. -/
def SquareFreeAt (_I : AlgorithmInput R Rc F Ω m k n)
    (_S : DiffSystem R) (_p : R) (_ω : Ω) : Prop := True

/-- **Algebraically simple** — ThomasDecomp Definition 2.2, restated at
`NewMethod.tex:1236-1250` as (A1)–(A3).

(A1) triangular: no two elements share a leader, and no element lies in the
base field.
(A2) non-vanishing initials: for each equation `p` with leader `x`, the
initial of `p` is non-zero at every solution of `S_{<x}`.
(A3) square-free (see `SquareFreeAt`). -/
structure AlgebraicallySimple (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) : Prop where
  /-- (A1) — ThomasDecomp Def. 2.2(1) `|S_{x_i}| ≤ 1`. -/
  triangular : ∀ p ∈ S.eqs, ∀ q ∈ S.eqs, I.leader p = I.leader q → p = q
  /-- (A1) — ThomasDecomp Def. 2.2(1) `S ∩ {c^=, c^≠ | c ∈ F} = ∅`. -/
  noConstants : ∀ p ∈ S.eqs, I.leader p ≠ 1
  /-- (A2) — ThomasDecomp Def. 2.2(2). -/
  nonVanishingInitials :
    ∀ p ∈ S.eqs, ∀ ω ∈ I.SolBelow S (I.leader p),
      I.ev ω (I.initial p) ≠ 0
  /-- (A3) — ThomasDecomp Def. 2.2(3). -/
  squareFree :
    ∀ p ∈ S.eqs, ∀ ω ∈ I.SolBelow S (I.leader p),
      AlgorithmInput.SquareFreeAt I S p ω

/-- (S2) involutivity.

Source: ThomasDecomp Definition 3.5: "A differential system `S` is
(Janet-)involutive if all non-reductive prolongations of `(S_T)^=` reduce
to zero modulo `(S_T)^=`."  Restated at `NewMethod.tex:1229-1231`.

NOTE (definitional stub): needs Janet division's multiplicative /
non-multiplicative split (`NewMethod.tex:1215-1226`) and the differential
`Reduce` of ThomasDecomp Alg. 3.3.  Neither is in Mathlib; building them is
out of scope.  Declared debt. -/
def Involutive (_I : AlgorithmInput R Rc F Ω m k n) (_S : DiffSystem R) : Prop :=
  True

/-- (S3) `S^=` is minimal.  ThomasDecomp Definition 3.5(3).
NOTE (definitional stub): declared debt, as for `Involutive`. -/
def MinimalEqs (_I : AlgorithmInput R Rc F Ω m k n) (_S : DiffSystem R) : Prop :=
  True

/-- (S4) no inequation of `S^≠` is reducible modulo `S^=`.
ThomasDecomp Definition 3.5(4); `NewMethod.tex:1262`.
NOTE (definitional stub): declared debt, as for `Involutive`. -/
def IneqsIrreducible (_I : AlgorithmInput R Rc F Ω m k n) (_S : DiffSystem R) :
    Prop := True

/-- **Differentially simple** — ThomasDecomp Definition 3.5, restated at
`NewMethod.tex:1253-1263` as (S1)–(S4). -/
structure DifferentiallySimple (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) : Prop where
  /-- (S1) -/ algSimple : I.AlgebraicallySimple S
  /-- (S2) -/ involutive : I.Involutive S
  /-- (S3) -/ minimalEqs : I.MinimalEqs S
  /-- (S4) -/ ineqsIrreducible : I.IneqsIrreducible S

/-- A *differential Thomas decomposition*: a decomposition into finitely
many differentially simple subsystems with pairwise disjoint solution sets.

Source: ThomasDecomp Definition 2.4 / Definition 3.5 (last line): "A
disjoint decomposition of a system into differentially simple subsystems is
called (differential) Thomas decomposition."  Restated at
`NewMethod.tex:1265-1266`. -/
structure IsThomasDecomposition (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) (Ss : List (DiffSystem R)) : Prop where
  simple : ∀ T ∈ Ss, I.DifferentiallySimple T
  cover : I.Sol S = ⋃ T ∈ Ss, I.Sol T
  disjoint : ∀ T ∈ Ss, ∀ U ∈ Ss, T ≠ U → I.Sol T ∩ I.Sol U = ∅

end AlgorithmInput

/-! ############################################################
# §4.  Constructible sets, `V(·)` and `Z(·)`

Source: `NewMethod.tex:1537-1592`; CTD Definition 6, read from
`~/project/papers/chen-golubitsky-lemaire-moreno-maza-pan-2007-
comprehensive-triangular-decomposition.pdf`, §2.4.
############################################################ -/

section Constructible

variable {n : ℕ}

/-- `V(𝔞) ⊆ ℂⁿ` — the variety of an ideal of `ℚ[c₁,…,c_n]`. -/
def Vz (a : Ideal (Cst n)) : Set (Pt n) := {c | ∀ p ∈ a, evalPt c p = 0}

/-- **Constructible subset of `ℂⁿ`** — CTD, Definition 6, quoted verbatim
at `NewMethod.tex:1537-1545`: "A constructible subset of `ℂⁿ` is any finite
union `(A₁ \ B₁) ∪ … ∪ (A_e \ B_e)` where the `Aᵢ, Bᵢ` are algebraic
varieties in `ℂⁿ`." -/
def IsConstructible (s : Set (Pt n)) : Prop :=
  ∃ L : List (Ideal (Cst n) × Ideal (Cst n)),
    s = ⋃ B ∈ L, (Vz B.1 \ Vz B.2)

/-- A pair `(𝔞, 𝔟)` of ideals, the paper's representation of a locally
closed set (`NewMethod.tex:1576-1580`). -/
abbrev CPair (n : ℕ) := Ideal (Cst n) × Ideal (Cst n)

/-- `Z(𝔞,𝔟) := V(𝔞) \ V(𝔟)` (`NewMethod.tex:1578`). -/
def Zp (B : CPair n) : Set (Pt n) := Vz B.1 \ Vz B.2

/-- `Z(𝒲) := ⋃_{B ∈ 𝒲} Z(B)` (`NewMethod.tex:1580-1581`).

Following CTD §2, `Z(·)` is reserved for the passage from a representation
to the set it represents (`NewMethod.tex:1581-1583`). -/
def Zw (W : List (CPair n)) : Set (Pt n) := ⋃ B ∈ W, Zp B

lemma Zw_isConstructible (W : List (CPair n)) : IsConstructible (Zw W) :=
  ⟨W, rfl⟩

lemma Vz_span (E : Set (Cst n)) :
    Vz (Ideal.span E) = {c | ∀ p ∈ E, evalPt c p = 0} := by
  ext c
  constructor
  · intro h p hp; exact h p (Ideal.subset_span hp)
  · intro h p hp
    have hle : Ideal.span E ≤ RingHom.ker (evalHom c) := by
      rw [Ideal.span_le]
      intro q hq
      simpa [RingHom.mem_ker, evalPt] using h q hq
    have hk := hle hp
    simpa [RingHom.mem_ker, evalPt] using hk

@[simp] lemma Vz_top : Vz (⊤ : Ideal (Cst n)) = ∅ := by
  ext c
  simp only [Vz, Set.mem_setOf_eq, Set.mem_empty_iff_false, iff_false]
  intro h
  have := h 1 Submodule.mem_top
  simp at this

@[simp] lemma Vz_one : Vz (1 : Ideal (Cst n)) = ∅ := by
  rw [Ideal.one_eq_top, Vz_top]

@[simp] lemma Vz_bot : Vz (⊥ : Ideal (Cst n)) = Set.univ := by
  ext c
  simp only [Vz, Set.mem_setOf_eq, Set.mem_univ, iff_true]
  intro p hp
  rw [Ideal.mem_bot.mp hp]
  simp [evalPt]

lemma Vz_antitone {a b : Ideal (Cst n)} (h : a ≤ b) : Vz b ⊆ Vz a :=
  fun _ hc p hp => hc p (h hp)

/-- A product of ideals vanishes exactly on the union of their varieties:
`V(𝔞₁ ⋯ 𝔞_m) = ⋃_k V(𝔞_k)` (`NewMethod.tex:1584-1586`). -/
lemma Vz_mul (a b : Ideal (Cst n)) : Vz (a * b) = Vz a ∪ Vz b := by
  ext c
  constructor
  · intro h
    by_cases ha : c ∈ Vz a
    · exact Or.inl ha
    · refine Or.inr ?_
      simp only [Vz, Set.mem_setOf_eq, not_forall] at ha
      obtain ⟨p, hp, hp0⟩ := ha
      show ∀ q ∈ b, evalPt c q = 0
      intro q hq
      have hpq : evalPt c (p * q) = 0 := h _ (Ideal.mul_mem_mul hp hq)
      rw [evalPt_mul] at hpq
      rcases mul_eq_zero.mp hpq with h1 | h2
      · exact absurd h1 hp0
      · exact h2
  · rintro (h | h) p hp
    · exact h p (Submodule.mem_inf.mp (Ideal.mul_le_inf hp)).1
    · exact h p (Submodule.mem_inf.mp (Ideal.mul_le_inf hp)).2

/-- The list form of `Vz_mul`.  The empty product is `(1)` and
`V((1)) = ∅` (`NewMethod.tex:1590-1592`). -/
lemma Vz_listProd (L : List (Ideal (Cst n))) :
    Vz L.prod = ⋃ a ∈ L, Vz a := by
  induction L with
  | nil => simp
  | cons a t ih =>
      rw [List.prod_cons, Vz_mul, ih]
      ext c
      simp only [Set.mem_union, Set.mem_iUnion, exists_prop, List.mem_cons]
      constructor
      · rintro (h | ⟨i, hi, hci⟩)
        · exact ⟨a, Or.inl rfl, h⟩
        · exact ⟨i, Or.inr hi, hci⟩
      · rintro ⟨i, (rfl | hi), hci⟩
        · exact Or.inl hci
        · exact Or.inr ⟨i, hi, hci⟩

/-- Evaluation is a ring hom, so it turns a list product into a product. -/
lemma evalPt_listProd (c : Pt n) (L : List (Cst n)) :
    evalPt c L.prod = (L.map (evalPt c)).prod := by
  induction L with
  | nil => simp
  | cons a t ih => simp [List.prod_cons, evalPt_mul, ih]

lemma Vz_sup (a b : Ideal (Cst n)) : Vz (a ⊔ b) = Vz a ∩ Vz b := by
  refine Set.Subset.antisymm ?_ ?_
  · intro c hc
    exact ⟨Vz_antitone le_sup_left hc, Vz_antitone le_sup_right hc⟩
  · rintro c ⟨ha, hb⟩ p hp
    obtain ⟨y, hy, z, hz, rfl⟩ := Submodule.mem_sup.mp hp
    have h1 : evalPt c y = 0 := ha y hy
    have h2 : evalPt c z = 0 := hb z hz
    have h3 : evalPt c (y + z) = evalPt c y + evalPt c z := by simp [evalPt]
    rw [h3, h1, h2, add_zero]

end Constructible

/-! ############################################################
# §5.  Cells

Source: `NewMethod.tex:1547-1567`.
############################################################ -/

namespace AlgorithmInput

variable {R Rc F Ω : Type} [CommRing R] [CommRing Rc] [Algebra ℂ Rc]
variable [CommRing F] [Algebra ℂ F]
variable {m k n : ℕ}

/-- **The cell of a differentially simple system.**

Source (verbatim), `NewMethod.tex:1547-1552`: "The *cell* of a
differentially simple system `(S^=, S^≠)`, computed with the constant
parameters `c = (c₁,…,c_n)` ranked below all other indeterminates, is the
constant locus in `ℂⁿ` that solves the equations and inequations that lie
in `ℚ[c₁,…,c_n]`." -/
def cell (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : Set (Pt n) :=
  {c | (∀ p ∈ I.constPart S.eqs, evalPt c p = 0) ∧
       (∀ q ∈ I.constPart S.ineqs, evalPt c q ≠ 0)}

/-- `E_i` — line 3 of Algorithm 1 / line 8 of Algorithm 2. -/
def Ei (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : List (Cst n) :=
  I.constPart S.eqs

/-- `h_i` — line 4 of Algorithm 1 / line 9 of Algorithm 2: the product of
the *constant-only* inequations. -/
def hi (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : Cst n :=
  (I.constPart S.ineqs).prod

/-- `C_i = (⟨E_i⟩, ⟨h_i⟩)` — line 10 of Algorithm 2
(`NewMethod.tex:1692`); the same pair is `W_i` at line 5 of Algorithm 1. -/
def Ci (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : CPair n :=
  (Ideal.span {p : Cst n | p ∈ I.Ei S}, Ideal.span {I.hi S})

/-- `h_i` does not vanish at `c*` exactly when no constant-only inequation
of the component does. -/
lemma hi_ne_zero_iff (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R)
    (c : Pt n) :
    evalPt c (I.hi S) ≠ 0 ↔ ∀ q ∈ I.constPart S.ineqs, evalPt c q ≠ 0 := by
  rw [hi, evalPt_listProd, Ne, List.prod_eq_zero_iff]
  constructor
  · intro h q hq hq0
    exact h (List.mem_map.mpr ⟨q, hq, hq0⟩)
  · intro h hmem
    obtain ⟨q, hq, hq0⟩ := List.mem_map.mp hmem
    exact h q hq hq0

/-- **The cell is exactly the locally closed set `Z(C_i)`.**

`NewMethod.tex:1640`: "`S_i`'s cell `C_i = V(E_i) \ V(h_i)`".  This is
proved, not assumed: it is the mechanical ideal-arithmetic content of
lines 3–5 of Algorithm 1. -/
theorem cell_eq_Zp (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) :
    I.cell S = Zp (I.Ci S) := by
  have h1 : Vz (Ideal.span {p : Cst n | p ∈ I.constPart S.eqs})
      = {cc : Pt n | ∀ p ∈ I.constPart S.eqs, evalPt cc p = 0} := Vz_span _
  have h2 : Vz (Ideal.span {I.hi S})
      = {cc : Pt n | evalPt cc (I.hi S) = 0} := by
    rw [Vz_span]; ext cc; simp
  ext c
  simp only [Zp, Ci, h1, h2, Set.mem_diff, Set.mem_setOf_eq, cell, Ei]
  constructor
  · rintro ⟨he, hq⟩
    exact ⟨he, fun h0 => (I.hi_ne_zero_iff S c).mpr hq h0⟩
  · rintro ⟨he, hne⟩
    exact ⟨he, (I.hi_ne_zero_iff S c).mp hne⟩

/-- **The projection fact.**

Source: ThomasDecomp Remark 2.3, read from
`~/project/papers/bachler-gerdt-lange-hegermann-robertz-2012-differential-
thomas-decomposition.pdf`, p. 3: "if `b ∈ Sol(S_{<x})` and `S_x` is not
empty, then … each solution `b ∈ Sol(S_{<x})` extends to a solution
`(b,a) ∈ Sol(S_{≤x})` … Conversely, if `(a₁,…,a_n) ∈ Sol(S)` … then
`(a₁,…,aᵢ) ∈ Sol(S_{≤xᵢ})`."

The paper's inductive consequence, `NewMethod.tex:1554-1561`: "so,
inductively, the projection of `Sol(S^=, S^≠)` onto `ℂⁿ` is exactly the
solution set of the constant-only part of `(S^=, S^≠)`."

NOTE (departure from source): Remark 2.3 is a one-level extension /
restriction statement about `S_{≤x}`; the form used here is the paper's
inductive corollary of it, specialized to the projection onto the constant
block.  It is stated as an axiom because it is a theorem *of ThomasDecomp*,
not of this development. -/
axiom projection_eq_cell (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) (hS : I.DifferentiallySimple S) :
    I.proj '' (I.Sol S) = I.cell S

/-- **Lemma `lem:specialization`** (`NewMethod.tex:1919-1925`): a
differentially simple system stays differentially simple when specialized
at a point of its cell.

The paper proves this (`NewMethod.tex:1927-1943`) from the fact that (A2)
and (A3) are quantified over solution sets containing the whole cell.
Stated here as an axiom of the source. -/
axiom specialization_simple (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) (hS : I.DifferentiallySimple S) (cs : Pt n)
    (hcs : cs ∈ I.cell S) :
    ∀ p ∈ S.eqs, I.spec cs (I.initial p) ≠ 0 ∧ I.spec cs (I.separant p) ≠ 0

end AlgorithmInput

/-! ############################################################
# §6.  The four subroutines

Table at `NewMethod.tex:1876-1911`.  Each is an opaque constant plus a
contract quoting its source's own numbering.  None is implemented here:
the paper takes all four from the literature and uses them unmodified
(`NewMethod.tex:1874-1876`).
############################################################ -/

namespace AlgorithmInput

variable {R Rc F Ω : Type} [CommRing R] [CommRing Rc] [Algebra ℂ Rc]
variable [CommRing F] [Algebra ℂ F]
variable {m k n : ℕ}

/-- **`DifferentialDecompose`** — a differential Thomas decomposition.

Source: ThomasDecomp, Algorithm 3.6 (`NewMethod.tex:1884`), read from
`~/project/papers/bachler-gerdt-lange-hegermann-robertz-2012-differential-
thomas-decomposition.pdf`, p. 22: "Input: A differential system `S′` with
`(S′)_T = ∅`.  Output: A differential Thomas decomposition of `S′`." -/
axiom DifferentialDecompose (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) : List (DiffSystem R)

/-- The contract of `DifferentialDecompose`: its output is a differential
Thomas decomposition of its input (ThomasDecomp Alg. 3.6, Output clause,
together with Definition 2.4 for what "decomposition" means). -/
axiom differentialDecompose_spec (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) :
    I.IsThomasDecomposition S (I.DifferentialDecompose S)

/-- **`FullReduce`** — Ritt's full reduction: the remainder `r`, fully
reduced with respect to the equations of the system.

Source: Ritt, *Differential Algebra*, AMS Colloq. Publ. 33 (1950), §I.6
("REDUCTION"), read via `djvutxt` from `~baccala/Books/Differential
Algebra/Ritt J.F. Differential algebra (AMS, 1950)(T)(189s).djvu`, p. 5.
Cited by the paper at `NewMethod.tex:1888`. -/
axiom FullReduce (I : AlgorithmInput R Rc F Ω m k n) (P : R) (S : List R) : R

/-- **Ritt's reduction identity**, Ritt §I.6, p. 5, quoted:

"Let `G` be any d.p.  There exist nonnegative integers `s_i, t_i`,
`i = 1,…,r`, such that, when a suitable linear combination of the `A_i`,
and of a certain number of their derivatives, with d.p. for coefficients,
is subtracted from `S₁^{s₁}⋯S_r^{s_r} I₁^{t₁}⋯I_r^{t_r} G`, the remainder
`R` is reduced with respect to `(4)`."

Restated by the paper at `NewMethod.tex:1148-1152` and used as
eq. `ritt reduction 1` at `NewMethod.tex:1963`: `h·P − r ∈ [S]`.

NOTE (departure from source): Ritt's `h` is a specific power product of the
`Sᵢ` and `Iᵢ`; here we only record that `h` lies in the multiplicative
closure of `H_S = I_S ∪ S_S`, which is what the paper's proof of
Thm. `completeness` uses (it needs only `h(c*) ≠ 0`). -/
axiom fullReduce_spec (I : AlgorithmInput R Rc F Ω m k n) (P : R) (S : List R) :
    ∃ h : R, h ∈ Submonoid.closure {y : R | y ∈ I.initSep S} ∧
      h * P - I.FullReduce P S ∈ diffIdealBy I.d {y : R | y ∈ S}

/-- **`minAss`** — the minimal primes of an ideal of `ℚ[c₁,…,c_n]`,
presented by generating sets; `minAss((1)) = ∅`.

Source: Gianni–Trager–Zacharias, *Gröbner Bases and Primary Decomposition
of Polynomial Ideals*, JSC 6 (1988), §9 ("Applications to Computing
Radicals and Associated Primes"), read from
`~baccala/Books/Algebra/Primary Decomposition/GianniTragerZacharias.pdf`.
Cited by the paper at `NewMethod.tex:1893-1896` as "Cor. 3.2(v), §9".

NOTE (suspected error in source): GTZ Corollary 3.2(v) is *not* about
minimal primes — it says that `I R[x]_f ∩ R[x]` is computable for a
nonzerodivisor `f`, i.e. it is the *saturation* primitive `I : f^∞`.  The
minimal-primes content is §9.  Cor. 3.2(v) is the right citation for the
`𝔞 : 𝔟^∞` of Algorithm 2 line 17, so the two citations appear to have been
attached to the wrong subroutine.

The Mathlib notion this corresponds to is `Ideal.minimalPrimes`. -/
axiom minAss (a : Ideal (Cst n)) : List (Ideal (Cst n))

/-- `minAss` computes exactly `Ideal.minimalPrimes`. -/
axiom minAss_spec (a : Ideal (Cst n)) :
    {p : Ideal (Cst n) | p ∈ minAss a} = a.minimalPrimes

/-- `minAss((1)) = ∅` (`NewMethod.tex:1895`).  This is *forced* by
`minAss_spec`, since no prime contains the unit ideal — proved below as
`minAss_top`. -/
lemma minimalPrimes_top : (⊤ : Ideal (Cst n)).minimalPrimes = ∅ := by
  ext p
  simp only [Set.mem_empty_iff_false, iff_false]
  rintro ⟨⟨hp, hle⟩, -⟩
  exact hp.ne_top (top_le_iff.mp hle)

lemma minAss_top : minAss (⊤ : Ideal (Cst n)) = [] := by
  have h := minAss_spec (⊤ : Ideal (Cst n))
  rw [minimalPrimes_top] at h
  cases hm : minAss (⊤ : Ideal (Cst n)) with
  | nil => rfl
  | cons a t =>
      exfalso
      have : a ∈ {p : Ideal (Cst n) | p ∈ minAss (⊤ : Ideal (Cst n))} := by
        rw [hm]; exact List.mem_cons_self ..
      rw [h] at this
      exact this

/-- `𝔞 : 𝔟^∞` for an ideal `𝔟` — needed by Algorithm 2 line 17
(`NewMethod.tex:1706`). -/
def colonIdeal {A : Type} [CommRing A] (I J : Ideal A) : Ideal A where
  carrier := {y | ∀ z ∈ J, z * y ∈ I}
  zero_mem' := by intro z _; simp
  add_mem' := by
    intro a b ha hb z hz
    have := I.add_mem (ha z hz) (hb z hz)
    simpa [mul_add] using this
  smul_mem' := by
    intro t a ha z hz
    have := I.mul_mem_left t (ha z hz)
    have h : t * (z * a) = z * (t • a) := by simp [smul_eq_mul]; ring
    rwa [h] at this

/-- `𝔞 : 𝔟^∞ = ⋃_k (𝔞 : 𝔟^k)`. -/
def satByIdeal {A : Type} [CommRing A] (I J : Ideal A) : Ideal A :=
  ⨆ k : ℕ, colonIdeal I (J ^ k)

/-- **`RefiningPartition`** — for a finite family of constructible sets, a
finite set of pairs `(B, ι)` in which the `B` are pairwise disjoint, `ι` is
the subset of the family `B` is contained in, and every member of the
family is the union of the `B` with that member in their `ι`.

Source: Chen–Lemaire–Li–Moreno Maza–Pan–Xie, *The ConstructibleSetTools and
ParametricSystemTools modules of the RegularChains library in Maple*,
ICCSA 2008, Example 8, read from
`~/project/papers/chen-lemaire-li-moreno-maza-pan-xie-2008-
constructiblesettools-regularchains.pdf`: "there is a set-theoretical
co-prime factorization problem: constructing another finite list of
pairwise disjoint constructible sets `out_lcs` out of `in_lcs` such that
every constructible set in `in_lcs` can be uniquely written as a union of
several constructible sets in `out_lcs`.  This task is achieved by the
command `RefiningPartition`."  Its output "is represented by a matrix in
which the first column are constructible sets and the second column are
indices showing where the constructible sets come from."

NOTE (suspected error in source): `NewMethod.tex:1897` cites "\cite{CSTools},
§3" for `RefiningPartition`, but §3 of that paper is the
`ParametricSystemTools` module; `RefiningPartition` is documented in §2.4
(`ConstructibleSetTools`, Example 8).

NOTE (departure from source): the paper's tabulated description
(`NewMethod.tex:1897-1902`) says the `B` are pairwise disjoint and `ι` is
"the subset of the family `B` is contained in"; CSTools only claims a
pairwise-disjoint refinement with provenance indices.  We formalize the
paper's stronger reading, since Algorithm 2 line 16 relies on it. -/
axiom RefiningPartition (family : List (CPair n)) : List (CPair n × List ℕ)

/-- Contract of `RefiningPartition`, part 1: the blocks are pairwise
disjoint. -/
axiom refiningPartition_disjoint (family : List (CPair n)) :
    ∀ B ∈ RefiningPartition family, ∀ B' ∈ RefiningPartition family,
      B ≠ B' → Zp B.1 ∩ Zp B'.1 = ∅

/-- Contract of `RefiningPartition`, part 2: every member of the family is
the union of the blocks whose index set contains it. -/
axiom refiningPartition_cover (family : List (CPair n)) (j : ℕ)
    (hj : j < family.length) :
    Zp (family.get ⟨j, hj⟩) =
      ⋃ B ∈ (RefiningPartition family).filter (fun B => B.2.contains j), Zp B.1

/-- Contract of `RefiningPartition`, part 3: `ι` really records containment
— a block is contained in exactly the members its `ι` names. -/
axiom refiningPartition_index (family : List (CPair n)) (j : ℕ)
    (hj : j < family.length) (B : CPair n × List ℕ)
    (hB : B ∈ RefiningPartition family) :
    (B.2.contains j = true) ↔ Zp B.1 ⊆ Zp (family.get ⟨j, hj⟩)

end AlgorithmInput

/-! ############################################################
# §7.  The three constant loci

Source: `NewMethod.tex:1327-1500` (`sec:constant-loci`).
############################################################ -/

namespace AlgorithmInput

variable {R Rc F Ω : Type} [CommRing R] [CommRing Rc] [Algebra ℂ Rc]
variable [CommRing F] [Algebra ℂ F]
variable {m k n : ℕ}

/-- The specialized differential ideal `[L(c*)] ⊆ Rc`. -/
def specDiffIdeal (I : AlgorithmInput R Rc F Ω m k n) (cs : Pt n) (L : List R) :
    Ideal Rc :=
  diffIdealBy I.dc {y : Rc | ∃ p ∈ L, I.spec cs p = y}

/-- The saturated ideal `[L(c*)] : Q(c*)^∞`. -/
def satIdealOf (I : AlgorithmInput R Rc F Ω m k n) (cs : Pt n) (L : List R) :
    Ideal Rc :=
  satBy (I.specDiffIdeal cs L) (I.spec cs I.Q)

/-- **The non-degeneracy locus**
`N_Q = {c* ∈ ℂⁿ : 1 ∉ [A(c*)] : Q^∞}` — `NewMethod.tex:1360-1362`,
eq. `nondegeneracy`. -/
def Nq (I : AlgorithmInput R Rc F Ω m k n) : Set (Pt n) :=
  {cs | Consistent (I.satIdealOf cs I.ansatz)}

/-- **The consistency locus**
`V_∃ = {c* ∈ ℂⁿ : 1 ∉ [A(c*), 𝒫^=] : Q^∞}` — `NewMethod.tex:1371-1373`. -/
def Vexists (I : AlgorithmInput R Rc F Ω m k n) : Set (Pt n) :=
  {cs | Consistent (I.satIdealOf cs (I.ansatz ++ I.Peq))}

/-- **The membership locus**
`V_∀ = {c* ∈ ℂⁿ : 1 ∉ [A(c*)] : Q^∞ ∧ 𝒫^= ⊆ √([A(c*)] : Q^∞)}` —
`NewMethod.tex:1387-1389`. -/
def Vforall (I : AlgorithmInput R Rc F Ω m k n) : Set (Pt n) :=
  {cs | Consistent (I.satIdealOf cs I.ansatz) ∧
        ∀ P ∈ I.Peq, I.spec cs P ∈ (I.satIdealOf cs I.ansatz).radical}

/-- **Theorem `thm:locus-containment`** (`NewMethod.tex:1404-1412`):
`V_∀ ⊆ V_∃`.

The paper's proof (`NewMethod.tex:1414-1432`) is followed step for step:
`√([A(c*)] : Q^∞)` is a differential ideal containing both `A(c*)` and
`𝒫^=`, hence contains `[A(c*), 𝒫^=]`; saturating that containment by
`Q^∞` and using `(I : Q^∞) : Q^∞ = I : Q^∞` gives the result. -/
theorem Vforall_subset_Vexists (I : AlgorithmInput R Rc F Ω m k n) :
    I.Vforall ⊆ I.Vexists := by
  rintro cs ⟨hcons, hmem⟩
  have hdiffJ : IsDiffIdeal I.dc (I.specDiffIdeal cs I.ansatz) :=
    diffIdealBy_isDiffIdeal
  have hdiffSat : IsDiffIdeal I.dc (I.satIdealOf cs I.ansatz) :=
    satBy_isDiffIdeal I.dc_isDeriv hdiffJ _
  have hdiffRad : IsDiffIdeal I.dc (I.satIdealOf cs I.ansatz).radical :=
    radical_isDiffIdeal I.dc_isDeriv hdiffSat
  -- `[A(c*), 𝒫^=] ≤ √([A(c*)] : Q^∞)`
  have hgen : {y : Rc | ∃ p ∈ I.ansatz ++ I.Peq, I.spec cs p = y}
      ⊆ ((I.satIdealOf cs I.ansatz).radical : Set Rc) := by
    rintro y ⟨p, hp, rfl⟩
    rcases List.mem_append.mp hp with hpa | hpe
    · exact Ideal.le_radical (le_satBy (subset_diffIdealBy ⟨p, hpa, rfl⟩))
    · exact hmem p hpe
  have hle : I.specDiffIdeal cs (I.ansatz ++ I.Peq)
      ≤ (I.satIdealOf cs I.ansatz).radical := diffIdealBy_le hgen hdiffRad
  have hmono := satBy_mono (q := I.spec cs I.Q) hle
  -- `(I : Q^∞) : Q^∞ = I : Q^∞`, applied to the radical
  have hrad : satBy ((I.satIdealOf cs I.ansatz).radical) (I.spec cs I.Q)
      = (I.satIdealOf cs I.ansatz).radical := by
    simp only [satIdealOf]
    exact satBy_radical_satBy
  rw [hrad] at hmono
  have hradcons : Consistent (I.satIdealOf cs I.ansatz).radical := by
    intro h
    obtain ⟨p, hp⟩ := h
    rw [one_pow] at hp
    exact hcons hp
  exact consistent_of_le hmono hradcons

end AlgorithmInput

/-! ############################################################
# §8.  The two algorithms

Every step is annotated with the pseudocode line number of the
corresponding `\item` in `NewMethod.tex`.
############################################################ -/

namespace AlgorithmInput

variable {R Rc F Ω : Type} [CommRing R] [CommRing Rc] [Algebra ℂ Rc]
variable [CommRing F] [Algebra ℂ F]
variable {m k n : ℕ}

/-! ## §8.1  Algorithm 1 — `ConsistencyLocus` (`NewMethod.tex:1594`) -/

/-- **Algorithm 1, `ConsistencyLocus`** — `NewMethod.tex:1594-1626`.

```
1  (S₁,…,S_s) ← DifferentialDecompose((A(c) ∪ 𝒫^=, {Q}), ≺)
2  for i = 1,…,s do
3     E_i ← S_i^= ∩ ℚ[c₁,…,c_n]
4     h_i ← ∏ (S_i^≠ ∩ ℚ[c₁,…,c_n])
5     W_i ← { (⟨E_i⟩, ⟨h_i⟩) }
6  end for
7  return ⋃_{i=1}^{s} W_i
```

Each `W_i` is a singleton, so the union of line 7 is the list of the
pairs. -/
def consistencyLocus (I : AlgorithmInput R Rc F Ω m k n) : List (CPair n) :=
  -- line 1
  let Ss := I.DifferentialDecompose ⟨I.ansatz ++ I.Peq, [I.Q]⟩
  -- lines 2–6; lines 3, 4, 5 are `Ei`, `hi`, `Ci`
  Ss.map (fun S => I.Ci S)
  -- line 7

/-- Output clause of Algorithm 1 (`NewMethod.tex:1615-1619`): the returned
representation is constructible. -/
theorem consistencyLocus_isConstructible (I : AlgorithmInput R Rc F Ω m k n) :
    IsConstructible (Zw I.consistencyLocus) :=
  Zw_isConstructible _

/-- `Z(𝒲) = V_∃` — the Output clause of Algorithm 1
(`NewMethod.tex:1618-1619`).

The paper's justification is at `NewMethod.tex:1628-1641`: each component
`S_i` is a differentially simple specialization of `A(c) ∪ 𝒫^=` together
with `Q`; by the cell argument (ThomasDecomp Remark 2.3) every `c*` in
`S_i`'s cell extends to a genuine solution of the whole specialized
system, "so `c* ∈ V_∃` exactly when it lies in the cell of *some*
component, and `V_∃ = ⋃_{i=1}^{s} C_i`."

`sorry` classification: **(b) cited theorem** — the paper proves it, from
`projection_eq_cell` plus the differential Nullstellensatz linking
`1 ∉ [·] : Q^∞` to the existence of a solution.  That last link is the step
the paper leaves implicit (`NewMethod.tex:1381-1384` explicitly declines to
assert existence of solutions in a function space), so mechanizing it would
require choosing a Nullstellensatz; out of scope here. -/
theorem consistencyLocus_correct (I : AlgorithmInput R Rc F Ω m k n) :
    Zw I.consistencyLocus = I.Vexists := by
  sorry

/-! ## §8.2  Algorithm 2 — `MembershipLocus` (`NewMethod.tex:1657`) -/

/-- The flag `μ ∈ {∀, ∃∀}` selecting the assembly
(`NewMethod.tex:1678-1679`).

Modelled as an inductive rather than a `Bool`, so that the two branches
carry their names. -/
inductive Assembly where
  | forAll        -- `μ = ∀`
  | existsForAll  -- `μ = ∃∀`
  deriving DecidableEq

/-- `J_i` — lines 3–7 of Algorithm 2 (`NewMethod.tex:1685-1689`):

```
3     J_i ← ∅
4     for j = 1,…,t do
5        r ← FullReduce(P_j, S_i^=)
6        J_i ← J_i ∪ Coeffs(r, ℚ[c₁,…,c_n])
7     end for
```
-/
def Ji (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : List (Cst n) :=
  I.Peq.foldr (fun P acc => I.coeffs (I.FullReduce P S.eqs) ++ acc) []

/-- `𝔥_i` — line 11 of Algorithm 2 (`NewMethod.tex:1693`):
`𝔥_i ← ∏_{q ∈ S_i^≠} ⟨Coeffs(q, ℚ[c₁,…,c_n])⟩`.

Note this is the product over *all* inequations of the component, of the
ideal generated by the constant-coefficients of each; it is a different
object from `h_i` (`hi`), which is the product of only the
*already-constant* inequations.  The paper's argument turns on the
distinction (`NewMethod.tex:1741-1760`).

NOTE (suspected error in source): the surrounding prose
(`NewMethod.tex:1738-1750`) defines `𝔥_i = ∏_{𝔮 ∈ K_i} 𝔮` where `K_i`
collects, over every `q ∈ S_i^≠`, the *minimal associated primes* of
`⟨Coeffs(q, ℚ[c])⟩` — and further says primes `𝔭` with some `𝔮 ∈ K_i`
contained in them are "discarded outright".  Neither `K_i` nor the
discarding step occurs in the pseudocode, which forms the product of the
coefficient ideals themselves.  `V(∏ 𝔮 over minimal primes) = V(∏ of the
ideals)` because `V(𝔞) = V(√𝔞) = ⋃ V(𝔭)`, so `Z(W_i)` is unaffected; but
the pseudocode and the prose are not the same computation, and the
`K_i`-based discarding is genuinely absent from the algorithm.  We follow
the pseudocode. -/
def hfrak (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : Ideal (Cst n) :=
  (S.ineqs.map (fun q => Ideal.span {p : Cst n | p ∈ I.coeffs q})).prod

/-- `W_i = (⟨J_i ∪ E_i⟩, 𝔥_i)` — line 12 of Algorithm 2
(`NewMethod.tex:1694`). -/
def Wi (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) : CPair n :=
  (Ideal.span {p : Cst n | p ∈ I.Ji S ++ I.Ei S}, I.hfrak S)

/-- The family handed to `RefiningPartition` at line 15
(`NewMethod.tex:1704`): `{((0),(1)), W₁, C₁, …, W_s, C_s}`.

We fix the concrete indexing used by `keepBlock` below: `((0),(1))` sits at
index `0`, and for `i = 0,…,s−1`, `W_{i+1}` sits at index `2i+1` and
`C_{i+1}` at index `2i+2`. -/
def refiningFamily (WC : List (CPair n × CPair n)) : List (CPair n) :=
  ((⊥ : Ideal (Cst n)), (⊤ : Ideal (Cst n))) ::
    WC.foldr (fun p acc => p.1 :: p.2 :: acc) []

/-- Line 16 of Algorithm 2 (`NewMethod.tex:1705`):
`D ← { B : (B,ι) ∈ ℬ, C_i ∈ ι ⟹ W_i ∈ ι for every i }`,
under the indexing fixed by `refiningFamily`. -/
def keepBlock (s : ℕ) (Bi : CPair n × List ℕ) : Bool :=
  (List.range s).all (fun i => ! (Bi.2.contains (2 * i + 2)) || Bi.2.contains (2 * i + 1))

/-- **Algorithm 2, `MembershipLocus`** — `NewMethod.tex:1657-1712`.

```
 1  (S₁,…,S_s) ← DifferentialDecompose((A(c), {Q}), ≺)
 2  for i = 1,…,s do
 3     J_i ← ∅
 4     for j = 1,…,t do
 5        r ← FullReduce(P_j, S_i^=)
 6        J_i ← J_i ∪ Coeffs(r, ℚ[c₁,…,c_n])
 7     end for
 8     E_i ← S_i^= ∩ ℚ[c₁,…,c_n]
 9     h_i ← ∏ (S_i^≠ ∩ ℚ[c₁,…,c_n])
10     C_i ← (⟨E_i⟩, ⟨h_i⟩)
11     𝔥_i ← ∏_{q ∈ S_i^≠} ⟨Coeffs(q, ℚ[c₁,…,c_n])⟩
12     W_i ← { (⟨J_i ∪ E_i⟩, 𝔥_i) }
13  end for
14  if μ = ∀ then
15     ℬ ← RefiningPartition({((0),(1)), W₁, C₁, …, W_s, C_s})
16     D ← { B : (B,ι) ∈ ℬ, C_i ∈ ι ⟹ W_i ∈ ι for every i }
17     return { (𝔭, 𝔟) : (𝔞,𝔟) ∈ D, 𝔭 ∈ minAss(𝔞 : 𝔟^∞) }
18  else
19     return ⋃_{i=1}^{s} W_i
20  end if
```
-/
def membershipLocus (I : AlgorithmInput R Rc F Ω m k n) (μ : Assembly) :
    List (CPair n) :=
  -- line 1
  let Ss := I.DifferentialDecompose ⟨I.ansatz, [I.Q]⟩
  -- lines 2–13
  let WC : List (CPair n × CPair n) := Ss.map (fun S => (I.Wi S, I.Ci S))
  match μ with
  | Assembly.existsForAll =>
      -- line 19
      WC.map Prod.fst
  | Assembly.forAll =>
      -- line 15
      let B := RefiningPartition (refiningFamily WC)
      -- line 16
      let D := B.filter (keepBlock WC.length)
      -- line 17
      D.foldr
        (fun bi acc =>
          (minAss (satByIdeal bi.1.1 bi.1.2)).map (fun p => (p, bi.1.2)) ++ acc)
        []

/-- Output clause of Algorithm 2 (`NewMethod.tex:1670-1677`): the returned
representation is constructible, for either flag. -/
theorem membershipLocus_isConstructible (I : AlgorithmInput R Rc F Ω m k n)
    (μ : Assembly) : IsConstructible (Zw (I.membershipLocus μ)) :=
  Zw_isConstructible _

/-- The intermediate locus `V_{∃∀}` (`NewMethod.tex:1673-1675`), defined as
the paper defines it at `NewMethod.tex:1866-1872`:
`V_{∃∀} = ⋃_i Z(W_i)`. -/
def VexistsForall (I : AlgorithmInput R Rc F Ω m k n) : Set (Pt n) :=
  ⋃ S ∈ I.DifferentialDecompose (⟨I.ansatz, [I.Q]⟩ : DiffSystem R), Zp (I.Wi S)

/-- The paper's assembly formula for `V_∀` (`NewMethod.tex:1866-1870`):
`V_∀ = ⋂_i (Z(W_i) ∪ (ℂⁿ \ Z(C_i)))`. -/
def VforallAssembled (I : AlgorithmInput R Rc F Ω m k n) : Set (Pt n) :=
  ⋂ S ∈ I.DifferentialDecompose (⟨I.ansatz, [I.Q]⟩ : DiffSystem R),
    (Zp (I.Wi S) ∪ (Set.univ \ Zp (I.Ci S)))

/-- **Corollary `cor:assembly`** (`NewMethod.tex:2100-2110`), set-theoretic
content: "A point lies in the membership locus when it lies in
`𝒱_i ∩ N_i` for *every* component whose cell contains it"
(`NewMethod.tex:1871-1873`).

This is proved, not assumed. -/
theorem mem_VforallAssembled_iff (I : AlgorithmInput R Rc F Ω m k n)
    (cs : Pt n) :
    cs ∈ I.VforallAssembled ↔
      ∀ S ∈ I.DifferentialDecompose (⟨I.ansatz, [I.Q]⟩ : DiffSystem R),
        cs ∈ Zp (I.Ci S) → cs ∈ Zp (I.Wi S) := by
  simp only [VforallAssembled, Set.mem_iInter₂, Set.mem_union, Set.mem_diff,
    Set.mem_univ, true_and]
  constructor
  · intro h S hS hC
    rcases h S hS with h1 | h2
    · exact h1
    · exact absurd hC h2
  · intro h S hS
    by_cases hC : cs ∈ Zp (I.Ci S)
    · exact Or.inl (h S hS hC)
    · exact Or.inr hC

/-- **`V_∀ ⊆ V_{∃∀}`** — the left half of eq. `intermediate-sandwich`
(`NewMethod.tex:1845`), under the hypothesis the paper attaches to it:
"whenever the ansatz is consistent throughout `ℂⁿ`"
(`NewMethod.tex:1675-1677`), which is exactly the statement that the cells
cover `ℂⁿ`.

NOTE: the hypothesis is *necessary*, and this Lean statement makes that
visible.  Without it, a point lying in no cell belongs to
`⋂_i (Z(W_i) ∪ (ℂⁿ \ Z(C_i)))` vacuously but to no `Z(W_i)`, so the
containment fails.  The paper states the hypothesis but does not remark
that it cannot be dropped. -/
theorem VforallAssembled_subset_VexistsForall
    (I : AlgorithmInput R Rc F Ω m k n)
    (hcover : ∀ cs : Pt n,
      ∃ S ∈ I.DifferentialDecompose (⟨I.ansatz, [I.Q]⟩ : DiffSystem R),
        cs ∈ Zp (I.Ci S)) :
    I.VforallAssembled ⊆ I.VexistsForall := by
  intro cs hcs
  obtain ⟨S, hS, hC⟩ := hcover cs
  have := (I.mem_VforallAssembled_iff cs).mp hcs S hS hC
  simp only [VexistsForall, Set.mem_iUnion₂]
  exact ⟨S, hS, this⟩

/-- **`V_{∃∀} ⊆ V_∃`** — the right half of eq. `intermediate-sandwich`
(`NewMethod.tex:1845`).

`sorry` classification: **(b) cited theorem** — asserted by the paper at
`NewMethod.tex:1838-1846` and explained at `NewMethod.tex:1826-1836`, but
not given a separate proof there. -/
theorem VexistsForall_subset_Vexists (I : AlgorithmInput R Rc F Ω m k n) :
    I.VexistsForall ⊆ I.Vexists := by
  sorry

/-- `Z(membershipLocus … ∃∀) = V_{∃∀}`.

This is immediate from the definitions: line 19 returns exactly the `W_i`
and `V_{∃∀}` is defined as `⋃ Z(W_i)`. -/
theorem membershipLocus_existsForAll_correct
    (I : AlgorithmInput R Rc F Ω m k n) :
    Zw (I.membershipLocus Assembly.existsForAll) = I.VexistsForall := by
  simp [membershipLocus, VexistsForall, Zw, Wi]

/-- `Z(membershipLocus … ∀) = V_∀` — the Output clause of Algorithm 2 for
`μ = ∀` (`NewMethod.tex:1672-1673`).

`sorry` classification: **(b) cited theorem** — this is
Corollary `cor:assembly` (`NewMethod.tex:2093-2110`) combined with the
`RefiningPartition` contract and the fact that
`Z(𝔞,𝔟) = ⋃_{𝔭 ∈ minAss(𝔞:𝔟^∞)} Z(𝔭,𝔟)` (line 17).  The paper proves the
first ingredient; the second is standard (minimal primes cut out the
irreducible components, `NewMethod.tex:1815-1822`, citing
Cox–Little–O'Shea Cor. 4.5.4).  Mechanizing it needs the Nullstellensatz
link `V(𝔞) = ⋃ V(𝔭)`, which over `ℚ[c] → ℂⁿ` is real content and out of
scope here. -/
theorem membershipLocus_forAll_correct (I : AlgorithmInput R Rc F Ω m k n) :
    Zw (I.membershipLocus Assembly.forAll) = I.Vforall := by
  sorry

/-- The paper's identification `Z(W_i) = 𝒱_i ∩ Z(C_i) ∩ N_i`
(`NewMethod.tex:1855-1865`), where `N_i = ℂⁿ \ V(𝔥_i)`.

`sorry` classification: **(b) cited theorem** — the paper argues it at
`NewMethod.tex:1858-1865`, resting on "`V(h_i) ⊆ V(𝔥_i)` since `h_i`'s own
factors are among `K_i`'s".  Since `K_i` is absent from the pseudocode (see
`hfrak`), the Lean statement uses the pseudocode's `𝔥_i`; `V(h_i) ⊆
V(𝔥_i)` still holds because an already-constant inequation `q` has
`Coeffs(q, ℚ[c]) = {q}`. -/
theorem Zp_Wi_eq (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R) :
    Zp (I.Wi S)
      = Vz (Ideal.span {p : Cst n | p ∈ I.Ji S}) ∩ Zp (I.Ci S)
        ∩ (Set.univ \ Vz (I.hfrak S)) := by
  sorry

end AlgorithmInput

/-! ############################################################
# §9.  Bridge to `~/project/Completeness.lean`

The declarations below are reproduced from `~/project/Completeness.lean`
(the earlier formalization of `sec:completeness`, `NewMethod.tex:1911`), so
that the two developments meet in one file.  Theorems 1 and 2 there are
*not* re-proved; they are restated and pointed at.
############################################################ -/

section CompletenessBridge

universe uu vv

/-- A differential ring: a commutative ring with a finite family of
commuting derivations indexed by `ι`.  From `~/project/Completeness.lean`. -/
class DifferentialRing (R : Type uu) (ι : Type vv) extends CommRing R where
  deriv      : ι → R → R
  deriv_add  : ∀ i x y, deriv i (x + y) = deriv i x + deriv i y
  deriv_mul  : ∀ i x y, deriv i (x * y) = deriv i x * y + x * deriv i y
  deriv_comm : ∀ i j x, deriv i (deriv j x) = deriv j (deriv i x)

variable {K : Type uu} [Field K]
variable {R : Type uu} {ι : Type vv} [DifferentialRing R ι] [Algebra K R]

/-- `[S]` — smallest ideal of `R` containing `S` and closed under every
derivation.  From `~/project/Completeness.lean`, where it was `sorry`;
here it delegates to the real definition of §1. -/
def diffIdeal (S : Set R) : Ideal R :=
  sInf {J : Ideal R | S ⊆ (J : Set R) ∧
    ∀ (i : ι) (x : R), x ∈ J → DifferentialRing.deriv i x ∈ J}

/-- Triangular set of differential polynomials with leader/initial/
separant data attached.  From `~/project/Completeness.lean`. -/
structure TriangularSet (R : Type uu) [CommRing R] where
  elems    : List R
  leader   : R → R
  initial  : R → R
  separant : R → R

namespace TriangularSet

/-- Submonoid generated by the separants of `A`. -/
def sepMonoid (A : TriangularSet R) : Submonoid R :=
  Submonoid.closure {s | ∃ a ∈ A.elems, s = A.separant a}

/-- Algebraic ideal `(A)`. -/
def algIdeal (A : TriangularSet R) : Ideal R :=
  Ideal.span {a | a ∈ A.elems}

/-- Differential ideal `[A]`. -/
def diffIdealOf (A : TriangularSet R) : Ideal R :=
  diffIdeal (ι := ι) {a | a ∈ A.elems}

/-- No derivative of any leader of `A` appears in `r`.
NOTE (definitional stub), as in `~/project/Completeness.lean`. -/
def PartiallyReduced (A : TriangularSet R) (r : R) : Prop := sorry

/-- Partially reduced, and for each leader `uᵢ`,
`deg_{uᵢ}(r) < deg(aᵢ, uᵢ)`.  NOTE (definitional stub). -/
def FullyReduced (A : TriangularSet R) (r : R) : Prop := sorry

/-- `A` is a regular differential system (blop, Definition 22).
NOTE (definitional stub). -/
def Regular (A : TriangularSet R) : Prop := sorry

end TriangularSet

/-- Saturation `I : M^∞` by a submonoid.  From
`~/project/Completeness.lean`. -/
def satByMonoid (I : Ideal R) (M : Submonoid R) : Ideal R where
  carrier   := {x | ∃ m ∈ M, m * x ∈ I}
  zero_mem' := ⟨1, M.one_mem, by simp⟩
  add_mem'  := by
    rintro x y ⟨mm, hm, hmx⟩ ⟨nn, hn, hny⟩
    refine ⟨mm * nn, M.mul_mem hm hn, ?_⟩
    have h : mm * nn * (x + y) = nn * (mm * x) + mm * (nn * y) := by ring
    rw [h]
    exact I.add_mem (I.mul_mem_left nn hmx) (I.mul_mem_left mm hny)
  smul_mem' := by
    rintro c x ⟨mm, hm, hmx⟩
    refine ⟨mm, hm, ?_⟩
    rw [smul_eq_mul, ← mul_assoc, mul_comm mm c, mul_assoc]
    exact I.mul_mem_left c hmx

/-- **Ritt's full reduction** (`~/project/Completeness.lean`; Ritt §I.6). -/
axiom rittReduce (P : R) (A : TriangularSet R) :
    ∃ h r : R, A.FullyReduced r ∧ h * P - r ∈ A.diffIdealOf (ι := ι)

/-- **Rosenfeld's Lemma** (blop Theorem 23), from
`~/project/Completeness.lean`.  Source: Boulier–Lazard–Ollivier–Petitot,
*Computing representations for radicals of finitely generated differential
ideals*, AAECC 20 (2009),
`~/project/papers/boulier-lazard-ollivier-petitot-2009-radical-differential-
ideal.pdf`. -/
axiom rosenfeld {A : TriangularSet R} {f : R} :
    A.Regular → A.PartiallyReduced f →
    f ∈ satByMonoid (A.diffIdealOf (ι := ι)) A.sepMonoid →
    f ∈ satByMonoid A.algIdeal A.sepMonoid

/-- **The Lazard gap lemma**, from `~/project/Completeness.lean`.

`sorry` classification: **(c) genuine gap** — the published argument
(`NewMethod.tex`'s `sec:completeness`, and `~/project/proof.tex §Gap in the
Lazard step`) claims fully-reduced monomials form a `K(N)`-basis of each
`Fᵢ = Frac(R/pᵢ)`; that is false — they span but need not be independent.
Counterexample `A = {u³ − u²}`, `r = u − 1` (which fails only regularity;
whether regularity restores the conclusion is open). -/
lemma lazard_gap
    {A : TriangularSet R} {r : R}
    (hReg  : A.Regular)
    (hRed  : A.FullyReduced r)
    (hrSat : r ∈ satByMonoid A.algIdeal A.sepMonoid) :
    r = 0 := sorry

/-- **Core completeness lemma** (`ritt_remainder_zero` of
`~/project/Completeness.lean`, = the kernel of Theorem
`thm:completeness`, `NewMethod.tex:1946-1994`).

Restated, not re-proved: see `~/project/Completeness.lean` for the proof
modulo the three blop axioms and `lazard_gap`.

The link to §8: Theorem `thm:completeness` is what makes
`membershipLocus`'s inner loop (lines 5–6) correct — it says that for `c*`
in the cell, `c* ∈ V(Coeffs(r, ℚ[c]))` iff `P ∈ [S(c*)] : q(c*)^∞`, i.e.
that `J_i` cuts out exactly the constants at which `P_j` is a member. -/
axiom ritt_remainder_zero
    {A : TriangularSet R} {P h r : R}
    (hReg  : A.Regular)
    (hRed  : A.FullyReduced r)
    (hRitt : h * P - r ∈ A.diffIdealOf (ι := ι))
    (hP    : P ∈ A.diffIdealOf (ι := ι)) :
    r = 0

end CompletenessBridge

/-! ############################################################
# §10.  Proposition `simple-membership` and the decomposition
        proposition

These are the two facts about simple systems the completeness argument
rests on (`NewMethod.tex:1276-1324`).
############################################################ -/

namespace AlgorithmInput

variable {R Rc F Ω : Type} [CommRing R] [CommRing Rc] [Algebra ℂ Rc]
variable [CommRing F] [Algebra ℂ F]
variable {m k n : ℕ}

/-- **Proposition `prop:simple-membership`** (`NewMethod.tex:1281-1297`).

Source (verbatim), Robertz 2018, Proposition 3.31 ([Rob14], Prop. 2.2.50),
read from `~/project/papers/robertz-2018-formal-methods-systems-pdes-cours-
cirm-published.pdf`: "Let `S` be a simple differential system, defined over
`R`, with equations `p₁ = 0, …, p_s = 0`.  Moreover, let `E` be the
differential ideal of `R` which is generated by `p₁,…,p_s` and define the
product `q` of the initials and separants of all `p₁,…,p_s`.  Then
`E : q^∞` is a radical differential ideal.  Given `p ∈ R`, we have
`p ∈ E : q^∞` if and only if the pseudo-remainder of `p` modulo `p₁,…,p_s`
and their derivatives is zero."

NOTE (departure from source): "the pseudo-remainder … is zero" is rendered
here as `FullReduce p S.eqs = 0`; the paper itself makes the same
identification (`NewMethod.tex:1988-1992`, where Prop.
`simple-membership` is applied to Ritt's full remainder). -/
axiom simple_membership (I : AlgorithmInput R Rc F Ω m k n) (S : DiffSystem R)
    (hS : I.DifferentiallySimple S) (P : R) :
    (satBy (diffIdealBy I.d {y : R | y ∈ S.eqs}) (I.qOf S.eqs)).IsRadical ∧
    (P ∈ satBy (diffIdealBy I.d {y : R | y ∈ S.eqs}) (I.qOf S.eqs)
      ↔ I.FullReduce P S.eqs = 0)

/-- **Proposition `prop:decomposition-membership`**
(`NewMethod.tex:1300-1324`).

Source (verbatim), Robertz 2018, Proposition 3.32 ([Rob14], Prop. 2.2.72),
same file: "Let a (not necessarily simple) differential system `S` be given
by `p₁ = 0,…,p_s = 0`, `q₁ ≠ 0,…,q_t ≠ 0`, and let `S₁,…,S_r` be a Thomas
decomposition of `S` with respect to any ranking on `R`.  Moreover, let `E`
be the differential ideal of `R` generated by `p₁,…,p_s` and define the
product `q` of `q₁,…,q_t`.  For `i ∈ {1,…,r}`, let `E^{(i)}` be the
differential ideal of `R` generated by the equations in `S_i` and define
the product `q^{(i)}` of the initials and separants of all these equations
in `S_i`.  Then
`√(E : q^∞) = (E^{(1)} : (q^{(1)})^∞) ∩ … ∩ (E^{(r)} : (q^{(r)})^∞)`." -/
axiom decomposition_membership (I : AlgorithmInput R Rc F Ω m k n)
    (S : DiffSystem R) (Ss : List (DiffSystem R))
    (hSs : I.IsThomasDecomposition S Ss) :
    (satBy (diffIdealBy I.d {y : R | y ∈ S.eqs}) S.ineqs.prod).radical
      = ⨅ T ∈ Ss, satBy (diffIdealBy I.d {y : R | y ∈ T.eqs}) (I.qOf T.eqs)

end AlgorithmInput

end NewMethod
