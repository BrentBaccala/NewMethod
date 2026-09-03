# The symmetry group of ansatz 25

*Goal.* Identify, as Lie groups, every transformation acting on the parameter
space of ansatz 25 (the Navier–Stokes system ansatz, `sec:NavierStokes`) —
separating those that change only the *description* of a solution from those
that change the solution — and read off which normalizations the 25.31–25.34
ladder spends, which it misses, and what each is worth.

*Status.* Complete for the `background = no` corner (25.2, 25.3 and the rungs
beneath them). Three results are actionable:

1. **Three group directions are unused** — a unipotent `𝔾ₐ` on the ODE's
   independent variable, the parabolic Navier–Stokes scaling, and the pressure
   gauge. Together they are worth 3 parameters and 12 constancy equations per
   rung, and — unlike the rotation normal form — each admits a *complete*,
   exclusion-free chart cover.
2. **One excluded locus is unrecorded.** The `w5 = 0` normalization at
   `norm ≥ 1` has its own isotropic exclusion, which
   `ansatz-library.sage:1018-1021` does not list.
3. **25.34's answer is a single orbit.** Its genuine component is
   3-dimensional, and those three free directions are *exactly* the three
   unused group directions — so there is no genuine modulus, and the group
   below is complete.

The trigger for this note was the observation that the ladder's design
(`ansatz-library.sage:914-926`) reasons entirely with *weights* — characters of
a torus. A weight table is the character lattice of a maximal torus and is by
construction blind to unipotent directions, which is precisely the shape of
what was missed.

---

## 1. The parameter space and the cost model

Ansatz 25 (`ansatz-library.sage:851-1031`) poses one unknown profile `Ψ` of one
linear inner variable

```
   s  =  s₁x + s₂y + s₃z + s₄t
```

obeying a second-order ODE whose nine jet monomials each carry a coefficient
*affine in s*,

```
   (a₀+a₁s)Ψ″² + (b₀+b₁s)Ψ″Ψ′ + (c₀+c₁s)Ψ″Ψ + (d₀+d₁s)Ψ′² + (e₀+e₁s)Ψ′Ψ
 + (f₀+f₁s)Ψ² + (g₀+g₁s)Ψ″ + (h₀+h₁s)Ψ′ + (i₀+i₁s)Ψ  =  0,
```

with the fields `u = u₅Ψ`, `v = v₅Ψ`, `w = w₅Ψ` (this is `background = no`) and
`p = p₁x + p₂y + p₃z + p₄t + p₅Ψ`.

The cost model is the reason any of this matters. `pconst` is built at
`ansatz-library.sage:1234` as one constancy equation per (parameter,
coordinate) pair, so

```
   #constancy equations  =  4 × #parameters,
```

while the ansatz contributes a fixed 14 equations at every rung. Each
dimension of symmetry that can be pinned is therefore worth exactly one
parameter and four equations, and task 473 established that the feasibility
cliff for this family sits between 9 and 19 parameters
(`~/project/reports/ansatz-25-ladder-edge.md`).

## 2. The group, in two layers

Write `P` for the parameter space. It carries an action of `G = Γ × H`, split
by whether the transformation moves the solution set or only its coordinates.

### 2.1 Γ — the gauge group

`Γ` is the kernel of `parameters ↦ solution`: a `Γ`-orbit is a set of parameter
points describing *the same* family of Navier–Stokes solutions. Quotienting by
`Γ` is lossless up to charts, since nothing about the solution moves.

| generator | action | type | used by |
|---|---|---|---|
| `κ` | multiply the whole ODE by a constant | `𝔾ₘ` | `g₀ = 1` (norm 2) |
| `λ` | `Ψ → λΨ` | `𝔾ₘ` | `v₅ = 1` (norm 2) |
| `ν` | `s → νs` | `𝔾ₘ` | `s₁ = 1` (norm 2) |
| **`τ`** | **`s → s + c`** | **`𝔾ₐ`** | **— unused —** |

`ν` and `τ` generate `Aff(1,ℂ) = 𝔾ₘ ⋉ 𝔾ₐ`, the affine group of the s-line,
which is the right group to expect: the ODE's coefficients are *affine
functions on that line*, so `Aff(1,ℂ)` is exactly what preserves the family.
`κ` and `λ` are central, and the torus acts on `τ` through the `ν` character,
so

```
   Γ  ≅  (𝔾ₘ)³ ⋉ 𝔾ₐ ,
```

a 4-dimensional connected solvable complex Lie group — three of whose
dimensions the ladder pins and one of which it does not see.

**Weights.** For the torus part, a parameter `p` has weight `(a,b,c)` when
`p ⟼ κᵃ λᵇ νᶜ · p`. For the coefficient of a jet monomial of Ψ-degree `dⱼ` and
total derivative order `wⱼ`:

```
   x₀ⱼ :  (1, −dⱼ,  wⱼ)          x₁ⱼ :  (1, −dⱼ,  wⱼ−1)
   s₁, s₄ : (0, 0, 1)            u₅, v₅, w₅, p₅ : (0, −1, 0)
   p₁, p₂, p₃, p₄ : (0, 0, 0)
```

The `g₀` row is `(1,−1,2)`, matching `ansatz-library.sage:918-919`. Two
consistency checks fall out. First, the parameters of *zero* gauge weight are
exactly `p₁ … p₄` — precisely the ones no scaling can normalize, which is why
`norm 3` had to justify `p₁ = p₂ = p₃ = 0` by a pseudo-division argument
rather than by symmetry, and why `p₄` needs a non-torus generator (§2.2).
Second, `τ` has no weight at all: unipotent elements are invisible to a
character table, which is the mechanism of the blind spot.

**τ in coordinates.** Under `s → s + c`, an affine coefficient transforms as

```
   (x₀ + x₁s)  ⟼  (x₀ + c·x₁) + x₁s ,     i.e.   x₀ ↦ x₀ + c·x₁,  x₁ ↦ x₁ ,
```

so the 9 coefficient pairs are 9 copies of the standard 2-dimensional
representation of `𝔾ₐ` (a single Jordan block), all sharing one parameter `c`.
Its invariants are the leading parts `x₁ⱼ` together with the brackets
`x₀ᵢx₁ⱼ − x₀ⱼx₁ᵢ` — the τ-quotient is cut out by those.

**Why τ is legitimate here, and only here.** Translating the solution in `x`
sends `u₅Ψ(s)` to `u₅Ψ(s+c)`, which is again of the ansatz form with shifted
ODE coefficients; the translation also adds a constant to `p`, which is
invisible because the momentum template (`ansatz-library.sage:1102-1106`) uses
only `p[x]`, `p[y]`, `p[z]`. With `background = yes` the fields carry a linear
part whose translate acquires a constant term the ansatz has no slot for, so
`τ` — like the Galilean boosts of §2.2 — belongs to the `background = no`
corner only.

### 2.2 H — the form-preserving Navier–Stokes symmetries

These are genuine point symmetries of incompressible Navier–Stokes that carry
the ansatz family to itself. They *permute* solutions, so pinning one computes
a slice whose result must be spread back over the orbit.

| generator | action | type | used by |
|---|---|---|---|
| `O(3,ℂ)` | rotate `(x,y,z)`, with `(u,v,w)`, `(s₁,s₂,s₃)`, `(p₁,p₂,p₃)` as vectors | dim 3 | `s₂=s₃=0`, `w₅=0` (norm 1) |
| **`δ`** | **`x→μx, t→μ²t, u→u/μ, p→p/μ²`** | **`𝔾ₘ`** | **— unused —** |
| **`π`** | **`p → p + f(t)`, i.e. the `p₄` direction** | **`𝔾ₐ`** | **— unused —** |

so `H° ≅ SO(3,ℂ) × 𝔾ₘ × 𝔾ₐ`, dimension 5. The parabolic scaling `δ` preserves
`ρ = μ_visc = 1`, so it survives the nondimensionalization of
`--pde navier-stokes-nd`.

Galilean boosts `x → x + ct, u → u + c` are *absent*: with `background = no`
the fields have no constant term to absorb the shift. The `background` flag is
therefore not merely a size knob — it is exactly the switch for boost-closure
of the family, and `background = yes` is the corner where the linear
strain-and-rotation flows live.

### 2.3 The ledger

`dim G = 4 + 5 = 9`. The ladder spends 6:

| rung | pins | generators spent | params | pconst |
|---|---|---|---|---|
| 25.3 | — | — | 28 | 112 |
| 25.31 | `s₂=s₃=0`, `w₅=0` | 3 of `O(3,ℂ)` | 25 | 100 |
| 25.32 | `s₁=1`, `v₅=1`, `g₀=1` | `ν`, `λ`, `κ` | 22 | 88 |
| 25.33 | `p₁=p₂=p₃=0` | *none* — not a symmetry quotient | 19 | 76 |
| 25.34 | drop the Ψ-quadratic block | *none* — a different family | 9 | 36 |

Unspent: `τ`, `δ`, `π`.

## 3. Two kinds of group need two kinds of slice

The reason the ladder's normalizations behave so differently is that `Γ` and
`H` are different kinds of group, and "clean slice" means a different thing for
each.

### 3.1 For the torus: unimodularity

Pinning three parameters kills `(𝔾ₘ)³` cleanly iff their `3×3` weight matrix
has determinant `±1`. The determinant of the ladder's choice is

```
   det [ (0,0,1); (0,−1,0); (1,−1,2) ]  =  1 ,
```

so the map `(𝔾ₘ)³ → (𝔾ₘ)³` given by that minor is an isomorphism and the slice
meets each orbit exactly once. Determinant `d` with `|d| > 1` would leave a
residual `μ_d` — the normalized family would be a `d`-fold cover of the true
quotient, and the decomposition would return `d` copies of every component.
Non-vanishing is not enough; unimodularity is the condition.

### 3.2 For SO(3,ℂ): orbit structure, via sl₂

`SO(3,ℂ) ≅ PGL(2,ℂ)`, and its standard 3-dimensional representation is the
**adjoint representation on `sl₂`**. Explicitly,

```
   (s₁,s₂,s₃)  ⟼  X = [ s₃          s₁ − i·s₂ ]        det X = −σ,
                       [ s₁ + i·s₂  −s₃       ]        σ = s₁²+s₂²+s₃² .
```

The orbits are therefore the `sl₂` conjugacy classes, of which there are three
kinds:

- `σ ≠ 0` — `X` regular semisimple (eigenvalues `±√(−det X)`); orbit = the
  smooth affine quadric `{σ = const}`, dimension 2, stabilizer a maximal torus.
- `σ = 0, X ≠ 0` — `X` regular **nilpotent**; orbit = the punctured null cone,
  dimension 2, stabilizer `𝔾ₐ`.
- `X = 0` — the origin.

`(s₁,0,0)` with `s₁ ≠ 0` is semisimple, so the rotation normal form reaches the
first stratum only: a nonzero isotropic wave covector cannot be rotated to
`(s₁,0,0)`, since that would force `s₁ = 0`. This is the exclusion recorded at
`ansatz-library.sage:1019-1021`, and task 471 found the corresponding component
(a 25.2 prime containing `s₁²+s₂²+s₃²` — the nilpotent orbit closure). Over `ℝ`
the nilpotent stratum is empty (`σ = 0 ⟹ s = 0`), which is exactly why the
exclusion is vacuous there and not over `ℚ̄`.

### 3.3 A second isotropic exclusion, unrecorded

The same argument applies a second time, to the *residual* rotation, and that
instance is not in the `excludes` list.

After `s₂ = s₃ = 0` the stabilizer of the wave covector is `O(2,ℂ)` acting on
the transverse plane, and the ladder spends it on `w₅ = 0`. But `SO(2,ℂ)`
preserves `v₅² + w₅²`, and on the isotropic line `(v₅,w₅) = (z, i·z)` the
rotated second component is

```
   w₅′  =  i·z·e^{−iθ} ,
```

which vanishes only at `z = 0`. So the orbit of a nonzero isotropic amplitude
never meets `{w₅ = 0}`, and the reflections in `O(2,ℂ)` do not help
(`w₅ → −w₅` maps `(z, i z)` to `(z, −i z)`, still isotropic). Rungs `norm ≥ 1`
therefore also give up

```
   { v₅² + w₅² = 0,  (v₅,w₅) ≠ (0,0) } ,
```

which `ansatz-library.sage:1018-1021` does not record. The locus is live rather
than vacuous: continuity forces `u₅ = 0` in the `s₁ ≠ 0` chart, so the
amplitude vector is purely transverse and lies in exactly the plane where the
isotropic lines sit. Over `ℝ` it is empty, for the same reason as the first
exclusion.

**Fix.** Add to the `norm in (1,2,3)` branch:

```python
excludes += ['v5^2 + w5^2 = 0, (v5,w5) != 0 '
             '(isotropic transverse amplitude; the residual SO(2,C) '
             'cannot rotate it to w5 = 0)']
```

## 4. The 25.34 answer is a single orbit

Task 473 found for 25.34 exactly one genuine prime,

```
   ⟨ u₅, p₅, i₀, i₁, h₀ + s₄, h₁ + g₁s₄ ⟩ ,
```

3-dimensional, with `s₄`, `p₄`, `g₁` free. Those three free directions are the
three unspent generators of §2.3, one each:

| free parameter | generator | action |
|---|---|---|
| `p₄` | `π` | `p₄` occurs in **no** equation at all — see below |
| `s₄` | `δ` | `s₄ ↦ μ·s₄`, with `s₁ = 1` preserved |
| `g₁` | `τ` | `g₁ ↦ g₁/(1 + c·g₁)` (after the `κ` compensation restoring `g₀ = 1`) |

`p₄` is the starkest: the momentum template uses `p[x]`, `p[y]`, `p[z]` and
never `p[t]`, and `p` never appears undifferentiated
(`ansatz-library.sage:1102-1106`). So `p₄` enters only its own defining field
equation, contributes 4 constancy equations, and the answer variety is a
product `V′ × 𝔸¹_{p₄}`.

Modulo `G`, then, the component is a **single point**. On it the ODE reads
`(1+g₁s)(Ψ″ − s₄Ψ′) = 0`, whose profile is

```
   v(x,t)  =  α + β·exp( s₄·(x + s₄t) ) ,
```

the Rayleigh–Stokes shear layer — with the whole three-parameter family it
appeared in being nothing but its symmetry orbit. That the counts match
exactly (3 free = 3 unspent) is the evidence that the group of §2 is complete:
no further symmetry is hiding in the family, and no free parameter in the
answer is a genuine modulus.

## 5. What is *not* a group: the `(1+g₁s)` multiplier

Task 473's compaction flagged multiplication of the ODE by the degree-1
polynomial `(1+g₁s)` as a fourth redundancy. It is real, and it is correctly
*not* a group action, so no normalization can pin it.

Two operators define the same solution family iff `L₂ = A·L₁` for a unit `A` —
what matters is the left ideal `D·L` in the Weyl algebra `ℂ[s]⟨∂⟩`, not `L`.
Restricted to a family with affine coefficients this relation is a groupoid
with jumping fibres, not a group action: over a constant-coefficient operator
the fibre `{(1+g₁s)·L : g₁ ∈ ℂ}` is 1-dimensional, over a generic affine
operator it is a point. `κ` is only its degree-0 part.

`τ` recovers half of it. The `τ`-orbits on the `g₁`-line are `{0}` and
`ℂ∖{0}`, so the shift moves `g₁` freely among nonzero values but fixes
`g₁ = 0`: the jump from a constant-coefficient operator to its `(1+g₁s)`
multiple is not realized by any group element. After pinning `g₁ ∈ {0,1}`
(§6) a residual 2:1 ambiguity remains — two parameter points, same
`D`-module — which is finite and harmless.

## 6. Three normalizations still on the table

All three are complete covers, not charts-with-exclusions: `τ` and `δ` each act
on their line with orbits `{0}` and `ℂ∖{0}`, so `{0,1}` is a full set of orbit
representatives, and `π` needs no split at all.

1. **`p₄ = 0`** (generator `π`). Lossless, no chart, no split — `p₄` appears in
   no equation. Costs nothing, saves one parameter and four equations.
2. **`g₁ ∈ {0,1}`** (generator `τ`). Two runs, covering everything. At
   `norm ≥ 1` where `g₀` is free the pin lands on a different coefficient, with
   the same effect.
3. **`s₄ ∈ {1,0}`** (generator `δ`). Two runs; `s₄ = 0` is the steady branch
   and is worth having separately in any case.

Effect, per rung:

| rung | now | with the three pins | runs |
|---|---|---|---|
| 25.33 | 19 params / 76 eqs | **16 / 64** | 4 |
| 25.34 | 9 / 36 | **6 / 24** | 4 |

The interesting target is 25.33 — the rung that actually diverged (MemoryError
at the 80 GB ulimit after 2:09:49, ~36 GB/h linear growth). Whether 4 runs at
16 parameters beat 1 run at 19 is a *hypothesis*, not a prediction: it rests on
cost growing fast enough in the parameter count for the trade to pay, which the
9-vs-19 cliff suggests but does not establish. Note also that task 473's own
reading of the cliff — that the Ψ-quadratic block rather than raw parameter
count is the blocker, since 25.34 is the rung that removes it — argues the
other way, and three parameters may simply not be enough to move a rung that
still carries the quadratic block. The two readings are cheap to separate:
run 25.33 with the pins, and if it still diverges the quadratic block is
confirmed as the driver.

A caution on the `w₅ = 0` interaction: §3.3's unrecorded exclusion is a defect
in the *bookkeeping*, not in the runs, but a `25.35` rung that adds these pins
should carry both isotropic entries in `excludes`, since it inherits `norm 1`.

## 7. What was checked, and how

Everything asserted here was verified rather than asserted from structure
alone; the symbolic checks were run under `python3` with sympy.

- **τ is a symmetry of the parameter variety.** Substituting `s → s+c` into the
  25.34 coefficients and renormalizing by `κ = 1/(1+g₁c)` gives
  `g₁ ↦ g₁/(1+cg₁)`, `h₀ ↦ (h₀+ch₁)/(1+cg₁)`, `h₁ ↦ h₁/(1+cg₁)`, and likewise
  for `i₀, i₁`. Substituting the genuine component `i₀=i₁=0, h₀=−s₄,
  h₁=−g₁s₄` into that map returns `h₀′ + s₄ = 0` and `h₁′ + g₁′s₄ = 0` — the
  component is τ-stable. Solving `g₁′ = 1` gives `c = (g₁−1)/g₁`, valid iff
  `g₁ ≠ 0`, confirming the orbit structure of §5.
- **The second isotropic exclusion.** `SO(2,ℂ)` preserves `v₅²+w₅²`; solving
  `w₅′ = 0` for `θ` in general returns `2·atan((v₅ ± √(v₅²+w₅²))/w₅)`, which
  requires `v₅²+w₅² ≠ 0`; on `(z, i·z)` the rotated component is `i·z·e^{−iθ}`,
  vanishing only at `z = 0`.
- **The profile solves the PDE, and δ acts as claimed.**
  `v = α + β·exp(s₄x + s₄²t)` gives `v_t − v_xx = 0` identically; under
  `x→μx, t→μ²t, u→u/μ` it becomes `(α + β·exp(μs₄(μs₄t + x)))/μ`, i.e. the same
  form with `s₄ ↦ μs₄` and the amplitudes absorbed into `Ψ`'s own integration
  constants — so `s₁ = 1` is preserved and no `λ` compensation is needed.
- **`p₄` occurs in no equation.** Read directly off
  `ansatz-library.sage:1102-1106`.
- **Parameter counts.** 25.3 = 4 (`s`) + 8 (fields) + 16 (ODE) = 28;
  25.34 = 1 + 3 + 5 = 9. Both match the code and the run reports.

## 8. Gaps

- **`background = yes` is not covered.** Boosts re-enter, `τ` drops out, and
  the group is a different (larger, non-reductive) object. Ansätze 25 and 25.1
  are not analyzed here.
- **The payoff of §6 is untested.** No run has been made with the three pins;
  the parameter counts are exact but the cost claim is a hypothesis.
- **Whether the ladder should quotient `H` at all** is a modelling question
  this note does not settle. `Γ` is free to quotient; every `H` pin trades
  parameters for an obligation to spread the answer back over an orbit, which
  is bookkeeping the paper must then carry. The `O(3,ℂ)` pin already costs two
  isotropic exclusions; `δ` and `π` cost none, which is an argument for
  taking those two first.
- **The `μ_d` question for the proposed pins.** Unimodularity was checked for
  the ladder's existing torus pins. `δ` acts on `s₄` with weight 1 and `π` on
  `p₄` by translation, so neither can leave a residual finite stabilizer; `τ`'s
  residual is the 2:1 of §5, which is understood. No analogue of the
  determinant check is needed, but a `25.35` rung should say so explicitly.

---

*Written by Claude (Opus 5) in an interactive session of 2 September 2026,
following the ansatz-25 normalization ladder of tasks 472/473; the symbolic
verifications of §7 were executed live during the session. Reviewed and edited
by the author — see `AI-CONTRIBUTIONS.md`.*

/claude-opus-5
