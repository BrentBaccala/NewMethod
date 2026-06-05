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
