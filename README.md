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

### `joca-rg.sage` — the general algorithm's regularize-then-reduce route (≈17 s)
Runs the general algorithm literally, *using* the Rosenfeld–Gröbner result for
the calculation. The constants are moved into the coefficient field
`ℚ(E,v₁,…,c₁)` via `BaseFieldExtension` (**required** for `RosenfeldGroebner`
to terminate — see the cliff below); RG returns a **single** regular component,
discharging hypothesis (3) (the paper's "B = ∅, one component" claim); and the
PDE is then Ritt-reduced **against the regular chain RG returned** (its
`.equations()` are the reduction set), projected, and prime-decomposed. The
change of reduction set is the entire difference from joca.sage. Output
(`joca-rg.out`): the chain, the Ritt denominator, and **four** minimal primes.

**What the output shows — regularize-then-reduce is not faithful.** RG is *not*
a no-op: it squares the ODE's initial — the chain's ODE element is
`(a₀+a₁v)·(original ODE)`, rationalized — and the chain represents the
*saturated* ideal `[C] : H_C^∞`, where `H_C = (a₀+a₁v)` is the
initial/separant. Reducing the PDE against the chain's bare `.equations()` does
**not** itself saturate, so the projection re-admits the bad locus
`H_C = 0 ⇔ a₀=a₁=0` as a **spurious** prime `(a₀,a₁)` — the PDE is *not*
actually redundant there (witness: the genuine redundancy ideal contains
`v₄·b₁`, nonzero on generic `a₀=a₁=0`) — while the two genuine strata *inside*
`a₀=a₁=0` (the paper's \eqref{ideal:3}, \eqref{ideal:4}, where the 2nd-order
ODE degenerates to 1st order) are absorbed into that plane and lost. The four
primes are: the classical branch, the `E = 0` parabolic branch, the `v ≡ 0`
component, and the spurious `(a₀,a₁)`; saturating away the artifact would also
delete the two genuine strata, leaving 3. The mechanism is visible in the
printed Ritt denominator, `8(x²+y²+z²)³·N` with `N` the rationalized norm of
the squared initial `(a₀+a₁v)²`: it vanishes identically on `a₀=a₁=0`, so this
route cannot resolve strata inside the locus where an initial vanishes — the
paper's bad-locus `B` discussion made concrete. Only the raw projection
(joca.sage) keeps the `a₀=a₁=0` information and yields the faithful **five**
primes.

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

For the full record of AI contributions to the paper (kept for the journal's
generative-AI declaration), see **`AI-CONTRIBUTIONS.md`**.
