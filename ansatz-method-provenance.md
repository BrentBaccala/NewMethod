# Provenance of the "Ansatz Method": what NewMethod inherits and where it departs

**Date:** 2026-07-18
**Context:** Companion note to `NewMethod.tex`, written to settle the paper's
related-work positioning around the word **"ansatz method"** — which appears as
a *keyword* in the closest prior art (Chaharbashloo et al. 2013) and is often
waved off as a centuries-old idea. It is not: in this literature it is a
specific, nameable term of art with a datable coiner and a precise technical
origin. This note traces that provenance and states, layer by layer, what
NewMethod inherits and where it departs. All citations are verified (CrossRef /
primary PDF); companion to the prior-art note
`~/project/papers/chaharbashloo-basiri-rahmany-zarrinkamar-2013-ansatz-qm-grobner.md`
and the litsearch report
`~/project/reports/newmethod-parametric-differential-elimination-litsearch.md`.

---

## TL;DR

- **"Ansatz" the generic word** (German: an assumed trial form) is untraceable
  to one origin and is not what the keyword means here.
- **"The ansatz method" as a named QM technique** was branded by **Shi-Hai Dong
  (2002)** — write ψ = (polynomial) × exp(f), substitute into the *radial ODE*,
  equate powers. Chaharbashloo (2013) then solves the resulting polynomial
  system with a **commutative** Gröbner basis.
- **Why the ansatz closes** into a finite algebraic system is the
  **quasi-exactly-solvable (QES)** "hidden Lie algebra" (Turbiner 1988;
  classified and lifted to two variables by González-López–Kamran–Olver).
- **Everything in this lineage is either a single separated ODE, or a PDE that
  is reachable only under a strong precondition** (separability, or a hidden
  finite-dim Lie algebra of differential operators).
- **NewMethod's contribution is orthogonal to all of it:** it mechanizes the
  *differential* step by differential-algebraic elimination on the
  **un-separated** system, needs **no** hidden algebra and **no** separability,
  and carries a **completeness certificate** — which is exactly why it reaches
  non-separable, no-hidden-symmetry problems (helium) where the QES machinery
  cannot even start.

---

## The four layers of provenance

### Layer 0 — "Ansatz," the generic term
German *Ansatz* = an assumed functional form for a solution, substituted and
then solved for its free parameters. Ubiquitous in math/physics, no single
origin. The only ansätze with a canonical first paper are the **Bethe ansatz**
(H. Bethe, *Z. Physik* 71:205, 1931), the **coupled-cluster exponential
ansatz** (Coester 1958 / Čížek 1966), and the **Ritz** variational trial
function (Crelle 135:1, 1909). Chaharbashloo's "ansatz method" is **none** of
these.

### Layer 1 — The named "ansatz method" for the Schrödinger equation: Dong (2002)
The specific technique — ψ = (node polynomial) × exp(y(r)) with y a solution of
a **Riccati** equation, substituted into the D-dimensional/radial Schrödinger
equation, powers equated to close a finite polynomial system in the
parameters — was made a **titled, named method by Shi-Hai Dong**:
- S.-H. Dong, *"On the solutions of the Schrödinger equation with some
  anharmonic potentials: wave function ansatz,"* **Phys. Scr. 65:289–295 (2002)**,
  DOI 10.1238/Physica.Regular.065a00289.
- S.-H. Dong, *"The Ansatz Method for Analyzing Schrödinger's Equation with
  Three Anharmonic Potentials in D Dimensions,"* **Found. Phys. Lett.
  15:385–395 (2002)**, DOI 10.1023/A:1021220712636.
- Consolidated as the chapter *"Wavefunction Ansatz Method"* (Ch. 8, *Wave
  Equations in Higher Dimensions*, Springer 2011,
  DOI 10.1007/978-94-007-1917-0_8).

Propagated as a keyword by the **Turkish school** (Ikhdair & Sever, *"…wave
function ansatz,"* Cent. Eur. J. Phys. 6(3):697, 2008,
DOI 10.2478/s11534-008-0060-y) and the **Iranian school**
(Hassanabadi–Rajabi–Zarrinkamar, *Mod. Phys. Lett. A* 27:1250057, 2012 — the
immediate source Chaharbashloo lifts its ansatz from). All Dong 2002 papers are
paywalled with no OA copy or arXiv preprint (verified 2026-07-18).

### Layer 2 — The QES trial form that Dong named: Flessas (1979) and the QES program
The "polynomial × exponential, equate powers" idea predates Dong's naming:
- **G. P. Flessas**, *"Exact solutions for a doubly anharmonic oscillator,"*
  **Phys. Lett. A 72(4–5):289–290 (1979)**, DOI 10.1016/0375-9601(79)90471-7;
  part II, **81(1):17–18 (1981)**, DOI 10.1016/0375-9601(81)90292-9; and the
  general-potential statement **"Exact solutions for anharmonic oscillators,"
  J. Phys. A 14(6):L209–L211 (1981)**, DOI 10.1088/0305-4470/14/6/001 (the
  most-cited of the cluster). CrossRef-verified. All paywalled, pre-arXiv.
  *Caution:* "Exact quantum-mechanical solutions for anharmonic oscillators,"
  Phys. Lett. A 81:116 (1981), is by **E. Magyari, not Flessas.**
- This sits inside the **quasi-exactly-solvable (QES)** program — the class of
  potentials for which a *finite* number of eigenstates is algebraically
  computable while the rest are not. Founding work: **A. V. Turbiner**,
  *"Quasi-exactly-solvable problems and sl(2) algebra,"* Comm. Math. Phys.
  118:467–474 (1988); **C. M. Bender & G. V. Dunne**, J. Math. Phys. 37:6–11
  (1996) (QES orthogonal polynomials); consolidated in **A. G. Ushveridze**,
  *Quasi-Exactly Solvable Models in Quantum Mechanics*, Institute of Physics
  Publishing, Bristol, 1994, ISBN 0-7503-0266-5.

### Layer 3 — Why the ansatz closes: the hidden Lie algebra (Turbiner; GKO)
The QES "hidden algebra" is the structural reason the poly×exp ansatz reduces
to a *finite* polynomial system. After a gauge rotation ψ = e^{−g}φ and a change
of variable, the Schrödinger operator becomes a **quadratic element of the
enveloping algebra of a finite-dimensional Lie algebra of first-order
differential operators**. In 1D that algebra is **sl(2,ℝ)**:
```
J⁺ = x²∂_x − n·x,    J⁰ = x∂_x − n/2,    J⁻ = ∂_x,
```
which preserves the (n+1)-dimensional module P_n = span{1, x, …, xⁿ}. The
operator therefore maps P_n → P_n; restricting to P_n yields an
(n+1)×(n+1) matrix whose eigenvalues are precisely the algebraically accessible
energies — the remaining infinite spectrum is invisible (hence *quasi*-exact).
"Hidden" because the algebra appears only after the gauge transform and **does
not commute with H** (it is a spectrum-generating structure, not a symmetry).

This is continuous with **Peter Olver's** work on Lie algebras of differential
operators (distinct from his classical Lie-*point-symmetry* theory): with
González-López and Kamran he **classified** the finite-dimensional Lie algebras
of first-order differential operators possessing a finite-dimensional invariant
module, and **lifted the construction to two variables** — the genuine
non-separated **PDE** QES branch:
- A. González-López, N. Kamran, P. J. Olver, *"Quasi-exactly solvable Lie
  algebras of first order differential operators in two complex variables,"*
  **J. Phys. A 24:3995–4008 (1991)**.
- A. González-López, N. Kamran, P. J. Olver, *"New quasi-exactly solvable
  Hamiltonians in two dimensions,"* **Comm. Math. Phys. 159:503–537 (1994)**,
  DOI 10.1007/BF02099982. (Note: this DOI is CMP **159** (1994), not CMP 153
  (1993) — the latter is GKO's distinct "Normalizability of one-dimensional QES
  Schrödinger operators.") Both PDFs are free on Olver's site and archived in
  `~/project/papers/`.

---

## The pipeline, and where NewMethod stands

Every method in this lineage has the shape:

> **(one radial ODE, or a separable/hidden-algebra-structured PDE)**
> → *hand* differential work (substitute, differentiate, equate powers)
> → **finite polynomial system in scalar parameters**
> → solve it (by hand: Dong, Flessas; by **commutative** Gröbner basis:
>   Chaharbashloo 2013).

Two structural commitments are shared by all of them and *not* by NewMethod:

1. **The differential step is manual.** The substitution/elimination that turns
   the differential equation into a polynomial system is done on paper. The
   computer (when used at all, as in Chaharbashloo) only touches the
   *commutative* residue.
2. **The problem must be pre-reduced.** Either it is a single separated ODE
   (Dong/Flessas/Chaharbashloo), or it is a PDE reachable only via a **strong
   precondition** — separability (Ushveridze's inverse method) or a **hidden
   finite-dim Lie algebra** of differential operators (González-López–Kamran–
   Olver). Generic potentials satisfy neither.

**NewMethod inverts both.** It mechanizes the differential step itself by
**differential-algebraic elimination** (Ritt reduction / Rosenfeld–Gröbner over
a differential ring tower) applied to the **un-separated** differential ideal;
it requires **no** separability and **no** hidden Lie algebra; and its
decomposition carries a **completeness certificate** (Theorem 3 of
`NewMethod.tex`) — a guarantee that *every* solution of the ansatz form is
found, which none of the hand-substitution methods can even state. Chaharbashloo
is the closest prior art precisely because it brings computer algebra into the
picture at all — but a *commutative* Gröbner basis on the post-ODE polynomial
system, not a *differential* one on the equation.

The one-line differentiator, sharpened against the strongest version of the
prior art (the QES-PDE branch):

> Prior algebraic methods for QM eigenvalue problems either work on a single
> separated ODE, or reach PDEs only under a separability or hidden-Lie-algebra
> precondition; **NewMethod performs general differential elimination on the
> un-separated system, with a completeness guarantee, and therefore applies
> where those preconditions fail** — as in the non-separable, no-hidden-sl(2)
> helium problem.

---

## Recommended citations for `NewMethod.tex`
Credit the tradition precisely rather than diffusely:
- **Dong 2002** (Phys. Scr. 65:289 and/or Found. Phys. Lett. 15:385) — the
  term-of-art source for "the ansatz method."
- **Flessas 1979/1981** and **Flessas, J. Phys. A 14:L209 (1981)** — the
  exact poly×exp trial-form origin.
- **Turbiner 1988; Bender–Dunne 1996; Ushveridze 1994** — the QES program.
- **González-López–Kamran–Olver 1991/1994** — the hidden-Lie-algebra
  classification and its multivariate (PDE) extension; cite when drawing the
  precondition-free / completeness line.
- **Chaharbashloo et al. 2013** (Z. Naturforsch. A 68:646) — the closest prior
  art (ansatz + commutative Gröbner), the paper to engage head-on.

*Provenance research and this note by Claude (Opus 4.8), reviewed and edited by
the author; see `AI-CONTRIBUTIONS.md`.*
