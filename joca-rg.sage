# -*- mode: python -*-
#
# Sage script to do the computations for my paper in the Journal of Computational Algebra.
#
# This is the variant of joca.sage that adds the general algorithm's
# Rosenfeld-Gröbner regularization of the ansatz (see the block below and the
# README).  joca.sage itself is the projection-only version whose output
# matches the paper's printed decomposition.
#
# Author: Brent Baccala
# Date: December 12, 2025
#
# Tested on Ubuntu 24 with Sage 10.6
#
# François Boulier's Differential Algebra package is required to run this script.
# To install it, run this command from the sage prompt:
#
# %pip install DifferentialAlgebra
#
# DifferentialAlgebra (sympy based) is used to perform the differential algebra reduction.
# Sage is used to compute a primary decomposition of the resulting system of equations.
# Some fiddling is required to juggle back and forth between the two.

import sys

import sympy

# Command-line options.  Run e.g. `sage joca-rg.sage --print-rg-system` to also
# print the regular differential chain that RosenfeldGroebner returns (the
# *regularized* ansatz).  By default we only use RG as the coherence check and
# project against the raw ansatz -- see the note at the RG call below -- so the
# regularized system is not normally shown.
print_rg_system = '--print-rg-system' in sys.argv

# `--rg-verbose` traces the Rosenfeld-Groebner computation itself: it passes
# verbose=True to RosenfeldGroebner, which prints its step-by-step operation
# (rounds, quadruples handled, reg_characteristic calls, new equations,
# reductions to zero, splits) to standard output.  This requires the patched
# DifferentialAlgebra build that wires BLAD's verbose splitting tree through to
# a Python-level `verbose` argument.
rg_verbose = '--rg-verbose' in sys.argv

# `--rg-dot` dumps the Rosenfeld-Groebner splitting tree in graphviz/dot syntax
# (a `digraph`) to standard output, by passing dot=True to RosenfeldGroebner.
# Pipe the output to `dot -Tpng -o tree.png` to render the splitting tree (the
# branch structure: factor/initial/separant splits, critical-pair edges,
# reductions to zero, reg_characteristic terminals).  Like --rg-verbose this
# requires the patched DifferentialAlgebra build that wires BLAD's splitting
# tree through to a Python-level `dot` argument; the two are independent and
# may be combined.
rg_dot = '--rg-dot' in sys.argv

try:
    import DifferentialAlgebra
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(ex.msg + "\nInstall it with '%pip install DifferentialAlgebra'")

# Claude Sonnet 4's solution to printing "DPsi" as "\Psi'" in LaTeX

def patch_latex_varify():
    """
    Patch sage.misc.latex.latex_varify to handle 'DPsi' specially.
    """
    from sage.misc.latex import latex_varify
    import sage.misc.latex

    # Save reference to the original function
    original_latex_varify = latex_varify

    # Define the custom function
    def custom_latex_varify(a, is_fname=False):
        if a == "DPsi":
            return r"\Psi'"
        else:
            return original_latex_varify(a, is_fname=is_fname)

    # Replace the function in the module
    sage.misc.latex.latex_varify = custom_latex_varify

patch_latex_varify()

# Declare our independent variables
x,y,z = sympy.var('x,y,z')

# Declare our constants
E = sympy.var('E')
v1,v2,v3,v4 = sympy.var('v1,v2,v3,v4')
a0,a1,b0,b1,c0,c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E,v1,v2,v3,v4,a0,a1,b0,b1,c0,c1]

# Declare our dependent variables
#
# `indexedbase` is DifferentialAlgebra's suggested way of declaring dependent variables
# if we want to use jet notation (i.e, v[x] for dv/dx) to write their derivatives.

Psi,DPsi,DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')

# Create the DifferentialRing and set the ranking on the variables

DiffRing = DifferentialAlgebra.DifferentialRing (derivations = [x,y,z],
                                                 blocks = [[DDPsi,DPsi,Psi,v,r], constants],
                                                 parameters = constants,
                                                 notation = 'jet')

# Define the PDE we're trying to solve
# sympy can't handle a Sage Integer, so use casts to make these Python integers

PDE = -int(1)/int(2)*(Psi[x,x] + Psi[y,y] + Psi[z,z])*r - Psi - E*r*Psi

print("PDE:", PDE)

# Define the ansatz, the parameterized function space in which we're looking for solutions

ansatz = [Psi[x] - DPsi * v[x],
          Psi[y] - DPsi * v[y],
          Psi[z] - DPsi * v[z],
          DPsi[x] - DDPsi * v[x],
          DPsi[y] - DDPsi * v[y],
          DPsi[z] - DDPsi * v[z],
          (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
          v - (v1*x + v2*y + v3*z + v4*r),
          r**2 - x**2 - y**2 - z**2]

# DifferentialAlgebra can't handle the parenthesized expressions directly, so expand them

ansatz = list(map(sympy.expand, ansatz))

print("\nAnsatz:", *ansatz, sep='\n')

# Run Rosenfeld-Gröbner to DISCHARGE HYPOTHESIS (3): confirm the ansatz forms a
# single coherent, squarefree regular differential system.  The constants must
# be moved into the coefficient field Q(E,v1,...,c1) for RosenfeldGroebner to
# terminate; left as ring variables it does not finish (see rg_basefield.py and
# the README).
#
# RG is used here ONLY as the coherence check, NOT to rewrite the reduction set.
# It is *not* a no-op: it squares the ODE's initial -- the regularized ODE is
# (a0+a1*v)*(original ODE) -- and the regular chain it returns represents the
# SATURATED ideal [C] : H_C^inf, with H_C = (a0+a1*v) its initial/separant.
# Reducing the PDE against that chain's bare .equations() does NOT saturate, so
# it re-admits the bad locus H_C = 0  <=>  a0 = a1 = 0 as a SPURIOUS prime
# (a0,a1); and the two genuine strata that live inside a0=a1=0 (where the
# 2nd-order ODE degenerates to 1st order) fall in the chain's bad locus, so
# saturating to remove the artifact loses them too.  Net: reducing against the
# regularized chain gives a wrong decomposition (4 primes incl. the spurious
# (a0,a1); saturating leaves only the 3 generic strata).  We therefore project
# against the RAW (coherent) ansatz below -- exactly joca.sage's reduction --
# which keeps the a0=a1=0 information and is the faithful decomposition (5
# primes).  See the README for the full diagnosis.

F = DifferentialAlgebra.BaseFieldExtension(generators=constants, ring=DiffRing)
components = DiffRing.RosenfeldGroebner(ansatz, basefield=F, verbose=rg_verbose, dot=rg_dot)
assert len(components) == 1, f"expected a single regular component, got {len(components)}"
print("\nHypothesis (3) discharged: the ansatz is a single coherent regular system.")

# Optionally print the regularized system RG returns -- the regular differential
# chain's equations.  This is NOT what we reduce against (see the note above and
# the next line): the chain represents the saturated ideal [C] : H_C^inf, with
# the ODE's initial squared, so reducing against it would drop the a0=a1=0
# strata.  We print it only to show what the general algorithm's hypothesis-(3)
# step produces.
if print_rg_system:
    rchain = components[0]
    rg_eqns = rchain.equations()
    print("\nRegularized system from Rosenfeld-Groebner (the regular chain):",
          *rg_eqns, sep='\n')

    # Try to factor each equation of the regularized system.  sympy.factor()
    # leaves the irreducible ones expanded; where the chain carries a repeated
    # initial (e.g. the squared ODE initial a0+a1*v) the factored form makes it
    # visible.
    print("\nFactored regularized system:",
          *[sympy.factor(eq) for eq in rg_eqns], sep='\n')

    # The component's inequalities.  The sympy binding does not wrap BLAD's
    # `Inequations` accessor, but that accessor is *defined* as the initials and
    # separants of the regular chain (Boulier): the chain represents the
    # saturated ideal [C] : H_C^inf, with H_C the multiplicative family they
    # generate -- a regular differential chain has no inequalities beyond these.
    # We reproduce `Inequations` faithfully: reduce each initial/separant modulo
    # the chain (normal_form -- a no-op here, since the ranking already keeps
    # x^2+y^2+z^2 rather than r^2), make it primitive, drop the numeric ones,
    # and dedupe.
    rg_inequalities = []
    seen = set()
    for h in list(rchain.initial()) + list(rchain.separant()):
        h = sympy.expand(rchain.normal_form(h))
        if h.is_number:
            continue
        h = h.as_content_primitive()[1]          # strip the numerical content
        if h in seen:
            continue
        seen.add(h)
        rg_inequalities.append(h)
    print("\nInequalities (initials and separants of the chain = BLAD 'Inequations'), must be nonzero:",
          *[f"{h} != 0" for h in rg_inequalities], sep='\n')

    # ... and factored, which exposes the bad-locus structure (e.g. the product
    # (x^2+y^2+z^2) * (ODE initial) and the ODE initial itself, whose vanishing
    # is the a0=a1=0 locus discussed above).
    print("\nInequalities, factored:",
          *[f"{sympy.factor(h)} != 0" for h in rg_inequalities], sep='\n')

# NB: deliberately do NOT overwrite `ansatz` with components[0].equations();
# we reduce against the raw (coherent) ansatz -- see the note above.

# Reduce the PDE modulo the (regularized) ansatz using Ritt's reduction algorithm

h,r = DiffRing.differential_prem(PDE, ansatz)

print("\nRemainder:", r)

# Convert the remainder to Sage

PolyRing = PolynomialRing(QQ, names=[str(indet) for indet in DiffRing.indets(selection='all')])
PolyRing_constants = list(map(PolyRing, constants))
PolyRing_r = PolyRing(r)

# Given an equation and a list of constants, factor each term into constant and non-constant factors,
# then group together all terms with identical non-constant factors and return the resulting
# list of equations (which will only involve constants).

def build_system_of_equations(eqn, constants):
    ring = eqn.parent()
    system_of_like_terms = dict()
    non_constant_sub = tuple(1 if ring.gen(n) in constants else ring.gen(n) for n in range(ring.ngens()))
    for coeff, monomial in eqn:
        non_constant_part = monomial(non_constant_sub)
        constant_part = coeff * monomial // non_constant_part
        if (non_constant_part) in system_of_like_terms:
            system_of_like_terms[non_constant_part] += constant_part
        else:
            system_of_like_terms[non_constant_part] = constant_part
    return tuple(set(system_of_like_terms.values()))

eqns = build_system_of_equations(PolyRing_r, PolyRing_constants)

print("\nSystem of equations:", *eqns, sep='\n')

# Build a polynomial ideal from the system of equations and construct its prime decomposition

I = ideal(eqns)
prime_decomposition = I.minimal_associated_primes()

# Sort this result (so it prints in the same order as in the JOCA paper) and print it

prime_decomposition.sort(key=lambda x:str(x))
print("\nMinimal associated prime ideals:", *prime_decomposition, sep='\n')
