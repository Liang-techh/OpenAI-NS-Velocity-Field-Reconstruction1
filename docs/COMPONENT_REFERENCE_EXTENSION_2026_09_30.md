# Component reference continuation â€” 2026-09-30

This interval continues the Section9.30 reshape endpoint Rsh=110 exp(400A) to Rz=exp(-8)Rref, where axial restoration begins. The reference radius uses the common source relation logRref=log110+10(logCstar+logPstar).

For x=log(R/Rsh), continue Utheta=Utheta_sh exp(x/10), F=Fsh exp(-2x/5), and Uz=Uz_sh. Continue their actual inherited first Z derivatives rather than resetting a derivative to a desired target. The reference endpoint shape was checked separately in the preceding receipt.

The five moments and raw quadratic integrals use exact power-law primitives on this interval. Original inner seeds, reshape increments and new reference increments are exposed as separate parts, including their first Z derivatives. Axis pressure P0 is preserved explicitly and total pressure receives only the new pressure-moment increment. No pressure gauge repair or terminal moment substitution is performed.

## Independent resolved evidence

A T=400, pressure3/width2 fixture independently integrates the physical moment derivatives at x=1 and x=5. Maximum scaled integral discrepancy is 3.62e-100; fourth-order first-Z stencil discrepancy is 1.02e-15; x=0 continuity discrepancy is zero at recorded precision. Separate part receipts agree to 2.94e-100. These checks establish the implemented finite-ring primitives on the resolved fixture, not actual-scale source or global error bounds.

Reproduce:

    python experiments/root_st073/lei_ren_part1_paper_pressure_width_reference_extension_fixture.py
    python experiments/root_st073/lei_ren_part1_paper_pressure_width_reference_extension_check.py

The actual check supports --resume only for the locally generated ignored R100 cache. A fresh run reconstructs the common source.

## Actual-source evidence

Degree18 pressure9/width2 common-source Z=.3 replay reaches one logarithmic unit beyond Rsh and the actual Rz axial-restoration entry. Recorded x=0 field errors are zero. Original seeds and reshape increments remain separately unchanged; nonzero Uz_Z pressure-tail atoms survive; P0 remains unchanged. These are local propagation checks, not functional terminal closure.

## Scope

The interval primitives introduce no new quadrature. Source, reshape quadrature and finite pressure/width-ring errors are inherited and remain unenclosed. Public summed quantities may lose small contributions within an atom; functional repairs must consume separate parts. This module does not establish global matching, finite energy, uniform C2, the stress cone, exact heat or temporal recursion.

## Next

Implement Section9.38 on Rz <= R <= e Rz: Uz=v1+(4Z-v1)sigma(log(R/Rz)). Continue Uz=4Z afterward to Rh=exp(-5)Rref, retaining all previous moment parts and new axial-restoration integrals. The correction interval [Rm,2Rm], Rm=exp(-6)Rref, must restore the five reference moment functions without changing axis pressure afterward.
