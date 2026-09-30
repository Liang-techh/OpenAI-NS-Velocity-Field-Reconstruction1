# Pressure/width components of flat reshape defects

This step lifts the flat-shape scalar kernel into the same rectangular finite pressure/width algebra as the core and annular connection. It uses the inherited B and B_Z component atoms, without first replacing them by a rounded total.

For n>=1 the coefficient sensitivity is

    I^(n)(B0) = integral_0^T exp(-k ell) (m q)^n exp(m B0 q) d ell,
    q = sigma(ell/T).

Each sensitivity has its own small-q source mode near (2 n T^2/k)^(1/3). Reusing the first-order mode would discard higher-order contributions.

For Delta=B-B0 and rectangular orders P,W, Delta^(P+W+1)=0 in the declared finite ring. Compose

    I(B) = sum_n=0..P+W I^(n)(B0) Delta^n/n!
    I_Z = [sum_n=0..P+W I^(n+1)(B0) Delta^n/n!] B_Z.

The adapter retains each derivative-order term separately before its public component sum. High-order contributions can be much smaller than another contribution to the same atom; public rounding does not remove the separate term receipt. Pressure/width orders here are parameter expansions, not the temporal coefficient recursion required later.

## Evidence

The actual serialized pressure9/width2 B/B_Z atoms pass all three kernels. Each public value and tangent retains all30 rectangular atoms. Terms n=0..11 remain separate, including a nonzero n=11 contribution to atom(9,2); derivatives through n=12 are used for the tangent. This is propagation of the serialized common-source atoms, not a newly reconstructed source or a global error bound.

Independent resolved P2/W1 coefficientwise integration of the full jet integrand gives maximum relative value discrepancy3.73e-60 and first-Z discrepancy3.39e-61. Independent scalar derivatives n=1,2,3 at T=400 differ by at most3.83e-59. Extreme-scale checks confirm order-dependent modes; no tail or quadrature enclosure is claimed.

Reproduce:

    python experiments/root_st073/lei_ren_part1_paper_flat_shape_derivatives_fixture.py
    python experiments/root_st073/lei_ren_part1_paper_flat_shape_component_fixture.py
    python experiments/root_st073/lei_ren_part1_paper_flat_shape_component_check.py

## Scope

The three angular reshape sources are component building blocks, not all five functional defects. Inner/axial/reference/restoration centered sources and pressure compatibility must still be combined. Finite quadrature windows, source errors, finite-ring remainders and functional closure are unenclosed. No stress cone, global finite energy, heat exterior or temporal recursion is certified by this lift.
