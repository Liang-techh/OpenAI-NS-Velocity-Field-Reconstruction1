# Component long-annulus reshape — 2026-09-30

The pressure/width branch now carries the Section 9.30 long reshape from the common R=110 switch endpoint. It uses the same analytic preheat pressure, core, five moments and first axial tangents as the preceding bridge. Axis P0 is required explicitly; total endpoint pressure is never substituted for the axis datum.

For y=log(R/110), T=400A, B=log(Cstar Utheta110(1+Z^2)), the module evaluates

    log Utheta = log Utheta110 + y/10 - sigma(y/T) B
    a = 4/5 + 2 sigma_prime(y/T) B/T
    Uz = Uz110

Normalized backward quadrature computes angular, pressure and swirl increments. Axial increments use analytic primitives. All five inherited moment seeds and newly accumulated increments, their first Z derivatives, and separate raw axial/swirl quadratic integrals remain exposed. A rounded public sum alone must not be used for subsequent tiny-defect repair.

## Evidence

- Actual degree18, pressure9/width2 source at Z=.3, A=1e150, logC=5e151 reaches phases .5 and 1. T=4e152.
- Phase-zero field agreement is zero at recorded precision. Inherited five-moment components remain separate; the nonzero pressure-tail Uz_Z survives.
- Endpoint log-Utheta reference error is zero at recorded precision; Utheta_Z/Utheta differs from -2Z/(1+Z^2) by about 4.90e-409 in the leading atom.
- Independent resolved T=400 scalar adaptive quadrature: maximum normalized integral discrepancy 3.13e-26; fourth-order independent first-Z discrepancy 6.11e-17. This is a resolved fixture, not an actual-scale error bound.

Reproduce with Python from the repository root:

    python experiments/root_st073/lei_ren_part1_paper_pressure_width_long_reshape_fixture.py
    python experiments/root_st073/lei_ren_part1_paper_pressure_width_long_reshape_check.py

The actual check supports --resume using the locally generated ignored R100 cache. A fresh run reconstructs that cache; no cache is required from GitHub.

## Remaining work

Quadrature, finite pressure/width and source errors are unenclosed. Tail envelopes are conditional nominal bounds, not atomwise certified bounds. First axial tangents do not establish uniform C2 estimates. This branch is not yet installed as the complete global field and does not close terminal moments, finite radial energy, the stress cone, exact heat, flat remainder or temporal recursion.

Next implement the original reference-power extension, then axial restoration and functional moment repair. Keep the same axis pressure and retain inherited seeds separately from each interval's increments. Independent Cartesian validation belongs after the joined branch is installed; no global divergence certificate follows from local formulas alone.
