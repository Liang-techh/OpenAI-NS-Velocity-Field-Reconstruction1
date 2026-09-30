# Analytic positive-epsilon continuation after the collar

The derivative-aware prescribed exit now continues beyond the collar toward R=100 without freezing its velocities. The Section 9.23 auxiliary field is frozen by definition after its transition; its pressure and cumulative moments still grow with radius. This auxiliary branch is distinct from the prescribed field, whose positive-epsilon shear must continue changing its velocities and moments.

With y=log(R/Rb), the auxiliary driver D has powers exp(y), 1 and exp(-y). J=sqrt(R/2)*I_z has powers exp(2y), exp(y) and 1. Their coefficients are recovered by constant linear algebra at three radii. An unused fourth radius checks this representation. The g=log(F/F_start) equation integrates directly. Because epsilon=W has strictly positive width degree, exp(g) terminates algebraically at the retained width order. Uz and all five moments are then integrated as exponential polynomials, including polynomial factors in y. No continuation quadrature or additional RK discretization is used.

The continuation retains endpoint pressure-tail and width atoms and their automatic Z tangents. Fixed-R queries use the component-valued y=log(R/Rb), preserving the tiny collar endpoint displacement. Pressure, radial velocity and stress use the resulting common field and moment data. Derivatives of the actual analytic Mz primitive supply a local similarity divergence check, separately from substituting moment ODE identities. The inherited collar RK and analytic core/data errors remain unbounded.

The actual-source run completed at R=1,10,100. Both F and Uz retain nonzero first-width changes at every radius; the pressure-tail Uz_Z atom survives. This replaces the legacy frozen prescribed-value policy on this local route. Independent resolved scalar RK4 agrees with the continuation fields and five moments to 2.93082671725e-13 maximum scaled difference. In that resolved fixture, directly differentiated analytic primitives give maximum atom defects about 3.57e-102 for the similarity divergence numerator, 1.36e-106 for dMz/dR-Uz, and 2.73e-107 for the axial slope/driver identity. These are local finite-ring checks, not independent Cartesian or global certificates.

## Reproduction

- `python experiments/root_st073/lei_ren_part1_paper_pressure_width_continuation_fixture.py`: independent resolved scalar RK4 replay of the prescribed continuation, checking F, Uz, P and all five moments. The auxiliary frozen adapter separately reproduces the original comparison at s=3.
- `python experiments/root_st073/lei_ren_part1_paper_pressure_width_continuation_check.py`: actual degree18, pressure9/width2 source, delta=1e-200 and h_b=exp(-100-1e154), with R=1,10,100 probes.
- `python experiments/root_st073/lei_ren_part1_paper_exponential_polynomial.py`: independent quadrature and derivative checks of the integration algebra, including tiny mixed atoms.

These are local finite-ring calculations. Finite pressure/width/core truncation, source integration, collar RK and complete axial bounds are not certified. They are not temporal coefficient recursion. Switching at R=100..110, long reshape, functional five-moment closure, finite global energy, exact/controlled heat and stress-cone closure still remain.
