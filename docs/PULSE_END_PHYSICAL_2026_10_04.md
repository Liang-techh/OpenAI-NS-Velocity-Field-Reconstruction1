# Original pulse-end physical tensor and three-component remainder — 2026-10-04

The original end chart s in [-4,0], Z in [-1,1], including both beta supports centered at -3 and -1 with width .15, now has a physical stress companion for the SAME full meridional source. The current five-moment, absolute-pressure and similarity-stress interface at s0 / flatten t0 is composed into the physical tensor, divergence and remainder interface.

Run the bounded layer:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulseendphysical
```

The producer/checker pair is `lei_ren_part1_paper_compliant_pulse_end_physical_C2.py` and `lei_ren_part1_paper_compliant_pulse_end_physical_C2_check.py`. Original source data, selected coefficients, signed incoming memory, complete future energy half, analytic absolute pressure, radii and positive log scales are preserved.

## Physical construction

For fixed nu > 0, use x_source = x / sqrt(nu), u_physical = sqrt(nu) u_source, p_physical = nu p_source, and unchanged time. The actual coordinate relation is

- lambda^2 - lambda^(2 delta) z^2 / nu = tau = 1-t;
- R = r^2 / (2 nu lambda^2);
- Z = z / (sqrt(nu) lambda^(1-delta)).

The symmetric completed tensor has Trtheta = Ttheta, Trz = Tz, and Ttheta_theta = r partial_z(Tz); Trr and Tzz are zero. Angular/axial stress use nu lambda^(-2-delta). The diagonal completion cancels radial stress divergence exactly.

Full background momentum is -div(T) + E. The three cylindrical error components are

- Er = partial_t(ur) + ur partial_r(ur) + uz partial_z(ur) - nu [partial_rr(ur) + partial_r(ur)/r - ur/r^2 + partial_zz(ur)];
- Etheta = -nu partial_zz(utheta);
- Ez = -nu partial_zz(uz).

The pressure/centrifugal cancellation in Er uses the actual same-source pressure FTC, not a fitted pressure. Er is generally nonzero inside the pulse. Reusing downstream pure-swirl zero Er/Ez would lose the original meridional dynamics.

The implementation retains four radial error sectors: fixed-position time, radial viscosity, nonlinear meridional transport and axial viscosity. Their lambda exponents are -3, -3, -3 and -3+2delta; angular/axial errors use -3+delta. All errors and stress divergence scale as sqrt(nu); tensor entries scale as nu. Spatial derivative viscosity powers are applied afterward.

## Source and derivative evidence

Focused checker PASS: 411 current source hashes, 900 finite signed physical rows and 260 structural-zero rows. The independent comparison tolerance is 1e-55; maximum positive enclosure miss is 1.18849154463436673139964949068e-88.

There are 28 symbolic identities replaying the original full equations (3.16)-(3.18): radial recovery, full angular/axial inertia, radial viscosity, centrifugal cancellation, incompressibility and tensor product derivatives. The source/physical interface bridge supplies 18 bindings and consumed functional identities.

The output provides stress mixed3, completed diagonal/divergence mixed2 and all three error mixed2 as signed sectors with exact source logs. Unmaterializable factors B, D, H and R are never chosen from cap intervals as field values. Whole-Z evaluation preserves 1+Z^2 >= 1 before reciprocal Taylor algebra, avoiding a spurious zero denominator from independent interval multiplication.

At s0 the actual original beta supports and backward meridional moments are empty: Er and Ez mixed2 vanish, and the generally nonzero Etheta equals the flatten angular axial-viscosity error. This is a source-functional physical join; numerical interval overlap is not its proof.

An independent smooth full-meridional fixture differentiates the original physical velocity, absolute pressure and completed tensor in Cartesian coordinates. It includes fixed-position time, convection, pressure gradient, full vector Laplacian, moving cylindrical basis and tensor divergence. Both nu=.01 and nu=.7 retain nonzero Ur/Uz and all three errors. The 114 comparisons cover stress mixed3, diagonal/divergence/error mixed2 and incompressibility. This fixture is local formula evidence, not the actual source cone, global error or final corrected NS accuracy.

The read-only GPT-5.6 Luna/max review confirmed original remainder equations and found no material sector-unit/source gap.

## Scope and next work

This completes the original end-chart physical companion and its flatten endpoint. It does not complete all pulse regions, prove a cone on either beta support, or establish flatness/energy/recursion.

1. Close the four original beta support boundaries -3 ± .15 and -1 ± .15 with source flat-tail derivative estimates and partial-primitive continuity. Bound all required ordinary-logR and axial derivatives uniformly.
2. Prove continuous whole-end admissibility through both supports and the intervening region, retaining signed moment/energy histories and the current selected coefficient correlations. Do not use a sample grid or raw interval subtraction as a cone proof.
3. Extend full moments, stress, physical errors and cone to gap, main and entrance charts with their original radii and functional joins.
4. Resolve upstream finite-width bridge feedback, then completed full-tensor/global admissibility.
5. Independently establish global flat error, physical-volume norms and required-domain energy; implement actual n-dependent recursive recovery, oscillatory cancellation, corrected Cartesian residual and measured vortex/particle dynamics.

All corresponding broader gates remain false.
