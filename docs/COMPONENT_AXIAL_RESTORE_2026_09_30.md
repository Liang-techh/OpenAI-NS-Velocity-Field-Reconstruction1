# Component axial restoration — 2026-09-30

Section9.38 restores the common Rz endpoint axial field v1(Z) to 4Z over t=log(R/Rz) in [0,1]. The angular reference remains Utheta=Utheta_z exp(t/10). After t=1 the axial field is 4Z, continued to Rh at t=3. Prior inner, reshape and reference moment parts are retained separately.

The three bounded kernels are

    j1(t) = integral_0^t exp(-(t-s)) sigma(s) ds
    j16(t) = integral_0^t exp(-1.6(t-s)) sigma(s) ds
    j2(t) = integral_0^t exp(-(t-s)) sigma(s)^2 ds

These determine the prescribed axial mass, mixed and quadratic moment increments. Angular, pressure and swirl increments use exact power primitives. Automatic first Z tangents differentiate the same fields and moments. Total pressure changes only by its angular pressure-moment increment; the inherited analytic axis P0 stays fixed.

## Actual-source evidence

Degree18 pressure9/width2 Z=.3 replay completes phases .5,1,3 (through Rh). All recorded entry field discrepancies are zero. At phases1 and3, Uz=1.2 and Uz_Z=4 within the actual check threshold1e-200. Original inner, reshape and reference moment parts and P0 are unchanged. The field now includes separate axial_restore and restored_reference moment contributions.

Reproduce the actual source with:

    python experiments/root_st073/lei_ren_part1_paper_pressure_width_axial_restore_check.py

The optional --resume uses the locally generated ignored R100 cache. Fresh mode reconstructs it from the common source.

## Independent resolved evidence

Independent scalar adaptive quadrature at phases .5,1,3 gives maximum physical integral discrepancy5.99e-11, first-Z stencil discrepancy4.16e-9 and separate-part receipt discrepancy2.54e-9. Join/entry discrepancy is zero at recorded precision. The radius is a fixed shared parameter; its first-Z tangent is checked structurally as zero rather than differencing enormous rounded scalar radii. These are resolved fixture discrepancies, not actual-scale enclosures.

    python experiments/root_st073/lei_ren_part1_paper_pressure_width_axial_restore_fixture.py

## Limits and next repair

Bounded Gauss quadrature and inherited core/reshape/finite-ring errors require enclosure. First Z tangents do not establish uniform C2. No complete global velocity field, finite-energy tail, stress cone, exact heat or temporal recursion follows from this interval alone.

For Section10 functional moment repair, compute defects from separate interval contributions rather than subtracting enormous full cumulative totals. In particular, the reshape receipt's normalized tail envelope is an error term, not a resolved terminal defect. Flat endpoint differences can lie below the working precision at the final quadrature nodes. A recorded zero from total-minus-reference subtraction cannot establish an exact moment identity. The defect provider must retain source terms, reference differences and unresolved quadrature/tail error contributions separately, then check the paper's quantitative smallness input before solving the five corrections on [Rm,2Rm]. Axis P0 must not be reset after that solve.
