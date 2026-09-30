# Five centered component defect functions

This module assembles the Section10 centered defect inputs from the common R110 field, its five moments, analytic axis pressure, and first Z tangents. It computes interval sources rather than terminal total-minus-reference differences.

Let g=Uz110-4Z, Am=exp(logPstar-.6)/(1+Z^2), rho0=110/Rm, alpha=Rsh/Rm. Rm=exp(-6)Rref, Rz=exp(-8)Rref. The flat angular kernels from the preceding component adapter contribute to d3, d5, the angular part of d4, and the mixed part of d2.

Constant axial-source contributions through Rz are

    d1 += g (exp(-2)-rho0)
    d2 += g [(exp(-3.2)-rho0^1.6)/1.6 + alpha^1.6 I_theta]
    d4 += (g/Am)^2 (exp(-2)-rho0)

Restoration adds

    d1 += g exp(-2) integral_0^1 exp(t)(1-sigma(t)) dt
    d2 += g exp(-3.2) integral_0^1 exp(1.6t)(1-sigma(t)) dt
    d4 += (g/Am)^2 exp(-2) integral_0^1 exp(t)(1-sigma(t))^2 dt

The reference mixed weight is rho^(8/5): rho^(3/2) must also include u_ref/Am=rho^(1/10). The axial sources vanish after restoration, and angular sources vanish once the reference swirl is reached. Inner seed/reference contributions at R110 remain separate, as do each flat derivative-order term and each axial source. First Z derivatives differentiate these same formulas, including the Am dependence.

The axis pressure is recorded and unchanged. Constructing these input defects does not itself change the velocity or set any terminal moment to its target.

## Actual common-source receipt

`lei_ren_part1_paper_centered_component_defects_check.json` evaluates the actual common R110 source at Z=.3, pressure order 9 and width order 2. All five baseline defects are nonzero, with 69 separate source labels. Approximate baseline values are d1=2.23437e-15, d2=5.68986e-16, and d4=5.91501e-41. The negative d3 and d5 have log absolute magnitudes approximately -1.6e152 and -2e151 respectively. The receipt retains arbitrary-exponent values and first-Z tangents rather than rounding these last two rows to zero.

This is a fresh source evaluation, not a target overwrite. These numbers are inputs for a correction, not evidence of terminal closure or of uniform analytic smallness. The next repair must retain the separate parts when multiplying defects and building nonlinear interaction terms.

## Scope

A complete finite-ring defect formula is different from certified functional closure. Inherited source/core errors, finite pressure/width truncation, flat-window quadrature and restoration quadrature are unenclosed. The paper's uniform analytic norms and quantitative smallness input still require verification before the five-bump repair. Stress cone, finite energy, exact heat, temporal recursion and full Cartesian NS validation remain open.

## Independent resolved density replay

The fixture independently integrates the centered densities through the reshape, reference and axial restoration intervals, using stable expm1 angular differences. Its maximum scaled row and first-Z errors are 2.96e-15 and 1.50e-14; relative errors for the tiny d3 and d5 rows are 6.17e-25 and 6.11e-25. The maximum relative row/tangent diagnostic is 2.04e-14. P0 and P0_Z are preserved exactly. This is a resolved surrogate with finite pressure/width data, not an enclosure of the actual extreme-parameter source.

A downstream solve must not discard separate defect responses when larger coefficients cancel in a tiny row. Ordinary fixed-precision total subtraction cannot validate closure relative to defects whose log magnitudes are of order -1e152. Preserve the source hierarchy or a formal defect expansion and independently replay moment changes.
