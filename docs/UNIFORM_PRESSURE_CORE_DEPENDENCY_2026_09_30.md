# Uniform pressure-to-core dependency — 2026-09-30

## Functional dependency completed for the finite core

The pressure-driven finite-core perturbation is now enclosed at every Z in [-0.8,0.8], at the core boundary R=4/Lambda. This extends the prior single-center error receipt and uses the same stored pressure datum, not a reset pressure anchor.

The actual stored-parameter degree18 run gives:

| Quantity | Uniform absolute perturbation upper |
| --- | ---: |
| Uz | 4.394371482111663e-36 |
| Uz_Z | 2.822870316102291e-35 |
| Uz_ZZ | 2.360245671266824e-34 |
| Ur | 4.516471968020449e-53 |
| Ur_Z | 4.010235089588441e-52 |

All five finite-core radial moments have value, first-Z, and second-Z pressure perturbation intervals in the JSON receipt. These are only the core contribution to the full defect vector. They are not terminal five-moment closure or complete core approximation error.

## How the interval conclusion is obtained

`lei_ren_part1_paper_uniform_axis_jets.py` encloses the real axis construction analytically. The anchor root lies in [-j,0]. On the real path, |L H/(H squared+sigma squared)|<=sup|L|/(2 sigma), so G is bounded by path length a+j times this envelope. F0=exp(-logC-Lambda G) is enclosed without evaluating G by quadrature. Rational interval jets recover its higher derivatives. The constant denominator coefficient retains the exact nonnegative square H squared rather than losing its dependency in a generic interval product.

`lei_ren_part1_paper_interval_difference.py` propagates a nominal interval and its perturbation separately. A product difference is x0 dy+y0 dx+dx dy; a quotient difference is (dx-(x0/y0)dy)/(y0+dy). Common uncertain axis factors have exactly zero pressure perturbation, avoiding an artificial width from subtracting two independently evaluated wide interval ranges.

`lei_ren_part1_paper_uniform_core_pressure_propagation.py/json` uses this paired arithmetic in the existing nonlinear radial recurrence and exact polynomial moment integrals. It consumes the 24-order uniform pressure error receipt. Axis constants are the stored MP values at precision473, and the pressure schedule is the one returned by the current `source_profile()` pressure-receipt helper. This calculation does not install a new accepted global field or reconcile a different coherent-waiting schedule automatically.

## Flatten coefficients

`lei_ren_part1_paper_uniform_flatten_error.py` now provides a separate directed coefficient-and-tail enclosure. The actual wrapper integrates normalized u in [0,1] with the exact stored length100, retaining the slope -1/2-mu and interval sigma values. It avoids converting huge radial origins or exact Decimal mu into rounded nominal inputs before enclosure.

At coefficient order96 and 512 endpoint panels, the actual relative uniform value error upper is approximately 0.8592 against the true Z0 mass lower bound. This is valid but coarse and does not establish relative-flat closure. The whole pressure absolute budget remains much tighter because the flatten mass is extremely small. A sharper enclosure integrating the affine exponential on each panel is the next refinement. The nominal MP-only coefficient API remains explicitly conditional on unbounded roundoff and is not a strict certificate.

## Checks and open dependencies

The paired arithmetic fixture verifies 36 independent perturbation samples, exact-zero propagation, denominator rejection, and recurrence compatibility. The axis fixture checks 80 independently constructed point jets inside the analytic uniform bounds. The actual core fixture checks 76 local field/moment perturbations at four axial positions inside the uniform result. Point checks are regressions; the uniform conclusion comes from analytic bounds and directed operations. Existing component-core and moment-derivative regressions pass.

Still open: errors in original parameter derivations; agreement with any separately constructed coherent-waiting global source; infinite radial-series remainder; complex analytic core norms; RK continuation; annular and source-integral approximation errors; full five-defect function bounds and inverse closure; exact heat velocity exterior; admissible stress cone; time-dependent higher-order recursion and oscillatory correction. These remain necessary for the full goal.
