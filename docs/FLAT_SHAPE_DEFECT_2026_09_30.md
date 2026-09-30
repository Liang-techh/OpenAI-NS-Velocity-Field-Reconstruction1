# Flat reshape defects in signed logarithmic coordinates

Terminal cumulative-moment subtraction is not a reliable defect calculation at the actual radii. For the reshape define ell=T-y, q=sigma(ell/T), and B=log(Cstar Utheta110(1+Z^2)). The exact ratio to the reference swirl is exp(q B). The three angular defect sources reduce to

    I(k,m) = integral_0^T exp(-k ell) expm1(m B q) d ell

with (k,m)=(8/5,1),(1/5,2),(6/5,2) for angular, pressure and angular energy terms. Their normalized prefactors are alpha^(8/5), alpha^(1/5)/2 and -alpha^(6/5)/2, where alpha=Rsh/Rm. The actual profile B is inherited, not fitted.

When q is tiny, log(abs(expm1(m B q))) is approximately log(abs(m B))+log(q). The log source is therefore approximately

    -k ell - T^2/ell^2 + 1 + log(abs(m B)).

It has a mode near (2 T^2/k)^(1/3), and width near sqrt(ell_mode/(3k)). This mode lies far beyond the previous ordinary moment quadrature cutoff when T=4e152. Using the derivative mB q_prime alone misses this mode: the derivative of the logarithm is asymptotically q_prime/q.

The new kernel should retain sign and log magnitude, use stable logistic log(q), locate the source mode, and integrate a normalized local coordinate. Its B sensitivity provides the first axial derivative when multiplied by inherited B_Z. Resolved full-interval scalar quadrature is required for comparison.

## Actual leading-atom receipt

The check reads the committed actual R110 reshape_B and reshape_B_Z leading atom from pressure_width_switches_check.json, not a fitted B. At T=4e152 the angular, pressure and angular-energy sources are all nonzero in signed-log representation. Source modes are approximately5.85e101,1.17e102 and6.44e101, far beyond the old cutoff1842. The independently derived leading saddle log formula differs from the numerical kernel by about2.92e-51,5.85e-51 and3.22e-51. These comparisons do not enclose the finite-window tail or full source error.

Reproduce:

    python experiments/root_st073/lei_ren_part1_paper_flat_shape_defect_check.py

The receipt covers only the serialized leading B/B_Z atom and three reshape source terms. It does not represent the completed five functional defects.

## Independent resolved evidence

Full-interval scalar mp.quad at T=400, both B=-4 and B=4, checks all three kernels and their direct derivative integrals. Maximum log-integral discrepancy is1.99e-17 and relative first-B sensitivity discrepancy1.24e-71. Signs are retained and the defects/sensitivities remain nonzero. Finite-T saddle location/curvature are additionally compared with the asymptotic model; this is not an exact finite-T identity.

    python experiments/root_st073/lei_ren_part1_paper_flat_shape_defect_fixture.py

For pressure/width component lifting, nth B sensitivities require separate kernels integral exp(-k ell)(m q)^n exp(m B q) d ell. Their small-q modes scale as (2 n T^2/k)^(1/3); reusing the n=1 center would miss higher-order contributions.

## Scope

This kernel is a building block for a functional defect provider. Scalar leading-atom evaluation alone does not recover all pressure/width components or all five defects. Finite-window and quadrature errors are not enclosed; source and core errors remain inherited. It does not certify terminal moments, stress admissibility, global energy, heat or temporal recursion.
