# Analytic parameter Jacobian and training-grid continuation

The cached residual is quadratic in the added velocity coefficients. Differentiate that polynomial directly, including both convective cross terms. The cached energy quadrature is also quadratic; its derivative supplies the energy-penalty rows. The fourth-power residual-vector transform uses the chain rule. These are parameter derivatives of the existing finite-difference training operator, not analytic spatial derivatives and not independent validation.

The new optimizer options are analytic_jacobian=True, coefficient_initial=<saved candidate>, and grid_resolution=<radial/axial count>. Warm starts verify parent identity. Two focused cache tests passed, comparing residual/energy parameter directional derivatives with central differences and cached energy with the original quadrature.

Starting at axial_swirl_dense on the same 24x24x9 collocation grid took7 function calls and stopped on ftol. Objective changed from .020337419715 to .020337401490; random development maximum changed from1.28545098 to1.28543627. This establishes local optimizer stagnation for the fixed basis/objective, not global optimality.

A separate64x64 cylindrical diagnostic at t=.748 exposed a peak1.62313085 near(r,z)=(.39952,.90683), dominated by the swirl component. The random maximum was not a continuum bound. Increasing training resolution to48x48x9, retaining the2048 uniform samples, reduced that grid peak to1.29826621. The random validator maximum slightly worsened to1.29013391; both facts must be reported. Structure samples pass, energy[.69829101,1.0], core drift .04949132. Neither candidate approaches the.001 acceptance threshold.

Use axial_swirl_grid48 as the next working reference because it controls the known missed peak better, while retaining axial_swirl_analytic as the lower random-sampled-maximum comparison. This choice changes no acceptance threshold. The64-grid informed the training-grid choice and is development evidence, not a blind audit.

Reproduce continuation with constrained_temporal_optimize.run(axial=True, fourth_power=True, dense_cylindrical=True, analytic_jacobian=True, coefficient_initial='artifacts/constrained/axial_swirl_dense/candidate.json', output='artifacts/constrained/axial_swirl_analytic'). For grid48 use grid_resolution=48 and coefficient_initial pointing to axial_swirl_analytic, with output axial_swirl_grid48. Validate using the original constrained_validation CLI and configs/constraints_axial_swirl.json.

Next: change the spatial/temporal representation around the core-exclusion collar or solve azimuthal evolution directly. Continuing the same36 coefficients at the same grid has demonstrated diminishing returns. Keep dense coverage and analytic parameter derivatives for subsequent candidates.
