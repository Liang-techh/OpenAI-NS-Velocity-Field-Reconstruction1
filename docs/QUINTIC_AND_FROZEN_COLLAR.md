# Quintic time extension and frozen-collar diagnosis

Cubic Bernstein coefficients are elevated exactly to degree five using positive binomial weights. All old candidates are contained in the new 60-coefficient family with the same [-1,1] bounds. The spatial modes, parent velocity parameters, pressure and restricted force remain fixed. Tests verify velocity and time-derivative equivalence of the elevated start and artifact round-trip. Two focused axial/quintic tests passed.

## Development fits

| Fit | Calls | Full independent sampled maximum | Structure probes |
| --- | ---: | ---: | --- |
| quintic_swirl | 19 | 1.26452244 | pass |
| quintic_swirl_scalar | 6 | 1.26947752 | pass |
| quintic_swirl_scalar_fourth | 12 | 1.27287148 | pass |

The full-momentum fourth-power fit is selected: maximum1.26452244 versus cubic1.29013391, with energy[.70555030,1.0] and unchanged core drift .04949132. Scalar comparisons optimize only azimuthal residual, using respectively mean-square and fourth-power objectives; neither improves the full validation maximum. All still fail .001. Training uses the same48x48x9 plus2048 uniform samples and analytic parameter Jacobians. These are development comparisons, not a blind audit.

## Decisive representation restriction

Every current correction mode contains the parent's outer basis. That basis is exactly zero for r²<=.125, so no coefficient in the temporal, localized, axial or quintic correction families changes velocity or its derivatives in that open region. Pressure and force are also frozen.

At (r,z,t)=(.32,.8857142857,.75), independent residual norms at h=.005,.0025,.00125 are1.15299872,1.15158477,1.15115944. Candidate and uncorrected parent residuals are identical at these stencils, all contained inside the zero-correction region. See quintic_swirl/frozen_collar.json. The support argument is exact; the residual values are numerical with refinement, not certified interval bounds. They demonstrate a large untouched defect that further fitting of these same corrections cannot repair.

Next CR003/CR005 task: introduce additional divergence-free pure-swirl modes that reach the annulus r>.1 while still vanishing on an open neighborhood of all core probes (maximum core radius .1*sqrt(.75)=.086603). For example use a smooth guard vanishing for r<=.1 and reaching full value by about.142, with compact outer support. Preserve the old modes and warm start. Orthogonalize each new mode's angular moment against the existing outer torque basis, retain coefficient bounds, force, energy/core thresholds, and dense collocation. Do not alter the parent field or merely increase the old polynomial degree again. If radial/axial residual then dominates, jointly optimize pressure/poloidal structure using full validation.

Reproduce the full fit with constrained_temporal_optimize.run(quintic=True, fourth_power=True, dense_cylindrical=True, analytic_jacobian=True, coefficient_initial='artifacts/constrained/axial_swirl_grid48/candidate.json', grid_resolution=48, output='artifacts/constrained/quintic_swirl'). Scalar comparisons add azimuthal_only=True and use the quintic candidate as their warm start. Use configs/constraints_quintic_swirl.json for independent validation.
