# Axially shifted nested swirl experiment

Retain all six successful radial rings and append six rings centered at z²=.8. Additional (radial center,width) pairs are (.16,.12),(.35,.2),(.7,.35); axial widths are .25 and .5. There are 36 bounded coefficients in [-1,1], with the original 18 coefficients retained as a warm start and the additional 18 initialized to zero. The parent field, pressure and force remain frozen. The shared outer cutoff, moment quadrature, serialization, energy integration and quadratic residual cache are reused.

## Development results

| Run | Calls | Independent sampled maximum | Termination |
| --- | ---: | ---: | --- |
| axial_swirl | 222 | 1.34996179 | `ftol` termination condition is satisfied. |
| axial_swirl_fourth | 500 | 3.63158700 | actual 500-call budget reached |
| axial_swirl_dense | 500 | 1.28545098 | actual 500-call budget reached |

All three pass structural probes and fail the .001 PDE threshold. The sparse fourth-power trial is a preserved failure: stronger pointwise weighting did not fix inadequate spatial coverage. The dense run appends 24 radial x 24 axial x 9 temporal collocation samples to 2048 uniform training samples, without using validation samples. Its objective is not a physical volume-weighted L2 norm; independent validation retains the original physical metrics.

Select axial_swirl_dense as the new development reference: maximum 1.28545098, energy [.69432588,1.0], core drift .04949132. This is about 5.2% lower than localized_swirl at the same derivative step. The optimizer hit its 500-call cap; convergence is not claimed. Final blind acceptance remains pending, and sampled maxima do not bound continuum maxima.

One focused regression verifies the nested basis reproduces the previous candidate, initial-field preservation, numerical zero added moment, serialization, and coefficient bounds. It passed with warnings treated as errors. The training cache is compared to direct residual evaluation on a nonzero coefficient vector before every fit.

Reproduce with constrained_temporal_optimize.run(output=..., axial=True, fourth_power=True, dense_cylindrical=True), then the standard independent validation CLI using configs/constraints_axial_swirl.json. Set the two comparison flags false as appropriate for the other saved trials.

Next: inspect the residual distribution on the denser collocation grid and determine whether the remaining bottleneck is azimuthal evolution or radial/axial pressure balance. Preserve the sample coverage improvement in subsequent fitting; do not relax bounds or acceptance thresholds.
