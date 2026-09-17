# Localized swirl and loss comparison

The previous global polynomial correction reduced mean residual but did not reduce the maximum. A separate 80-by-80 cylindrical training grid at t=.748 (r in [.03,1.95], z in [0,1.95]) located the largest residual near r=.4675, z=0; azimuthal residual was about1.855 for temporal_swirl. This is a new development diagnostic, not final blind validation.

Use six Gaussian ring factors times the existing outer swirl basis, minus each mode's quadrature angular-moment projection. The outer cutoff preserves compact support and core exclusion; each nonzero cubic Bernstein temporal factor preserves the initial field. All18 coefficients stay in[-1,1]. Pressure and restricted force stay fixed.

The candidate has smooth axisymmetric pure-swirl corrections and hence exact symbolic divergence cancellation. Zero angular moment is approximate to the quadrature used in the stored construction. Full PDE residual and energy are checked separately.

Compare three fits on identical uniform training points and times: previous polynomial basis with fourth-power loss, localized rings with mean-square loss, and localized rings with fourth-power loss. The extra polynomial comparison distinguishes effects of changing the objective from effects of spatial localization. Fourth-power means minimize mean(||R||^4), implemented as a least-squares residual vector R*||R||; validation thresholds and metrics remain unchanged.

Reproduction uses constrained_temporal_optimize.run(output=..., localized=True/False, fourth_power=True/False), then the standard constrained_validation CLI. Each run has a500 actual-call cap and uses the independent development validator. No run is final blind acceptance.

## Results

| Fit | Actual calls | Independent sampled maximum | Structure probes |
| --- | ---: | ---: | --- |
| Polynomial, fourth power | 248 | 2.03411666 | pass |
| Localized, mean square | 114 | 1.35622595 | pass |
| Localized, fourth power | 350 | 1.39127608 | pass |

Localized mean-square is selected for subsequent development: maximum improves about30.2% from continued_pressure (1.94296483). Its training mean-square is0.06195655; energy range[.72793324,1.0]; initial velocity unchanged and core drift .04949132. All candidates still fail .001. Losses with different powers are not directly comparable.

Derivative refinement on the existing development samples gives maximum1.35407291 at h=.0025 and1.35345289 at h=.00125, supporting a plateau near1.353 rather than a residual explained by finite differences. This is not final blind acceptance.

A separate80x80 cylindrical training grid at t=.748 now peaks near r=.39456,z=.88861, with residual(-.43747,-1.38432,.34311) and norm1.49179. Random sampled maxima do not bound the continuum maximum. No fitted coefficient is at its bound. This supports adding an axially shifted ring or direct azimuthal evolution, rather than widening coefficient bounds. Preserve the current trial as the reference and do not close CR006/CR012.

Artifacts: temporal_swirl_fourth/, localized_swirl/, localized_swirl_fourth/. Localized result includes diagnosis.json and refinement.json. Four focused localized candidate tests passed with warnings treated as errors.
