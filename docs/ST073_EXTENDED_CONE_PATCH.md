# Relaxed-cone location on the extended ST073 mean

The [OpenAI Navier–Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
requires a strict relaxed stress cone on a mean transition profile before
its periodic shear loop and oscillatory corrections can realize the
needed flux (Section 4.2, Appendix C, and Sections 7–9). This experiment
applies the paper's normalized snapshot formulas to our current extended,
time-dependent physical candidate. It is a **diagnostic analog**, not a
proof that our field is a paper-admissible leading profile.

`extended_relaxed_cone_screen.py` integrates all five radial moments using
the candidate's actual radial interfaces: the inner edge, the old radial
bridge edge, and the moving radius beyond which the field is pure axial
heat swirl. The old cone script used the previous fixed interfaces and
could not be applied directly here.

At `k=5.5`, `eta=.2,.3` and each sampled
`X=1.5,1.75,2,2.25,2.5`, the relaxed inequality has positive sampled
margin. The smallest of these ten margins is about `0.00584`. At
`X=2, eta=.2,.3`, it remains positive at `k=3,5.5,6`; the smallest of
those six margins is about `0.0483`. Increasing radial quadrature order
from 8 to 20 at `X=2,eta=.3,k=5.5` changed the margin only from
`0.0510` to `0.0514`. These facts suggest a local open cone patch,
but no continuous-space/time certificate exists.

The cone patch starts at the old radial bridge's outer edge, not inside
its main momentum defect. At `eta=.2,k=5.5`, fourth-order physical
finite differences give complete residual norms about `2379` at
`X=.5859375`, versus `8.71e-4`, `9.39e-4`, and `9.95e-4` at
`X=1.75,2,2.5`. Adding a wave only to the currently cone-admissible
outer patch cannot directly cancel the bridge-local residual and could
push the nearly passing outer patch above the `1e-3` target.

Inside the bridge at `X=1,eta=.2`, the baseline normalized shear value
`v_s=13.62` exceeds the relaxed upper bound `3.37`. A symmetric
`64 y^3(1-y)^3` azimuthal bubble does not produce a pass in the tested
amplitude range. An outer-biased azimuthal bubble with amplitude `.76`
reduces the radial swirl-slope measure `a` from `7.16` to `2.03`, but
`v_s=5.83` still exceeds its new upper bound `2.88`; at `eta=.3` the
gap is larger. On four physical momentum points, the maximum residual
increases from about `3947` to `4419`. Both shapes preserve two radial
derivatives at the bridge endpoints and analytic incompressibility, but
neither is an accepted correction.

Holding the outer-biased trial's other cone quantities fixed as a local
diagnostic, its `v_s=a+b_s^2/a` would require the magnitude of `b_s`
to fall by at least about **53%** at `eta=.2` and **82%** at `eta=.3`
to reach the current upper bounds at `X=1`. Because `b_s=2X U_X/E`, this
points to a substantial meridional-shear change in addition to the swirl
slope repair. A real poloidal change will also alter the stress vector
and upper bound, so these percentages are not a standalone solution.

The next construction must jointly alter the bridge's swirl slope,
meridional shear, and pressure/moment data while preserving the five
outgoing radial moments. Only after a strict cone exists where the
residual actually needs cancellation can the paper's shear loop and
nonaxisymmetric amplitude dynamics be tested. Reproducible reports are
`extended_relaxed_cone_screen.json`, `swirl_bridge_cone_screen.json`,
and `outer_swirl_bubble_screen.json` under `experiments/root_st073`.

The subsequent solenoidal poloidal trial opens the cone at two isolated
bridge points but makes complete momentum much worse; see
[`ST073_COUPLED_BRIDGE_FAILURE.md`](ST073_COUPLED_BRIDGE_FAILURE.md).
