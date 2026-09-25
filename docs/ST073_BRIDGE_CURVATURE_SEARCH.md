# ST073 bridge curvature and three-mode search

The extended core still leaves a physical momentum defect in its radial
bridge. The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
requires the mean transition to have suitable cumulative moments and a
strict relaxed stress cone before oscillatory flux can correct this defect.
This search is a finite-dimensional diagnostic toward that construction,
not a substitute for the paper's dynamically evolved mean or waves.

`bridge_poloidal_mode.py` now includes a `minimum_curvature` streamfunction
shape. Within the chosen polynomial family it minimizes
`integral_0^1 |phi_yyy|^2 dy`, holding `phi_yy` at the reference `X=1`
equal to the prior curvature mode and preserving the radial endpoint
jets. Its third-derivative L2 falls from about `300.45` to `35.37`.
On the same eight physical bridge points, replacing the old poloidal
shape cuts the sampled maximum complete momentum residual from about
`43927` to `5972`. The minimum-curvature shape alone loses the relaxed
cone because its cumulative stress vector has `P_c<2` at the two
reference slices. These results are in `minimum_curvature_bridge_screen.json`.

A second streamfunction mode has a triple zero at `X=1`: it leaves the
local `U`, `U_X`, and cumulative `M` there unchanged while altering other
upstream radial moments. Combining it with the minimum-curvature mode
can restore a positive two-point cone, at a momentum cost. The complete
residual is quadratic in the three mode amplitudes because the field and
all linear PDE terms depend affinely on them. The sampled quadratic
response in `bridge_quadratic_collocation.py` screened 405 grid candidates
and cone-checked the 30 with lowest sampled residual.

The best cone-positive grid candidate in that limited search has
`(outer_swirl, minimum_curvature, moment)=(1,-2.75,4.5)`. Direct physical
finite differences reproduce its eight-point maximum `6061.5`; the
unmodified extended candidate is `3946.7` on those same points. At
`X=1, eta=.2,.3`, the sampled relaxed cone remains positive at
`k=3,5.5,6`, but at `k=5.5` it fails by `X=.75` and `X=1.1` on both
slices. The full record is in `bridge_collocation_holdout.json`. This
is **not** a continuous cone region or a momentum improvement over the
baseline, and the five outgoing moments are not restored.

The concrete remaining design problem is a constrained mean-flow solve:
reduce the bridge's radial-viscous curvature and complete momentum at
every relevant radius and time while maintaining a strict cone over an
interval, matching five radial moments at the outgoing interface, and
preserving the core and pure heat exterior. A pressure/mean correction
must be part of that solve. Only then does the paper's nonaxisymmetric
oscillatory-stress stage have a viable background. No maximum or volume
L2 gate below `1e-3` has been demonstrated.

## Fifth-moment obstruction and internal balancing experiment

`bridge_swirl_moment_balance.py` compares the three-mode candidate with
the unmodified extended field at the old radial exit `X=1.5`, at
`tau=.5*2^-5.5`. The fifth moment is the nonnegative integral
`Cp(X)=integral_0^X E^2/(2x) dx`. The outer-swirl-1 candidate adds
`0.5294`, `0.5242`, and `0.4902` to `Cp` at `eta=.2,.3,.65`.
The **entire** unmodified downstream `Cp` tail, from `X=1.5` through
the attached heat exterior to infinity, is only about `0.02246`,
`0.02247`, and `0.02249`. Thus a correction restricted to `X>1.5`
cannot restore the unmodified outgoing fifth-moment target while keeping
the same exterior: even deleting all downstream swirl leaves a positive
defect greater than `0.46`. This is a necessary obstruction for that
specific downstream-only repair, not an obstruction to a redesigned
bridge or a different mean/exterior target.

The experiment next adds the existing symmetric bridge-swirl bubble with
amplitude `-0.2731162731` alongside outer-swirl amplitude `1`. This
single coefficient, chosen between the small `Cp`-neutral roots at
`eta=.2` and `.65`, reduces the three sampled fifth-moment changes to
`-0.01693`, `-0.01267`, and `+0.01454`. The positive change is now below
the downstream capacity on these slices. The sampled relaxed cone at
`X=1`, `eta=.2,.3` remains positive. The eight-point complete momentum
maximum is about `6111` versus `6062` for the prior three-mode candidate
and `3947` for the unmodified field. The other outgoing moments still
move substantially (for example `I` by `0.5231` and `S` by `3.026`
at `eta=.2`), and no continuous cone or PDE gate is met.

Next solve the bridge modes as a *coupled five-moment and momentum*
problem. Internal sign-changing swirl has removed one simple fifth-moment
impossibility but does not cure the radial-viscous defect. Any downstream
restoration must also account for the remaining four moments and its own
physical momentum cost.

## Pressure projection and angular velocity screen

`balanced_bridge_pressure_screen.py` fits six additional smooth pressure
modes to the balanced four-velocity-mode candidate, using complete physical
momentum on eight bridge points. The direct maximum rises from `6111` to
`6310` on training points, from `8425` to `8454` on disjoint spatial
points, and from `6495` to `6518` at a different time. The analytic
pressure-gradient prediction agrees with finite differences to below
`8e-7`, so the failure is not a gradient implementation error. The
pressure modes cannot affect the large axisymmetric angular momentum
component: its spatial holdout maximum is `4745`.

`angular_bridge_mode_screen.py` then treats the two swirl and two poloidal
amplitudes as a quadratic physical residual response. Its deterministic
search enforces the necessary downstream fifth-moment capacity on three
axial slices. It reduces the eight-point angular maximum from `3016` to
`2067`, with a disjoint spatial reduction from `4745` to `3721`.
However, its complete momentum maximum **increases** from `6111` to
`7634` on training points and from `8425` to `9962` on holdout points.
Both sampled relaxed-cone margins at `X=1` become negative. Thus the
existing four-mode family trades angular transport against meridional
viscosity/cone geometry; angular-only optimization is not an acceptable
mean-field correction.

The next mean-flow construction needs additional radial/axial degrees of
freedom that directly control the pressure-independent angular transport
while keeping poloidal curvature low. Optimize complete momentum,
five outgoing moments, and a *continuous* relaxed cone together. These
screens are diagnostic only; no full-domain maximum or volume-L2 target
has been met.

## Value-zero local shear mode

The [OpenAI paper, Appendix B.8 and C](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
separates local shear control from cumulative-moment restoration. In
particular, Appendix C uses rapid radial variation to make an order-one
change in shear with a small profile/moment change, then restores all five
moments on a reserved interval. This is a more relevant mechanism than
fitting the mean's complete momentum with a handful of low-order bubbles;
the paper's nonaxisymmetric pulses cancel the leading annular residual.

`slope_swirl_bridge_screen.py` tests a low-frequency precursor. It adds
`64y^3(1-y)^3(y-y0)` to normalized swirl, where `y0` maps to `X=1`.
The mode leaves swirl *value* at `X=1` unchanged while changing its radial
slope and vanishes to order three at both bridge endpoints. The seven
amplitude samples show the expected shear/moment conflict. At amplitude
`-2`, the sampled angular maximum falls from `3016` to `2274` and both
reference cones pass, but the outgoing fifth moment grows by about
`1.34`, compared with only `0.0225` of downstream deletion capacity.
At amplitude `+1`, the fifth-moment change is negative and the cones
still pass, but the eight-point complete momentum maximum increases from
`6111` to `6255` and angular maximum from `3016` to `3387`.
The best of the seven under the one-sided fifth-moment capacity and two
reference-cone checks is the original zero-amplitude candidate.

This rules out **this low-order value-zero shape**, not shear modulation
in general. Its fifth-moment first variation is roughly `-0.61` per unit
amplitude at `eta=.2`; the allowed exterior capacity is about `0.022`.
Next test a compact high-frequency radial modulation whose field amplitude
falls as `1/N` while its slope remains order one, and include its growing
viscous cost in the physical residual. A sampled cone alone will not
establish the paper's all-phase admissible cone or moment restoration.

## High-frequency shear scale separation

`high_frequency_shear_screen.py` adds
`2 B(y) sin(2 pi N (y-y0))/(2 pi N)` to normalized swirl, for
`B(y)=64 y^3(1-y)^3`, `N=8,16,32`, and `y0` corresponding to `X=1`.
At that reference point the added swirl value is zero and its slope is
independent of `N`. The maximum *field* modulation falls from `0.0393`
to `0.0195` to `0.00942`; the added fifth moment relative to the same
balanced reference falls from about `3.20e-4` to `8.26e-5` to `2.07e-5`.
Thus the experiment reproduces the paper's local shear versus cumulative
moment separation numerically, with positive swirl across the sampled
radial bridge.

At `N=16`, `X=1`, the sampled relaxed-cone margins at `eta=.2,.3`
increase from about `1.57,1.02` to `21.65,10.42`. This is only a
two-point snapshot; a sine perturbation does not implement the paper's
all-phase admissible shear loop. The eight-point complete physical
momentum maximum simultaneously rises from about `6111` to `9954`,
and its angular maximum from `3016` to `9509`. The change is not a
momentum correction. This outcome is consistent with the paper's
[Section 2 and Appendix C](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf):
the mean annular residual is subsequently canceled by nonaxisymmetric
oscillatory momentum flux, followed by higher-order and mean corrections.

The next indispensable implementation is a realizable stress target and
localized nonaxisymmetric pulse, with an explicit residual budget.
Do not promote a stronger sampled cone or small moment perturbation to
`pde_validated`; the full momentum maximum is worse and no volume-L2
gate has been checked.

## Physical stress and wave preflight

The existing `radial_peak_cone.py` Kelvin/covariance source generator now
resolves the extended field's `compact.joined.inner` and its actual
`join_X=.09375`, rather than assuming the older `3/64` bridge edge. It
also accepts an explicit stress quadrature order. This lets
`extended_wave_source_screen.py` generate candidate-specific sources for
the balanced and `N=16` shear backgrounds at `X=1`, `eta=.2`,
`tau=.5*2^-5.5`. The physical target is the radial primitive of the
**complete finite-time residual**. It is distinct from the paper's
normalized leading-order target used by `extended_relaxed_cone_screen.py`.

This distinction is decisive. The normalized snapshot cone at this point
passes, but the complete-residual target has positive projection onto the
local shear direction: about `+4.833` for the balanced field and `+4.857`
for the modulated field. The existing strict physical-wave cone requires
that projection to be negative. The `N=16` source's nonnegative covariance
fit has zero local error, but its two selected pulses both use angular
mode `1`; this is not the separated two-mode construction used by the
compact curl-wave prototype. No wave has been added to the candidate.

High-frequency stress integration must resolve the carrier. At the same
modulated point, Gauss order `16` gives a spurious transverse target and
ratio; orders `32`, `64`, and `96` agree on the positive shear projection
near `+4.857` and ratio near `0.0127`. The balanced field's sign is stable
from order `16` through `64`. `extended_physical_cone_map.py` therefore
uses order `64`; none of its ten sampled nodes (`X=.75,.9,1,1.1,1.25`,
`eta=.2,.3`) passes the strict physical cone. This finite map is not a
continuous impossibility result and does not rule out a differently
constructed paper leading profile.

`radial_flux_closure_screen.py` tests an even narrower ansatz: only
`Q_{r theta}` and `Q_{r z}` carry the correction and both vanish at the
bridge edges. Cylindrical divergence then requires the weighted bridge
integrals of the complete residual to vanish. At `eta=.2`, the balanced
field gives approximately `-3.95758e-4` and `+0.0996378` for those two
integrals; at `eta=.3`, `-3.73997e-4` and `+0.143986`. Gauss orders
`16` and `32` agree to the displayed precision. These force nonzero
outer-edge fluxes in that restricted model. The denser quadrature also
finds a local complete-residual maximum around `6.7e4`, well above the
earlier eight-point maximum. The paper uses more stress components,
mean corrections, and exact five-moment restoration, so this restricted
closure failure cannot be generalized to its full construction.

The next gate is a mean-profile solve that restores all five outgoing
moments and produces the appropriate strict stress direction on a radial
interval. Only after that should the existing Kelvin/curl-wave modules be
coupled to the extended candidate. A zero-error pointwise NNLS fit or a
positive normalized snapshot cone is insufficient for a supported wave
or the requested complete-momentum threshold.

## Two-stage radial moment control and physical-cone window

`radial_moment_step.py` adds a solenoidal streamfunction step that rises
across the old bridge, stays constant past its `X=1.5` exit, and falls to
zero on `1.75<X<3`. Its perturbation is exactly absent for `X>=3`, so the
same physical velocity and pressure are recovered there. At the old exit
an amplitude `+1` changes normalized `M` by `0.9799959036379546` at
`eta=.2` and `0.9543891407651217` at `eta=.3`, agreeing with the
streamfunction prediction to roundoff. It also changes nonlinear `J` and
`S`; this is a control direction, not a repaired mean profile.

Combining shear amplitude `+2` at 16 radial cycles with moment-step
amplitude `-.5` produces a first strict **physical full-residual** cone
pass at `X=1`, `eta=.2`: `lambda_squared=823.16`, target projection
`N=-2.7204`, ratio `0.3241` (Gauss order 64). The holdout passes at
`eta=.3` and a second sampled time, but fails at `X=.95` and `1.05`.
The denser `moment_shear_band_map.json` finds passes at `X=.995,1,1.005`
for both `eta=.2,.3`, and failures at `.985,1.015` for both. At `eta=.2`
the sampled passing span corresponds to physical radii about
`0.01513364`–`0.01520950`, only `7.59e-5` wide for `nu=.01` and
`tau=.5*2^-5.5`. These finite samples neither certify a continuous
space-time cone nor give a viable wave-support width. The local full
residual is still thousands, and no nonaxisymmetric correction is present.

The exact outer-field equality is weaker than moment restoration.
Piecewise Gauss-64 evaluation after the return step, at `X=3.5`, gives
`delta[M,I,J,S,Cp]` relative to the same shear background of approximately
`[0,0,-.65084249,+.11752822,0]` at `eta=.2` and
`[0,0,-.62348152,+.26373130,0]` at `eta=.3`. A fixed-slice diagnostic
with one derivative `U` bump and three separate `E` bumps on `1.76<X<2.98`
has a full-rank initial four-row Jacobian, but a deterministic 15-start
nonlinear solve leaves maximum moment defects `0.11677` and `0.14008`;
its best corrected swirl also becomes negative. Thus that particular
compensation basis fails even before it is lifted to a divergence-free
three-dimensional field. This is a numerical failure of the tested
basis/search, not an impossibility result.

The next construction must broaden the strict physical cone while
preserving positive swirl, then solve the five outgoing moments with
smooth `eta`/time-dependent coefficients in a solenoidal representation.
Only a supported wave with controlled remainder and independently checked
complete momentum may be promoted toward the `1e-3` gates. This follows
the distinction between [Appendix C's shear-loop and five-moment repair,
and Sections 7–9's oscillatory/mean residual corrections](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

## Return-geometry and bounded-repair follow-up

The default streamfunction return is not the only smooth solenoidal
geometry. `RadialMomentStep` now also permits its falling transition to
overlap the original rise after the tested `X=1` cone point. All six
geometries in `moment_step_return_screen.json` leave that pointwise field
and its inner radial stress primitive unchanged, so both sampled axial
cone points still pass. Returning over `1.05<X<1.3` reduces the outgoing
`J` defect to about `+0.0136` at `eta=.2`, but raises `S` to `+0.8736`;
the default return over `1.75<X<3` gives `J=-0.6508`, `S=+0.1175`.
The same tradeoff appears at `eta=.3`. These are alternatives in the
tested return family, not an optimum or a general lower bound.

The separate bounded slice solve in `moment_shear_bridge_repair.py`
allows two meridional derivative bumps, including one inside the active
bridge after the cone point, and three azimuthal bumps. It constrains
the sampled swirl to remain positive. Twenty-one deterministic starts
per axial slice did not repair the four nontrivial outgoing moments:
the best maximum defects were `0.2304` and `0.2757` at `eta=.2,.3`.
The optimizer also hit its evaluation limit. This is a failed numerical
search for this five-mode basis, not evidence that all positive-swirl
repairs are impossible.

`moment_step_amplitude_refine.json` shows that, with the same 16-cycle
shear, reducing the streamfunction amplitude from `-.50` to `-.45`
already makes `lambda_squared` negative at `X=1,eta=.2`; tested amplitudes
through `-.25` also fail. Hence the observed cone pass relies on an
order-one mean change and is not the small `O(N^-1)` profile modulation
whose moment defects [Proposition C.2](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
repairs. The current evidence favors redesigning the leading mean
profile for a broader strict cone and small outgoing defects before
attempting a wave, rather than tuning this single step further.

## Constructive five-moment slice repair after the cone point

The failed four/five-bump solves hid a useful structure. With `E` fixed,
any downstream `U` correction preserving outgoing `M` and `J` has a
minimum possible `integral U^2`: project onto the span of `1` and
`H=sqrt(2X)E` on its support. `meridional_moment_lower_bound.py`
evaluates this shape-independent bound by piecewise Gauss quadrature.
At `eta=.3`, a correction restricted to `1<=X<=3` needs at least
`3.57715896` of total normalized `U^2`, versus the baseline budget
`3.55791762`; a `U`-only repair there cannot also match `S`.
At `eta=.2` the corresponding slack is positive `0.06611706`.
The bound is numerical, assumes unchanged `E` and fixed radial support,
and says nothing against joint `U,E` reconstruction.

The joint route is now executable at fixed axial slices. Five compact
azimuthal bumps on `1.005<X<2.98` restore outgoing `I` and `Cp` while
keeping sampled `E` positive. A constrained quadratic solve with a
12-dimensional smooth meridional basis supported on `1<X<3` then
restores `M` and `J`, and uses a null direction to set `S` exactly.
`coupled_five_moment_slice.json` records successful `eta=.2` and
`eta=.3` solutions. At `eta=.2`, the smaller azimuthal correction
suffices; at `eta=.3`, only the larger tested positive-swirl correction
works in this finite basis, with `S` slack about `0.00588` before the
last null-direction adjustment. The maximum five-moment defect is
below `5e-13` on the construction quadrature. Independent Gauss-128
and Gauss-160 rules in `coupled_five_moment_holdout.json` keep the
maximum below `5e-14`. The smallest sampled ratio `E_corrected/E_base`
is about `0.10965` on the finer rule at `eta=.3`.

This result closes **two fixed normalized radial slices only**. The
coefficients are not yet functions of `eta` and time; the `U` profiles
have not been lifted through a physical streamfunction, and no new
three-dimensional field or pressure has been evaluated. The smooth
`U` basis uses a narrow `0.02`-wide `X` entrance/exit, so its derivative
and viscous cost must be measured before it is accepted as a mean
correction. The local physical cone remains a narrow finite sample,
and complete momentum is still orders of magnitude above `1e-3`.
The next necessary construction is a smooth parameter-dependent
five-moment branch followed by a solenoidal physical lift and complete
residual check. [Appendix A, Lemmas A.1–A.2](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
motivates independent bump moments and small-discrepancy nonlinear
repair; our order-one correction is a numerical experiment outside
that lemma's smallness guarantee.

## Physical lift and measured narrow-taper cost

`coupled_moment_physical_lift.py` now turns the two maximum-slack slice
solutions into one time-dependent Cartesian field. It smoothly blends
their `eta` coefficients on `.2<=eta<=.3`, tapers them to zero outside
`.1<eta<.4`, and preserves the zero total `U` correction by projecting
the interpolated coefficients onto the radial mass-nullspace. The swirl
correction is axisymmetric; the meridional correction is generated by
the explicit streamfunction `psi=nu^(3/2) q^(1-A) G(X,eta)` with
`G_X=delta U`. Thus the added velocity is solenoidal by construction,
compact in the selected radial/axial band, and exactly absent for
`X<=1` and `X>=3`. The pressure is still inherited from the prior
field, so this is a diagnostic lift, not a dynamically repaired mean.

The first complete-momentum screen on 15 physical nodes gives a sampled
maximum about `1.2898e7` for the lifted field, versus `8648` for the
same base points. The largest sampled hotspot is `X=2.99,eta=.2`, in
the narrow `0.02`-wide falling transition. On that point, refining the
fourth-order spatial FD factor from `.001` to `.000125` changes the
momentum norm from `1.2898e7` to `1.2982e7`, while the divergence
error falls from `170` to `0.0442` in approximately fourth-order ratios.
The divergence behavior is consistent with the exact streamfunction
identity and under-resolution at the original spacing. The converged
large momentum is a genuine obstruction for this narrow-taper field;
pressure and time interpolation are also unfinished. No spatial
volume-L2 or `1e-3` gate is claimed.

`taper_width_capacity.py` optimizes the five azimuthal correction
coefficients for each width against the **actual twelve-mode smooth U
basis**, imposing exact sampled `I,Cp` and a positive-swirl floor.
At `eta=.3`, the best tested `S` slack is `+0.00756` for width `.02`,
then `-0.02966`, `-0.09848`, and `-0.15606` for widths `.05`, `.10`,
and `.20`. At `eta=.2`, width `.05` remains feasible, but `.10` and
`.20` fail. Therefore a wider transition that might reduce the
physical curvature is incompatible with exact five-moment matching
in this tested finite basis; this is numerical evidence about the
chosen modes and supports, not a general no-go theorem. The next
mean-profile design needs more radial freedom or a leading profile
whose moment defects are small before the stress/wave stage, as in
the paper's [Appendix A/C construction](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).

## Wider modal repair and a lower physical residual

The twelve-mode width failure was a basis limitation. With fixed
azimuthal coefficients from the width screen, increasing the smooth
meridional Legendre basis to degree 19 makes width `.05` feasible at
both fitted axial slices (`S` slacks `+0.1168` and `+0.0160`). Degree 19
also barely permits width `.10`; degree 31 permits width `.20` at these
same slices. The corresponding raw-basis conditioning rises with degree,
so none of these slice facts alone is a stable physical construction.

`wide_taper_five_moment_slice.py` constructs a width-`.05`, degree-19
two-slice repair with all five radial moments below `3e-13` in the
piecewise quadrature. Its first physical lift has a sampled complete
momentum maximum `1.1423e7` on 18 points, still dominated by the
return region. Rather than changing those moments, the constrained
optimization in `wide_taper_curvature_optimize.py` minimizes a radial
second-derivative proxy over the same exact `M,J,S` manifold while
holding `E,I,Cp` fixed. That proxy falls by factors about `452` at
`eta=.2` and `2.18` at `eta=.3`; five-moment defects remain below
`5e-13`.

The optimized physical lift reduces the **same 18-point** complete
momentum maximum from `1.1423e7` to `1.3716e6` (about `8.3` times).
At the new hotspot `X=1.01,eta=.3`, spatial FD refinement stabilizes
the norm at `1.371617e6`, while its divergence estimate falls from
`0.0209` to about `5.4e-6` as the step factor goes from `.001` to
`.000125`. Thus the reduction is not a coarse-stencil artifact, but
the remaining residual is still around nine orders of magnitude above
the `1e-3` maximum target. The pressure is still inherited rather than
reconstructed for the new mean profile; the coefficient interpolation
has not been five-moment-matched over a continuous `eta`/time range;
and no nonaxisymmetric stress correction or spatial-volume L2 gate
has been run for this diagnostic candidate.

## Pressure-only diagnostic on the curvature-optimized lift

`curvature_pressure_screen.py` keeps that velocity fixed and fits 18 compact,
axisymmetric pressure modes on 30 physical nodes at the same late time.
The sampled complete-momentum maximum changes from `1.371562e6` to
`1.371590e6`, so this pressure fit does not improve the controlling peak.
At ten disjoint spatial nodes, it worsens the maximum from `8.32194e5` to
`3.26524e7`; at ten nodes at a nearby time, from `6.41250e5` to
`2.51566e7`. The pressure derivative agrees with a direct finite-difference
residual to about `1.5` in absolute residual units on one training node,
far smaller than these changes. The fit is therefore rejected, and its
coefficients are saved only to reproduce the diagnostic.

The pressure-independent azimuthal residual alone reaches `2.85728e5` on
the training nodes, `1.71833e5` on the spatial holdout, and `1.32412e5`
on the time holdout. Thus an axisymmetric pressure correction cannot meet
the momentum gate for this velocity. The next constructive step must change
the velocity/stress design, including its azimuthal equation, before another
pressure fit is useful. These are sampled obstructions for this candidate,
not a global impossibility claim.

## Compact swirl transport response

`curvature_swirl_transport_screen.py` tests eight compact axisymmetric
swirl modes against the **pressure-independent angular equation** while
holding the meridional velocity fixed. The operator is evaluated in
physical cylindrical coordinates with fourth-order differences; one
direct Cartesian residual check per split agrees with the angular
response model within `7.2e-5` absolute units. The modes preserve
analytic divergence, but change the swirl-dependent `S` and `Cp`
moments and the radial centrifugal term.

Weak regularization (`ridge=.1`) reduces the 30-point training angular
maximum from `2.85728e5` to `2.03193e5`, but raises the spatial
holdout maximum from `1.71833e5` to `1.93281e5` and the nearby-time
holdout maximum from `1.32412e5` to `1.48940e5`. Selecting the best
ridge by the worse of the two holdout maxima chooses `ridge=10`, with
only marginal changes (`1.71812e5` spatial, `1.32396e5` temporal).
At directly re-evaluated sample nodes, the full residual is slightly
worse. Thus this fixed-scale compact swirl basis does not provide a
robust angular correction; its coefficients are a diagnostic, not a
promoted field. A dynamically evolved swirl/stress correction with
simultaneous moment and radial-momentum accounting is needed.

## Physical stress cone after the moment-preserving lift

`curvature_physical_cone_screen.py` applies the complete-residual stress
primitive to the current lifted field. The radial integral is now split
at the narrow correction edges `X=1`, `1+width`, `3-width`, and `3`;
otherwise a low-order Gauss panel can miss the new support just outside
`X=1`. At the reference time and `eta=.2,.3`, the strict physical cone
still passes exactly at `X=1`. It fails at every sampled `X=1.005`,
`1.01`, and `1.025` at both slices. The sensitive failures at
`X=1.005` and `1.01` persist with Gauss order 24 after order 12.
For example, at `X=1.005,eta=.3`, the cone ratio is `1.34` at order 24
(it must be below `1`), and at `X=1.01,eta=.3` the target projection
has the wrong positive sign (`+49.25`).

This links the momentum hotspot to a loss of the local positive-stress
representation needed by the paper's wave mechanism. A wave cannot be
attached across this sampled return region using this background and
the same strict local cone. The corrective mean profile must first
restore a cone margin over a finite region while controlling curvature;
the current two-slice moment repair alone does not do that. This is a
sampled physical analogue, not the paper's normalized continuous cone
or a theorem that no alternative field can work.

## Delaying and widening the moment-preserving U return

The narrow return starts its U correction immediately at `X=1`, where
the strict physical stress cone is positive. `delayed_taper_capacity_screen.py`
tests a later start while retaining the existing five E-bump coefficients.
At `eta=.3`, degree-19 modes cease to have positive `S` slack by
`start_X=1.01`; degree-31 modes remain feasible there, but lose slack
by `start_X=1.02` for width `.05`. With `start_X=1.005`, degree 31,
and width `.4`, the two fixed-slice slacks are `+0.11302` and
`+0.00704`. The exact M/J/S repair followed by radial-curvature
optimization leaves all five fixed-slice defects below `7e-13`.

The width-`.4` delayed physical lift is constructed from the same
streamfunction formula, with its primitive now beginning at `X=1.005`.
On the **same 20 dense physical nodes** at the reference time, the
complete-momentum maximum falls from `1.371562e6` for the previous
curvature lift to `1.035288e6` (about `24.5%`). A separate 24-node
middle/return screen peaks at `6.73344e5`; it does not reveal a larger
sampled return hotspot. At `X=1.015,eta=.3`, spatial FD refinement
stabilizes the new residual norm near `1.035321e6`, while the divergence
estimate decreases from `1.10e-2` to `2.98e-6` between step factors
`.001` and `.000125`. Analytic solenoidality follows from the
streamfunction and axisymmetric swirl construction, not from the FD
number alone.

At Gauss order 24, the strict physical cone now passes at `X=1.005`
for both fitted slices (ratios `.541` and `.934`, below `1`). It still
fails at `X=1.01` and beyond. Starting instead at `X=1.01` creates a
new `1.866e6` residual spike near `X=1.015` that a coarser sampling
missed. Widening the `X=1.005` taper from `.05` through `.1` and `.2`
to `.4` reduced the dense sampled maximum from `1.166e6` to
`1.035e6`; width `.6` rose to `1.043e6`, so `.4` is the best in this
small scan. These values are local, one-time finite-difference
diagnostics. The candidate is still far above `1e-3`, has no volume-L2
gate or continuous cone, and is not a promoted PDE solution.

An independent piecewise Gauss48 integration of the **actual lifted
velocity** (rather than its coefficient algebra) gives five-moment
defects below `4.8e-14` at `eta=.2` and `2.5e-13` at `eta=.3` at the
reference time. The sampled minimum swirl profile remains positive
(`E>0.0174`, relative to the baseline `E` above `.109`).

At a nearby time `k=5.25`, the same 18 similarity nodes give a full
momentum maximum of `1.057105e6` for the previous lift and
`7.97930e5` for the delayed width-`.4` lift (again about `24.5%`
lower). This is a time holdout, not a uniform-time bound.

`curvature_mean_patch_screen.py` separately fits twelve compact
axisymmetric vector-potential/pressure mean modes to the **previous**
curvature lift. It improves same-time train maximum `1.234e6` to
`1.088e6` and disjoint spatial maximum `1.328e6` to `1.083e6`, but
the nearby-time maximum worsens from `1.178e6` to `2.144e6`.
Its time-independent local amplitudes are therefore rejected as a
dynamic correction. The code remains a reusable exact-curl mean basis,
not an accepted lift or a substitute for the paper's time-dependent
mean/wave equations.

## Reallocating swirl capacity and screening a momentum tangent

`delayed_wide04_component_decomposition.py` isolates the entrance
hotspot at `X=1.015, eta=.3`: the stepped baseline is about `9.2e3`,
the E-only correction about `2.5e4`, and the U-only correction about
`1.04e6` in complete momentum norm. The U return therefore dominates
the present defect. `radial_step_cone_amplitude.py` shows that weakening
the radial step from `-.5` to `-.45` already loses the sampled strict
cone at `eta=.2`; that simple amplitude reduction does not preserve
the wave-admission precondition.

`delayed_e_capacity_optimize.py` reallocates the five E bump
coefficients while retaining the same fixed-slice I/Cp constraints.
For `start_X=1.005`, width `.4`, and degree 31, the two S slacks are
`+.11558` and `+.00938`. With the reoptimized E and re-solved U,
the same 20-node entrance maximum falls from `1.035288e6` to
`9.565792e5`. A separate 24-node middle/return maximum is
`6.288559e5`; the 18-node nearby-time maximum is `7.372641e5`.
Independent integration of the actual physical lift gives fixed-slice
five-moment defects below `1.9e-13`. This improves sampled momentum
without providing a continuous stress cone or a global PDE bound.

`delayed_momentum_tangent_screen.py` then varies only the `eta=.3` U
profile in three projected directions tangent to the M/J/S constraints,
repairs the quadratic S term exactly, and selects by the full Cartesian
residual. The materialized experimental candidate is
`delayed005_wide04_momentum_tangent.json`. Its entrance maximum is
`9.481077e5` on the same 20 nodes, down about `0.9%` from the
reoptimized-E field. On disjoint spatial nodes the maximum changes
from `9.766330e5` to `9.528168e5`; on disjoint nearby-time nodes,
from `8.650167e5` to `8.573629e5`. Independent physical five-moment
integration remains below `2.4e-13` at the two reference slices.

The same tangent worsens the 24-node middle/return maximum from
`6.288559e5` to `7.571412e5`. On the 18-node `k=5.25` screen the
overall maximum improves only from `7.372641e5` to `7.307352e5`,
while the far return point `X=2.975, eta=.3` worsens from about
`1.03e5` to `1.91e5`. This makes the candidate a useful diagnostic,
not a promoted field: fixed-time slice moments do not control the
dynamic far-return residual. A subsequent correction must include
entrance and return nodes at several times in one objective, and must
restore a finite positive cone region before invoking oscillatory
stress cancellation. No volume-L2 or `1e-3` claim follows from these
local screens.

## Two-time entrance/return objective

`delayed_multitime_momentum_screen.py` reuses the three projected
fixed-M/J/S directions, but fits all three Cartesian residual components
at ten `eta=.3` nodes at each of two times (`k=5.5` and `5.25`). The
nodes include the entrance, middle, and far return at `X=2.975`; the
far return is normalized by its own baseline residual in the surrogate
so it is not silently sacrificed to the much larger entrance. The
candidate criterion requires lower sampled overall maximum, lower
middle/return maximum, and no far-return increase at either time.

The first larger surrogate step is infeasible on the exact S manifold
or raises the entrance maximum. A small step, `0.025` times the
linearized fit, produces the experimental
`delayed005_wide04_multitime_tangent.json`. On the 20 fit nodes, its
overall maximum changes from `9.565792e5` to `9.553478e5`, its RMS
from `5.208011e5` to `5.197927e5`, and its middle/return maximum
from `6.288559e5` to `6.255864e5`. The `k=5.25` maximum changes from
`7.372641e5` to `7.363150e5`. Its reference-time two-slice physical
five-moment defects remain below `1.8e-13` by independent quadrature.

`delayed_multitime_holdout.py` uses unseen `eta=.29,.31`, different
spatial nodes, and two other nearby times (`k=5.4,5.15`). The entrance
maxima improve by about `0.043%` at both times; far-return maxima
improve by about `2.4%`. The middle maxima **worsen by about `2.8%`**
at both times. The small fit-point gain does not survive as an
all-region improvement. The three fixed-E, single-slice U directions
are therefore insufficient for promoting a corrected field. The next
constructive route must give the correction a coupled time/axial shape
and optimize the physical stress cone and complete momentum jointly,
then measure a spatial volume norm rather than extrapolating from
these selected points.

## Axial support width is not a free cure

`CoupledMomentPhysicalLift` now accepts explicit `axial_rise_start` and
`axial_fall_end` in a slice artifact; absent those keys, the previous
`(.1,.4)` support is unchanged. The same septic interpolation and
streamfunction formula preserve analytic divergence freedom, and the
profiles at `eta=.2,.3` stay exactly fixed when the end changes.
`delayed_axial_support_screen.py` compares fall ends `.4`, `.45`, and
`.5` on three radial nodes, six axial coordinates, and two times.

Widening the fall end does **not** improve the complete-momentum maximum:
at the reference time it is `9.566187e5`, `9.572940e5`, and
`9.574215e5`, respectively. More importantly, the residual at
`X=1.015,eta=.4` rises from about `1.02e4` to `2.20e5` and `5.44e5`.
The same ordering persists at the nearby time. The wider taper simply
moves a large correction into an axial region that previously matched
the base field. No widened-support candidate is retained. This screen
also shows why checking only the two fitted `eta` slices would miss an
important spatial hotspot.

The cone diagnosis remains more fundamental. The strict physical cone
uses a radial primitive of the complete tangential/axial residual.
At `X=1.01`, the reoptimized-E field still fails the sampled cone at
both reference slices; at `eta=.3` its target projection along the
local shear normal is positive (`+4.44`), violating the required
negative sign. The component decomposition attributes the large
entrance residual chiefly to the U streamfunction return, whose radial
curvature is amplified by the onset. A useful next experiment is a
broad, endpoint-flat U return basis with five-moment projection,
screened jointly for cone margin and full momentum across entrance,
middle, and far return. It must be judged on off-slice and nearby-time
points before any wave or PDE claim.

The obvious lower-degree alternative also has a measured capacity
barrier. `delayed_e_capacity_optimize.py` now accepts a U degree and
reoptimizes the five E bumps under the same I/Cp equalities and E floor
for each degree. At width `.4`, start `X=1.005`, the `eta=.3` finite-basis
S slack is `-.17371`, `-.05876`, `-.01882`, and `-.00356` for degrees
11, 19, 23, and 27; it becomes `+.00938` only at degree 31. Direct
five-moment construction also fails for degree 11 and for degree 19
at `eta=.3`. These are multistart finite-basis numerical results, not
a proof of global infeasibility, but they explain why simply dropping
the high-degree modes cannot retain the current five-moment repair.
The radial step/return basis or the baseline moments must change to
obtain both low curvature and sufficient S capacity.

## Alternate radial basis and swirl-floor capacity

`delayed_bump_basis_capacity.py` checks overlapping compact C4 radial U
bumps of widths `.3`, `.5`, `.8` and counts 16, 24, 32, both alone and
combined with the first twelve global tapered Legendre modes. With the
reoptimized E profile fixed, every tested basis has negative S slack
at `eta=.3`. The best hybrid in this grid still has slack `-.03210`
(width `.5`, 32 local bumps) and an original-basis condition above
`3e5`; the pure local basis is further from feasibility. Hence a
localized smooth return cannot simply replace the global degree-31
space under the present moment debt.

`radial_restore_capacity.py` moves the outer streamfunction restoration
start from `X=1.55` through `2.4`, keeping its entrance rise and final
zero trace at `X=3`. The degree-27 `eta=.3` S slack remains near
`-.0035` throughout; degree 31 remains near `+.0094`. This shows that
the outer restoration location, by itself, does not remove the
capacity bottleneck for the fixed E profile.

The positive swirl floor was then varied in the E optimizer. At a
relative E floor of `.05`, degree 27 becomes only barely feasible:
`eta=.3` S slack is `+.000223`; at `.01` it is `+.00256`. Degree 19
remains infeasible even at `.01` (`-.05304`). The exact five-moment
degree-27, floor-`.05` candidate was physically lifted and screened.
Its 20-node complete-momentum maximum is `1.010374e6`, **worse** than
the degree-31 floor-`.11` value `9.565792e5`. The sampled cone at
`X=1.01` still fails both slices: the `eta=.2` ratio is `1.410`, and
the `eta=.3` target-normal projection is positive (`+3.329`). The
weaker swirl therefore buys too little smoothness and loses physical
momentum performance. This field is not promoted. The next structural
degree of freedom should change the entrance rise of the radial mean
step and co-design its five-moment return, rather than further lowering
the swirl floor or moving only the far restoration.

## Changing the entrance rise of the radial mean step

`RadialMomentStep` now accepts an optional `rise_end`; the default
retains the original end `X=1.5`. `radial_rise_capacity.py` varies it
while leaving the step amplitude `-.5`, outer restoration, E profile,
and delayed U support fixed. Moving the rise end to `1.3` gives much
more degree-27 S slack at `eta=.3` (`+.0560`) but fails the sampled
strict cone at `X=1`. The original `1.5` passes that point but has
degree-27 slack `-.00356`. Fine scans find a narrow crossover: at
`rise_end=1.46`, degree-27 slack is `+.000176` and both `X=1` cone
rows pass at Gauss24. A Gauss64 repeat gives cone ratios `.611` and
`.839` there; these are necessary pointwise conditions only.

The degree-27 physical lift at `rise_end=1.46` is exactly repaired at
the two reference slices but raises the 20-node entrance momentum
maximum to `1.018855e6`. Retaining degree 31 at the same rise end
gives a larger S slack (`+.01318` at `eta=.3`) and a better physical
candidate, `delayed005_rise146_degree31.json`. Its independent
piecewise physical-velocity integration finds five-moment defects
below `3.7e-13` at both slices and minimum relative swirl above `.11`.

Against the previous rise-`1.5`, degree-31 reoptimized-E field,
the new lift improves complete Cartesian momentum on three separate
local screens:

| Screen | Previous max | New max |
| --- | ---: | ---: |
| 20 entrance nodes, `k=5.5` | `9.565792e5` | `8.595515e5` |
| 24 middle/return nodes, `k=5.5` | `6.288559e5` | `5.727185e5` |
| 18 nodes, `k=5.25` | `7.372641e5` | `6.624820e5` |

At disjoint `eta=.29,.31` spatial nodes and two other nearby times,
the entrance maximum is about `10.5%` lower, the middle maximum about
`21.5%` lower, and the far-return maximum about `9.8%` lower at **both**
times. This is the first rise-shape change here that improves all three
sampled regions and the nearby-time holdout together. The screens are
still sparse and do not establish a spatial volume norm or a global
maximum.

The cone remains the next obstruction. At `X=1.005`, Gauss64 gives
strict ratios `.784` and `.989` at `eta=.2,.3`; the latter converges
to `.989130` at orders 96 and 128. Gauss24 had incorrectly reported
`1.119` for that sensitive node, so the higher-order result is used.
At `X=1.01`, both slices fail, and the `eta=.3` target-normal
projection is positive. The new rise shape therefore gives a real
sampled momentum improvement and a tiny positive onset cone margin,
but no certified cone region spanning the return, supported wave,
volume-L2 gate, or `1e-3` PDE acceptance. The field remains an
experimental candidate, not a promoted solution.

## Delaying the return onset and testing a dynamic pressure patch

With the improved `rise_end=1.46` step, moving the U correction start
directly from `X=1.005` to `1.01` loses fixed-slice feasibility at
`eta=.3`: even after E reoptimization under the same `.11` relative
swirl floor, the degree-31 S slack is `-.002443`. An intermediate
`start_X=1.008` is feasible (S slack `+.003822`) but its 20-node
complete-momentum maximum increases from `8.595515e5` to
`1.162895e6`. At `X=1.01,eta=.3`, the sampled cone target-normal
projection becomes negative, but the strict ratio is `2.31` and still
fails. Moving the onset alone trades one cone defect for a much larger
momentum peak; this candidate is rejected.

`delayed_similarity_pressure_screen.py` tests a different mechanism:
six compact pressure-only modes in the physical similarity coordinates,
scaled by `q^(-2A)`. Because they leave velocity untouched, they also
leave analytic solenoidality and all five moments untouched. Their
pressure gradients are computed analytically; comparison with the
independent finite-difference full residual differs by at most `11.1`
against residuals of order `1e6`. A two-time fit on the retained
rise-`1.46` velocity lowers training RMS from `3.92766e5` to
`3.83769e5`, but its maximum barely changes (`8.48799e5` to
`8.48321e5`). On disjoint space/time nodes it creates a severe new
maximum (`7.73588e5` to `1.38712e6`), especially near
`X=1.07,eta=.22`. The pressure patch is rejected. This result only
rules out the tested six-mode fit; it does not show that all matched
pressure corrections are ineffective. The present dominant momentum
error still calls for a dynamically coupled velocity/stress correction.

## Similarity-scaled exact-curl velocity screen

`delayed_similarity_curl_screen.py` adds three compact poloidal
streamfunction modes in the physical similarity coordinates. Their curl is
analytically divergence-free and the modes follow the changing similarity
scale. The unconstrained two-time full-momentum fit lowers the 12-node
training maximum from `848798.65` to `121859.83` and a separate space/time
holdout maximum from `773588.41` to `391616.60`. These are local screens,
not the volume-L2 or global maximum gates. More importantly, independent
physical-velocity integration finds third/fourth outgoing moment defects
`(-.009736,+.013798)` at `eta=.2` and `(-.040835,+.054180)` at `eta=.3`.
The unconstrained correction therefore cannot be accepted as a matched
mean field. See `delayed_similarity_curl_screen.json` and
`delayed_similarity_curl_moments.json`.

`delayed_similarity_curl_constrained.py` tests six modes, using two axial
supports, while enforcing the third and fourth moments at both reference
slices. Compact streamfunction support preserves the first moment, and the
swirl profile is unchanged so the second and fifth moments stay fixed.
Independent order-64 integration in
`delayed_similarity_curl_constrained_moments.json` confirms all five
reference-slice defects below `1.8e-13`. Yet the constrained peak-focused
fit changes the training maximum from `848798.65` to `859415.44` and the
holdout maximum only from `773588.41` to `764132.93`. The RMS-focused fit
also raises the training maximum (`883449.49`). The strong unconstrained
momentum reduction is therefore mostly incompatible with exact two-slice
moment closure within these six modes. This is a limitation of the tested
basis and fit, not a general obstruction theorem.

At the intervening `eta=.25`, the constrained patch alters the fourth
moment by only `5.6e-5`, while the full field's fourth-moment defect is
`-.01170`; the baseline already has essentially that defect. A viable
next mean-flow solve must impose moment and cone conditions across an axial
interval, not only at `eta=.2,.3`, while coupling velocity and pressure
or stress corrections. Neither curl candidate has passed a continuous
strict cone, global momentum maximum, or volume-L2 test; both remain
experimental (`accepted:false`).

## A third slice does not close the axial interval

`delayed_midband_moment_patch.py` adds three compact physical
streamfunction-curl modes and two axisymmetric swirl modes, all supported
between `eta=.2` and `.3`. Fitting the previously unconstrained middle
slice at `eta=.25` drives its second through fifth outgoing moment
defects from approximately `(0,-.000775,-.011756,-.0000185)` to
roundoff; the first moment remains near `1.2e-8` on the order-48
quadrature and is unchanged by the compact curl. The swirl stays above
the `.11` relative floor in the sampled rows. The endpoint slices are
unchanged because the axial bump vanishes there.

The interval is still open: at `eta=.225` the fourth-moment defect grows
from `+.027638` to `+.028585`, and at `eta=.275` it changes from
`-.024703` to `-.024002`. On nine disjoint space/time momentum nodes,
the maximum is essentially unchanged (`684006.34` to `684001.36`) and
RMS rises (`310962.23` to `311412.87`). The defects change sign along
the axial band, whereas the one-bump patch has a fixed axial profile.
This result in `delayed_midband_moment_patch.json` rules out promoting
the three-slice fit as interval closure. The next construction needs
axially varying correction coefficients (or a coupled mean/stress solve),
with physical curl derivatives and momentum included in the fit.

## Three interior axial collocation slices

`delayed_axial_lagrange_patch.py` tests axial Lagrange profiles over
`eta=.2` to `.3`, with the curl derivative included in the physical
velocity. Each profile isolates one interior slice (`.225,.25,.275`)
while vanishing at the two reference endpoints. Twenty compact poloidal
radial modes and five swirl modes supply the finite-dimensional freedom.
At `.225`, the original three- and seven-mode poloidal bases have
negative fixed-slice S capacity (`-.01987` and `-.01337` with the small
two-mode swirl repair). Twenty modes reduce the capacity gap to about
`-.00496` under a generic E optimization. Optimizing E for the actual
twenty-mode U basis changes it to a positive `+.00526`; the other two
interior slices have positive capacity too. Thus the earlier optimizer
singularity was not a general impossibility of this compact construction.

The resulting *physical* field satisfies the second through fifth
outgoing moments at all three fitted slices to numerical precision.
The first-moment differences there are between roughly `2e-8` and
`2.3e-7` on the order-48 quadrature; each added streamfunction is
analytically compact and has zero radial mass. Nevertheless, the
interpolated field overshoots between nodes: at `eta=.2625` the fourth
moment defect is `+.031916` versus the base `-.025999`, and the sampled
relative swirl minimum falls to `.108897`, below the `.11` floor. On
fifteen independent space/time momentum nodes the maximum barely changes
(`684006.34` to `683102.24`) while RMS rises from `305561.19` to
`333164.21`. The full record is
`delayed_axial_lagrange_patch.json`.

Three exact collocation slices therefore do not imply interval matching.
The next solve must include interleaved eta nodes and the swirl floor
directly, with axial regularity and complete momentum in the objective.
This candidate is rejected (`accepted:false`); no continuous cone,
volume-L2 gate, or `1e-3` PDE gate is established.

## Joint axial-interval moment/floor fit

`delayed_axial_interval_fit.py` replaces the three isolated slice solves
with one least-squares fit of the same compact physical modes. Seven eta
positions share all coefficients; each contributes the second through
fifth moments and a positive-swirl-floor penalty. Analytic moment
Jacobians make this a reproducible finite-dimensional screen. The fit
reaches the 200-evaluation limit rather than satisfying all conditions.
For example, at `eta=.2625` it cuts the fourth-moment defect from
`-.025999` to `+.000118`, but at `eta=.2875` the defect remains
`-.013222`; the third-moment defect there is `-.004158`. The sampled
relative swirl dips to `.1100066` at `eta=.2375`, below the `.11005`
fit margin (though above the older `.11` floor).

Five interleaved eta holdouts in `delayed_axial_interval_fit.json` show
partial generalization, not closure: at `eta=.26875` the fourth-moment
defect falls from `-.027034` to `+.007599`, while at `eta=.29375` it
stays near `-.006957`. On the independent full-momentum space/time
holdout the maximum increases from `684006.34` to `684699.04` and RMS
from `305561.19` to `370825.54`. This candidate remains rejected.

The [OpenAI paper's physical description](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
explicitly assigns cancellation of the annular singular momentum
residual to nonaxisymmetric oscillatory momentum flux, followed by mean
and higher-order corrections; matching the radial moments is needed to
keep the exterior unchanged. Our mean-only interval fit is not that
mechanism. The next constructive stage must jointly seek an interval of
strict *physical* stress-cone feasibility and a supported wave/mean
correction, rather than infer PDE improvement from moment fitting alone.

## Pressure-only physical cone window is geometrically blocked

`delayed_pressure_cone_window.py` applies the six compact similarity
pressure modes to the retained rise-`1.46` physical field. Since pressure
does not change velocity or shear, the two strict-cone inequalities at
each point are linear in the six coefficients. At `tau=.5*2^-5.5`, the
eight-point window `X=1.008,.012,.016,.02` and `eta=.2,.3` is linearly
infeasible for cone ratio at most `.8` with stress margin `.1`. The
two-eta radial prefix through `X=1.012` is feasible, but adding
`X=1.016` makes the system infeasible. Individual-node checks show that
the last two radii fail at both eta values; this is not merely conflict
between six pressure basis functions.

`delayed_pressure_cone_geometry.py` removes even the finite-basis
restriction and lets the *axial* stress component vary arbitrarily while
holding tangential stress and shear fixed. At `X=1.02` it still cannot
enter the cone at either eta, even with ratio bound `.999` and zero
stress margin. At `X=1.016`, arbitrary axial stress could reach the
near-boundary `.999` cone but not the more robust `.8` cone. Letting
both stress components vary does make the pointwise wedge accessible:
the smallest two-component stress shifts at `X=1.02` have norms about
`5.30` (`eta=.2`) and `58.26` (`eta=.3`) under the `.8`/`.1` margins.
These are necessary target changes, not realized wave covariances or
mean-flow corrections. The raw geometry and quadrature results are in
`delayed_pressure_cone_window.json` and
`delayed_pressure_cone_geometry.json`.

Pressure-only tuning therefore cannot provide a supported physical-cone
window through `X=1.02` on this fixed velocity. The next candidate must
change the tangential stress target and/or the underlying shear, while
preserving the moment and exterior conditions, then verify a continuous
space-time cone before invoking the paper's oscillatory realization.

## Compact solenoidal velocity response of the physical cone

`delayed_swirl_cone_response.py` tests one axisymmetric toroidal mode at
a time over three compact radial supports, with a common axial envelope.
Its full Cartesian momentum and radial stress primitive are quadratic
in the mode amplitude; the scan also enforces the sampled `.11` relative
swirl floor. The best short support `(1.005,1.04)` opens only one of eight
nodes, `X=1.008,eta=.2`, at amplitude `.11`. Its cone ratio is `.999715`
and `lambda_squared` becomes negative at the other three radii on both
eta slices. Wider supports open no nodes. This is not a robust wave
window and the selected amplitudes create outgoing moment defects.

`delayed_coupled_cone_response.py` adds a compact poloidal streamfunction
mode and retains its full quadratic cross-advection with the swirl mode.
Among 6561 sampled amplitude pairs, the best eight-node choice
`(swirl,poloidal)=(.125,-.0325)` opens two nonadjacent points on the
`eta=.3` slice, while all four `eta=.2` points fail. Its sampled maximum
momentum residual falls from `868264` to `670805`, but the outgoing
third/fourth moments shift by roughly `(-.00775,+.01439)` at `.2` and
`(-.00631,+.02304)` at `.3`. Even allowing each eta slice to choose its
own pair within the scanned box opens at most two of four radii per
slice. The best pair with positive `lambda_squared` at all eight nodes
still has a negative minimum cone margin (`-56.90`). See
`delayed_swirl_cone_response.json` and
`delayed_coupled_cone_response.json`.

The two-mode screen confirms that changing tangential and meridional
velocity can improve local momentum and alter the cone, but a single
radial shape per component does not supply a contiguous strict region.
The next construction needs multiple radial and axial controls optimized
against the cone *and* outgoing moments, followed by an independent
space-time and wave-support check. Neither screened field is promoted.

`delayed_coupled_cone_holdout.py` performs that disjoint space/time
momentum check on the two-mode point with two sampled cone passes. The
complete Cartesian maximum rises from `759459.69` to `2012315.55`, and
RMS from `422178.21` to `1203253.64`. The largest new errors occur at
intermediate `eta=.22,.28`, where the common axial envelope changes
rapidly. Finite-difference divergence reaches `.545` on this stencil;
the toroidal plus streamfunction formula is analytically solenoidal.
`delayed_coupled_divergence_refine.json` halves the spatial step at the
worst residual node (`X=1.025,eta=.22`): divergence falls approximately
by sixteen each time (`.23164,.01452,.000908,.0000568`), as expected
for fourth-order truncation, while the residual norm stabilizes near
`2012359.96`. The large momentum defect is therefore real, not a
finite-difference artifact. The local eight-node momentum
reduction does not generalize and cannot justify a PDE claim. See
`delayed_coupled_cone_holdout.json`.

## Eight-mode joint cone, moment, and momentum screen

`delayed_multimode_cone_fit.py` uses two overlapping axial bands and two
radial shapes each for swirl and poloidal velocity (eight exact-curl or
axisymmetric-toroidal modes total). The full physical stress primitive
and Cartesian momentum include every quadratic cross-advection term.
Its 12-node training grid is `X=1.008,.012,.016,.02` and
`eta=.2,.25,.3`. The frozen quadratic response and three-slice moment
bases are saved in `delayed_multimode_cone_model.npz`, so optimization
weights can change without recomputing the expensive physical stencil.

A first search, preserved as
`delayed_multimode_cone_fit_unregularized.json`, lowers the sampled
maximum residual from `868264` to `345583` but raises the independent
nearby-time maximum from `759460` to `1010482`. Its moment defects and
physical cone remain unacceptable. Adding a precomputed nine-node
nearby-time full-momentum response to the objective improves the
selected local candidate: training maximum `207310`, nearby-time
maximum essentially unchanged from the base (`759460`), and nearby-time
RMS `365894` versus `422178` for the base. The quadratic response
agrees with direct finite differences on that holdout to `2.1e-7` in
Cartesian residual components. The optimizer reaches its iteration
limit; these are sampled improvements, not a converged solution.

Only three of twelve training cone points pass. The worst margin is
`-19.46`, and the outgoing third/fourth moment defects at `eta=.3`
are about `(-.02170,+.02087)`; the field cannot be promoted.
`delayed_multimode_pressure_admission.py` tests whether compact pressure
could finish its cone window without changing velocity. Arbitrary axial
stress can enter the `.8`/`.1` cone at nine of twelve points and the
near-boundary `.999` cone at ten, but `X=1.016,.02` at `eta=.2` fail
even the near-boundary free-axial test. The six pressure modes are
linearly infeasible across all twelve points. The next mean correction
must change tangential target or shear at those two radii, while
restoring moments and maintaining the improved nearby-time momentum.

## Disjoint momentum and eight-mode tradeoff

`delayed_multimode_disjoint_holdout.py` checks the selected eight-mode
velocity at nine points excluded from both its cone fit and nearby-time
objective: `X=1.011,1.017,1.023`, `eta=.215,.265,.315`, and
`tau=.5*2^-5.35`. Its full Cartesian momentum maximum drops from
`745020` to `466919`, and RMS from `415572` to `276076`. This establishes
that the local momentum gain extends to one additional sampled window;
it does not imply a uniform or volume-integrated bound. The sampled
fourth-order finite-difference divergence is at most `.1162` here;
analytic toroidal and curl perturbations are solenoidal.

`delayed_multimode_moment_feasibility.py` uses the frozen eight-mode
response for twelve bounded local least-squares starts on the outgoing
moments at `eta=.2,.25,.3`. The best scaled L2 defect is `.04869`, with
largest raw moment defect `.000294`. This is a numerical feasibility
screen, not a proof of exact closure or global infeasibility. Crucially,
its training momentum maximum rises to `2208277` and none of the twelve
cone nodes pass. Moment repair in this family cannot be chosen in
isolation from momentum and cone geometry.

`delayed_multimode_tradeoff_fit.py` therefore runs three joint cached
searches with stronger moment penalties. At weight `.05`, a second
Pareto candidate has six of twelve cone nodes passing (versus three for
the earlier selection), a largest raw moment defect `.00606` (versus
`.02170`), and nearby-time maximum `628508` (versus `759460`). Its
training maximum `559320` is worse than the earlier candidate's
`207310`, but below the base `868264`. The passing points form the
`X=1.008,.012,.016` rows at `eta=.2,.25`; no `eta=.3` row passes.
On the same disjoint nine-point check, this candidate lowers the base
maximum `745020` to `579141` and RMS `415572` to `339029`, but is worse
than the previous candidate there. Both candidates remain available;
neither is promoted to a PDE solution.

For the second candidate, even freely changing axial stress cannot
enter the `.8`/`.1` cone at `X=1.02` for any of the three sampled eta
values; it can at the other nine points. More radial velocity control
near this outer part of the transition is required before a pressure or
oscillatory correction can close this window. Any new search must also
repair the moments across an interval, retain independent space-time
momentum improvement, and verify a continuous strict cone rather than
only the twelve fit nodes. The two disjoint momentum reports and three
moment/tradeoff reports are reproducible JSON artifacts in
`experiments/root_st073/`.
