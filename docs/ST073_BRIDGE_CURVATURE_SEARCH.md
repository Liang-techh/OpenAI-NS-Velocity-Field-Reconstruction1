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
