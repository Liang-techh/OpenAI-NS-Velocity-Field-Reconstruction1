# ST073 continuation: local kernel to matched scale recursion

The current user objective requires nonzero divergence-free finite-energy 3D
velocity, shrinking/elongating/winding cores, dynamically matched inner,
transition and outer regions, and complete physical momentum max and spatial
volume L2 below 1e-3 on an explicit domain/forcing contract. Python/MATLAB and
interactive delivery remain required. This is not the old visualization-only goal.

Imported frozen user delivery under experiments/root_st073; all 173 manifest
payloads verified and original ZIP SHA retained. No bundled sync script executed.
Current main integration at 6a293b3c does not contain this complete bundle; work
is isolated on codex/st073-transition-next without disturbing prior local edits.

Host progress: archived <f16 NPZ arrays cannot load on this Windows NumPy.
host_interface.py regenerates float64 interface traces from unchanged ST073-V
model at k=0,3,6 (17 eta nodes each). Maximum computed complete momentum residual
is 5.315e-7; pressure-gradient difference from supplied JSON <=2.701e-13.
This is same-model replay, NOT a new independent operator or whole-domain test.
Original files stay byte-identical. Windows longdouble epsilon is2.22e-16.

Next construction tasks:
1. Use actual corrected velocity/pressure jets at X=1/64, |eta|<=.5; rebuild
   transition compatibility rather than reusing ST068 leading-only moments.
2. Derive moving-interface normal velocity, mass flux and stress from source
   coordinate map, including time-dependent boundary motion; include axial caps.
3. Construct a solenoidal transition with matched velocity and stress; evaluate
   its full residual, not only continuity or divergence. No zero-extension claim.
4. Match an outer finite-energy solution with explicit forcing; do not select
   arbitrary residual-canceling force. Add nonaxisymmetric corrections only with
   an explicit closure/realizability check.
5. Distinguish prescribed k-dependent sampling from dynamically generated scale
   recursion; test material winding and shrinking/aspect metrics over the full
   registered window. No extension toward tau=0 without new evidence.

Local scope remains X<=1/64, |eta|<=.5, tau in[.5/64,.5], nu=.01 unforced.
No global finite-energy field, exterior match or scale-recursion acceptance yet.

## Moving interface mass budget completed

`experiments/root_st073/moving_interface.py` integrates the curved side and both
flat caps at k=0,3,6 using orders12 and24. Side inflow balances cap outflow;
net fluid flux at order24 is below2.8e-20. The relative outward flux is NOT zero:
it equals minus the shrinking domain volume rate. At k6 it is1.702657768e-5.
Order12/24 results agree; this is numerical integral consistency, not a rigorous
uniform bound or momentum acceptance. Interface labels move at b_r=-r/(2tau),
b_z=-(.5-h)z/tau. Treating this interface as a material no-through-flow surface
would contradict the frozen kernel. The transition must transport the measured
side/cap flow and match momentum stress; a closed impermeable shell is unsuitable.
Report: `experiments/root_st073/moving_interface/mass_flux.json`.

## Stress and moving momentum interface data

`interface_stress.py` constructs the physical Cartesian gradient directly from
frozen radial coefficient jets, sigma=-pI+nu(grad u+grad u^T), and moving flux
u((u-b).n)-sigma n. Both side and caps are included, orders12/24, k0/3/6.
Portable nodewise arrays include geometry, gradients, stress and oriented flux.
At k6/order24, total axial momentum outward flux=2.515468846e-7 and angular
momentum outward flux=1.080789491e-8. These are integrated fluxes, NOT residual
norms or admission thresholds. Gradient trace max5.68e-14, curl/archived-evaluator
vorticity replay difference1.14e-13; same-model algebra consistency only.
Next use these stresses together with mass/velocity traces to build a transition;
check volume momentum-rate plus boundary flux before selecting an exterior.
Artifacts: `experiments/root_st073/moving_interface/stress_flux.json` and
`stress_k0.npz`, `stress_k3.npz`, `stress_k6.npz`.

## Radial continuation constructed

An explicitly separate order10 continuation keeps ST073-V axis data and nu,
extends Xmax from1/64 to3/64 (radius factor sqrt3), and preserves solenoidality
through the same full recurrence. Frozen model unchanged. At the original
interface k6, order8/10 velocity difference4.37e-13 and pressure6.02e-12.
The new annulus X in[1/64,3/64], |eta|<=.5 has k6 sampled boundary max5.748e-5;
12x18 volume quadrature max5.463e-5 and physical L2=4.668e-9, volume1.780e-7.
Orders8x12 and12x18 were evaluated at k0/3/6. This is a finite local extension,
NOT decay to an exterior or proof across all points/times. Tiny volume explains
part of L2. Independent FD and separate directional convergence remain pending.
At X=1/16 (twice original radius), order10 boundary max1.0296e-3 fails the gate;
recorded wider-radius failures must remain. Next use the new outer boundary
for a dynamical transition, or assess higher-order radial continuation with
independent derivatives; do not mask outer errors by cutting the field off.
Reproduce: radial_continuation.py then annulus_check.py in experiments/root_st073.
Artifacts and explicit new model: experiments/root_st073/radial_continuation/.

## Independent annulus derivative check

annulus_fd.py reuses the independent Cartesian/time operator, not recurrence
residual assembly, at six off-calibration points (X=.026,.044; eta=-.22,.37;
k=.4,2.7,5.5; nonzero azimuth). Spatial and temporal steps varied separately.
All FD momentum maxima remain below1.031e-5, with finest late value1.021e-5,
well below1e-3 at these points. This supports using the extended boundary as a
matching datum; it does not prove a uniform bound or an exterior solution.
Artifact: radial_continuation/independent_fd.json under experiments/root_st073.
The next unresolved construction is still dynamical outer decay/axial closure,
not further repetition of these local residual checks.

## Compact solenoidal localization control rejected

compact_control.py constructs velocity by tapering the actual poloidal
streamfunction (degree9 C4 cutoff) and swirl, not multiplying velocity alone.
It is bounded and compactly supported at every registered time, so finite
energy follows; analytic solenoidality follows from the streamfunction.
Only X<=1/64, |eta|<=.3 is preserved (NOT the full frozen axial domain).
Pressure is explicitly tapered. No residual-canceling force is added.
Independent FD at radial/axial/corner points gives late k5.5 momentum norms
1.104e5,4.045e4,2.111e5, stable after halving steps. Divergence FD defects
shrink about16x, consistent with fourth-order truncation; they are not blanket
numerical divergence acceptance. This localization is rejected as an NS field.
The failure quantifies the missing transition dynamics. Do not repeat simple
cutoff sweeps or present compactness as matching. Next solve a stress/flux-driven
transition correction or a decaying exterior with evolving interface data.
Artifact: experiments/root_st073/compact_control/report.json.

## Updated source priority

Read docs/ST073_PAPER_ROUTE.md before new construction. Official PDF was reopened;
actual heat-exterior direct-join screen rejects value-only matching.
Evidence: experiments/root_st073/heat_join/screen.json.

## C2 radial swirl component constructed

swirl_bridge.py uses actual corrected coefficient jets at X=3/64 and the frozen
heat-exterior amplitude to build an autonomous quintic on r in[ri,2ri]. Values,
first and second physical radial derivatives match at both ends. At nine k/eta
probes endpoint absolute errors are <=4.22e-15,3.09e-12,3.05e-9 respectively;
101 radial samples per bridge retain positive swirl. This removes the tangential
viscous traction jump at those interfaces without fitting force. The polynomial
is OUR interpolation, not a formula claimed from the paper. It is only a swirl
component, not the complete transition: no poloidal velocity, pressure, dynamic
residual, axial caps or finite-global-energy acceptance is supplied by this step.
Next construct compatible poloidal streamfunction transport and derive the full
transition residual stress; do not interpret endpoint matching as NS acceptance.
Artifact: experiments/root_st073/swirl_bridge/coefficients.json.

## Poloidal streamfunction bridge constructed

poloidal_bridge.py builds a septic psi(r,z,t) matching corrected inner radial
jets through order3, and zero jets at outer radius2ri. Combined u_r=-psi_z/r,
u_z=psi_r/r is solenoidal when moving endpoints are differentiated consistently.
This is an autonomous interpolation, not a paper-derived dynamic solution.
Its zero outer streamfunction forces return transport: at k6,eta=.4 core axial
flux5.22367e-5 is balanced by annular flux-5.22367e-5. Sampled annular axial
velocity ranges -2.3740 to2.1197; return flow cannot be ignored in stress fitting.
Endpoint residuals are recorded; full z/time derivatives, pressure joining,
complete momentum and axial closure remain pending. No arbitrary force added.
Next evaluate the combined swirl/streamfunction bridge and integrate its actual
stress defect, then test the paper-inspired moment/covariance correction.
Artifact: experiments/root_st073/poloidal_bridge/coefficients.json.

## Callable full radial joining field

joined_field.py now evaluates all velocity components and pressure across the
inner field, Hermite transition and heat swirl exterior inside |eta|<=.5.
Poloidal psi_z analytically differentiates coefficient jets AND moving endpoints;
pressure matches value and first two radial derivatives. No forcing is fitted.
Independent Cartesian FD at four transition points, two times and two steps:
early momentum norms38.88--62.31, late7851.6--12562.9, stable under refinement.
Divergence errors decrease about16x on step halving (late finest<=1.054e-6).
This is a callable divergence-structured background with a large dynamic defect,
NOT a completed transition, axial closure, finite-energy global field or theorem.
Next compute actual integrated tangential residual stresses of THIS background
and test realizability/correction following paper Sections4,7-9; do not transfer
leading-only cone values or call smooth endpoint joins a momentum solution.
Artifact: experiments/root_st073/joined_field/report.json.

## Actual transition residual moment targets

transition_stress.py integrates the full FD tangential residual at fixed physical
z,t. With zero inner stress, the radial inverses give nonzero terminal stresses:
at late k5.5 rtheta is about-1.52 to-1.66 and rz ranges-2.24 to2.19 across three
axial locations. Thus this restricted inverse cannot vanish at both interfaces
without moment repair. Orders8/12 quadrature are recorded, not a continuum bound.
The sampled maximum reaches6.05e4 on the denser radial nodes; previous four-point
checks were not maxima over the whole transition. Both records are retained.
The pointwise PSD trace lower bound2*sqrt(T_rtheta²+T_rz²) is3.32--5.42, only an
algebraic bound. PSD completion also contributes diagonal/radial momentum and
is NOT a wave realization. No arbitrary forcing or post-hoc stress cancellation
has been counted as a solution. Next solve these measured moment defects with
boundary-jet-preserving background corrections before wave/cone admission.
Artifact: experiments/root_st073/transition_stress/moments.json.


## Boundary-preserving swirl moment repair: rejected as a PDE improvement

Implemented `swirl_moment_repair.py`: an axisymmetric pure-swirl bubble with
cubic zeros at both interfaces, preserving velocity and first/second jets.
Pressure, poloidal velocity and frozen core stay fixed. This is an autonomous
basis choice inspired by the moment-repair route, not a formula from the paper.
Two bounded coefficients fit at k=3, eta=-.3,0,.3 are approximately
[2.94846774, 1.25961520]. Independent holdouts use k=.7,5.5 and eta=+/-.2;
training quadrature uses 8 points, reporting uses 12. No force is introduced.
At late holdouts terminal angular stress magnitude falls from about1.60 to
.0105--.0123 (over99% reduction), but full sampled momentum grows from
4.36e4--4.55e4 to7.02e4--7.08e4. Axial moment remains unchanged.
Thus integral compatibility alone is not a residual reduction. Do not promote
this optional candidate to the default JoinedField or count it as PDE progress.
The final held-out point is also evaluated with half the spatial FD step.
Artifact: experiments/root_st073/swirl_moment_repair/report.json.

Next implementation tasks:
- Replace moment-only fitting with a constrained angular PDE collocation solve:
  use multiple radial cubic-endpoint bubble modes and axial modes, preserve
  interface jets, and minimize the full theta residual while bounding moment
  defects. Freeze core data and retain all rejected baselines.
- Include multiple training scales and disjoint scale/axial holdouts; one frozen
  two-parameter fit does not remove scale-dependent moment defects exactly.
- Couple poloidal streamfunction modes to axial moment and axial PDE repair;
  swirl alone cannot change the current axial residual at fixed poloidal field.
- Recompute radial pressure compatibility after velocity corrections. Angular
  residual reduction alone does not control centrifugal/radial momentum.
- Only then evaluate actual correction-stress compatibility and realizability;
  retain the separate axial closure, finite-energy and scale-recursion work.
Full momentum max and physical-volume L2 gates remain1e-3; none has passed.


## Angular PDE collocation: component gain, no full-field gain

`angular_collocation.py` adds 12 endpoint-preserving swirl modes (four radial
Legendre modes, three axial polynomial modes). It trains at k=1,4 and
eta=-.3,0,.3, using normalized pointwise theta residual plus a soft moment
penalty, bounded coefficients and mild regularization. Base field evaluations
are cached without changing numerical differentiation or the fitted fields.
Disjoint holdouts k=2.5,5.5, eta=+/-.2 use18 radial nodes. Theta maxima fall
about13%, but full momentum maxima rise about0.35--0.38%; terminal angular
stress magnitudes also rise about0.3%. Candidate is NOT promoted to baseline.
At late eta=.2, full max45530.7485 agrees with spatial-step-halved45530.7480.
The coarse divergence estimate8.78e-4 drops to1.79e-6 on refinement; stencils
near the C2 join affect coarse finite-difference errors. No supremum or volume
L2 certificate is claimed. Axial/poloidal correction is now the priority.
Artifact: experiments/root_st073/angular_collocation/report.json.


## Coupled poloidal / pressure collocation

`poloidal_collocation.py` fits six streamfunction and six pressure bubbles to
FULL nonlinear momentum at k=1,4 and eta=-.3,0,.3. Quartic streamfunction
endpoint zeros preserve velocity and all spatial jets through order2; physical
psi_z includes moving-interface and source-to-physical derivative factors.
The compact correction adds zero net axial flux through each annular slice.
It is analytically divergence-free; finite-difference divergence is separately
reported. Pressure bubbles have cubic endpoint zeros. Frozen core, original
swirl and heat exterior are retained, and no force is fitted.

`affine_momentum.py` precomputes affine velocity/derivative jets; evaluating
advection retains the exact quadratic coefficient interactions. The resulting
surrogate differs from a fresh full-field Cartesian FD evaluation by at most
2.50e-8 at the checked training slice. This is a discretization consistency
check, not independent proof of the PDE. The bounded nonlinear fit converges
in10 evaluations. These are autonomous basis/optimizer choices, not a direct
implementation of the paper's oscillatory stress realization.

At disjoint k=2.5,5.5 and eta=+/-.2 holdouts, full sampled maxima decrease
about8.8--9.1%. Late eta=.2 falls45373.40 ->41259.86; halving spatial FD step
gives41259.8573, divergence2.96e-6. Training decreases12--29%.
Angular moment compatibility slightly worsens; this is an improvement candidate,
NOT an accepted NS solution. Finite global energy, axial closure, scale
recursion, non-axisymmetric stress realization and global gates remain open.
Artifacts: experiments/root_st073/poloidal_collocation/report.json;
physical-volume subdomain audit: poloidal_collocation/volume_audit.json.

Derivative convention: q_z and eta_z from source_coordinates must be divided
by sqrt(nu). For y=r/ri-1, y_z=-(1+y)q_z/(2q). For increasing physical time,
q_t=-q_tau and eta_t=-eta_tau; the full FD evaluator already implements this
minus sign. Do not confuse remaining time tau with physical time t.


Physical-volume audit over eta in[-.4,.4], ri<r<2ri, all azimuths, uses
12 radial by8 axial Gauss nodes and the exact coordinate volume Jacobian.
At k2.5, sampled maximum2989.51 ->2606.71 and volume L2 5.48547 ->5.12251.
At k5.5, sampled maximum68136.50 ->59387.68 and volume L2 26.3981 ->24.6533.
These are approximately12.8% maximum and6.6% L2 improvements on this subdomain;
BOTH remain far above1e-3. Quadrature convergence and continuum bounds are
not established. These denser axial checks exceed the earlier slice maxima,
which must not be reported as whole-transition maxima.
Next combine poloidal, pressure and angular modes in one full-momentum solve,
with actual moment penalties and scale-dependent modes; preserve disjoint
holdouts and volume metrics before considering any candidate promotion.


## Joint full-momentum and moment fit: tradeoff, not a new accepted baseline

`joint_collocation.py` composes the12 poloidal/pressure and12 swirl modes,
retaining nonlinear cross interactions in the affine-jet momentum evaluator.
It uses an analytic optimizer Jacobian, bounded coefficients, six training
slices at k=1,4 / eta=-.3,0,.3, and soft weight3 penalties on normalized
r-squared angular and r-weighted axial residual moments. Core and interface
jets stay fixed; no external force is fitted. Source files include direct
full-field FD holdouts and the physical-volume audit for reproducibility.

At disjoint k=2.5,5.5 / eta=+/-.2, angular terminal stress magnitude improves
about16--17% against the original JoinedField; full maxima improve about9%.
Against the previous poloidal/pressure candidate, the additional slice maximum
improvement is only about0.07--0.10%. Axial moment changes are mixed by parity.
On the volume subdomain at k5.5, joint max59244.28 is below prior59387.68,
but L2 24.7024 is ABOVE prior24.6533. Thus the joint candidate does not dominate
the previous candidate on the user's two gates. Preserve both, promote neither.
Original baseline on the same quadrature: max68136.50, L2 26.3981.
All remain far above1e-3 and no finite-energy whole-space field exists yet.
Artifact: experiments/root_st073/joint_collocation/report.json.

Next avoid endless fixed-coefficient tuning: introduce explicit scale-dependent
streamfunction/swirl/pressure amplitudes, differentiating those amplitudes in
space and physical time, then test a cross-scale solve. Inspect actual moment
transport and outer matching freedoms before imposing a radial stress inverse.
If gains stay small, advance the paper-inspired stress realization route rather
than equating increasingly complex Hermite bubbles with dynamical matching.


## Explicit log-time amplitude comparison

`scale_collocation.py` extends24 joint coefficients to a0+s*a1, where
s=(k-3)/3 and k=-log2(2*tau). The amplitude depends only on time, so it
preserves the spatial divergence-free streamfunction/swirl construction and
all endpoint jets. The optimizer explicitly includes s_t*delta_u in momentum,
s_t=1/(3*ln(2)*tau), for increasing physical time. Product-rule jets are
compared against full-field FD with varying amplitudes at every time stencil.

A24-coefficient constant model and48-coefficient varying model are fitted on
IDENTICAL k=1,3,5 and eta=-.3,0,.3 data, with the same pointwise/moment objective,
bounds and regularizer. Initializing the varying model from the fitted constant
model avoids attributing extra training data to amplitude variation. Intercepts
and slopes are individually bounded[-4,4]; effective coefficients can exceed
that interval, which is reported explicitly. Holdouts remain k=2.5,5.5 and
eta=+/-.2, with a separate physical-volume subdomain audit. This is a finite
log-time ansatz, NOT a closed scale-recursion law or singularity construction.
Artifact: experiments/root_st073/scale_collocation/report.json.


Matched-data outcome: at k5.5, varying coefficients slightly worsen subdomain
max59244.12 ->59266.70 while slightly reducing L2 24.7044 ->24.6954.
At k2.5, max2600.444 ->2600.241 while L2 5.133039 ->5.133393 worsens.
No consistent improvement on both gates; do not promote this variant or
continue increasing temporal polynomial degree without a new mechanism.
The next route is the actual shear/phase/amplitude dynamics of paper Section7,
with the approximation limits recorded in ST073_PAPER_ROUTE.md. Global energy,
axial closure, recursive scale control and the1e-3 gates remain unfulfilled.


## Frozen phase/amplitude implementation

`frozen_pulse.py` extracts physical angular velocity F=u_theta/r and radial
tangential shear g=(d_r u_theta-F,d_r u_z) from the actual JoinedField and
poloidal/pressure candidate. It implements a frozen cylindrical local amplitude
ODE with the pressure projection needed to preserve n dot a=0, and viscosity.
The frozen phase has constant tangential wavevector and n_r_dot=-g dot n_tan.
This is an autonomous physical-coordinate approximation inspired by paper
Section7, not substitution of full ST073 into the leading normalized theorem.
Radial strain, axial gradients and background variation along trajectories are
omitted; angular mode1 is not a high-frequency justification.

Viscous damping is factored exactly as exp(-nu*(|n0|^2*t+(n0 dot n_dot)*t^2+
|n_dot|^2*t^3/3)) before integrating the undamped matrix ODE, avoiding relative
transversality errors once damped solutions fall below absolute solver tolerance.
Each sampled time uses the matrix singular value for optimal amplification;
the covariance uses the initial seed that maximizes the sampled peak gain.
Nonnegative covariance fitting is recorded separately and is NOT PDE admission.
Artifact: experiments/root_st073/frozen_pulse/report.json.


Screen result:156/192 sampled points have positive inviscid frozen growth
parameter, but only8/192 pass the diagnostic stress-direction inequalities.
The late candidate point with largest growth rate has the WRONG target-dot-N
sign. Its45 phase choices (angular modes1,2,4; five axial ratios; three radial
ratios) achieve peak optimal gain1.49859 over0.1*tau with viscosity retained.
The initial restricted radial-direction family missed this transient growth;
no absence-of-growth claim is retained. Nonnegative time-averaged covariance
fits the two-component target algebraically, but this cannot certify compact
pulse construction, mean cancellation, or full residual reduction.

Next integrate rays/amplitudes at points that actually satisfy the directional
screen, compare with full local velocity-gradient evolution, and construct a
supported vector-potential field only after checking its own full residual.
Keep compact moment tails and global axial/energy closure as separate failures.


## Full-gradient rays and actual collar residence

`affine_pulse.py` uses the full physical Cartesian velocity gradient J at the
two late poloidal/pressure points passing the previous directional screen.
It solves n_dot=-J^T*n and the incompressible projected amplitude equation,
retaining viscosity through an integrated scalar damping factor. The45 initial
wavevectors per point include angular-like modes1,4,16,64,256 and axial/radial
ratios-1,0,1. None of these sampled frozen-affine evolutions amplifies above its
initial norm over0.1*tau. Transversality relative errors stay below1.84e-9.
This is not a no-growth theorem for the changing background or other phases.
The average covariance fits algebraically with nonnegative weights, again
showing why an algebraic fit alone is not a dynamical construction.

The points sit only1.788e-5 physical length from the moving inner interface.
The dimensional cutoff rate nu/distance^2 is3.128e7, about3.57--3.59e5 times
the reduced inviscid growth rate; this is a scaling warning, not an operator
bound. A0.1*tau frozen horizon also predicts radial travel about30 times that
distance, making that local approximation unsuitable for a supported pulse.

`ray_residence.py` therefore integrates particles in the ACTUAL time-dependent
candidate. Both cross the moving inner interface in about3.84e-5 time units,
only0.00348*tau, far shorter than the assumed0.1*tau pulse interval. The event
uses q(z,tau), not a fixed cylindrical boundary. The y=.02 outer probe was
explicitly defined and is not claimed to be a cone boundary; neither path
reached it. Artifacts: affine_pulse/report.json and affine_pulse/residence.json.

Next redesign the transition to admit an interior region with appropriate
stress direction and transport residence BEFORE adding localized oscillations.
Do not attach the earlier short-time1.50-gain wave to this narrow collar: it
came from a different point that fails the direction screen. Preserve the
core, moment accounting, full residual gates and global-closure requirements.


## Variable-width transition construction

`width_field.py` generalizes the radial join to ro=ratio*ri while preserving
the same ST073 inner velocity/pressure jets and fixed heat amplitude. It
rebuilds all Hermite polynomials using width=(ratio-1)*ri. For
Y=(r-ri)/width, Y_z=-(1+(ratio-1)*Y)*ri_z/width; width_z/width=ri_z/ri.
These terms are retained in analytic psi_z, so changing width does not replace
the solenoidal streamfunction construction with a velocity cutoff.
The ratio2 evaluator is compared directly with the original JoinedField.

`width_screen.py` compares ratios2,3,4,6 at k=2.5,5.5 using16 radial and6 axial
Gauss nodes over eta in[-.4,.4]. It records physical-volume L2 on each ACTUAL
annulus, full sampled maxima, divergence, velocity magnitude and the same
approximate stress-direction screen. Larger annuli have larger volumes;
no volume-normalized score or smaller evaluation region hides that cost.
No fitting is performed, so this is an exploratory geometry comparison.
Artifact: experiments/root_st073/width_screen/report.json.


Width comparison result: at late k5.5, ratio2 -> ratio6 lowers sampled max
66741.68 ->4514.92 (93.2%) and physical-volume L2 26.3981 ->5.7862 (78.1%),
INCLUDING the larger volume. At k2.5 max2928.42 ->200.61 and L2 5.48547 ->1.21796.
Ratio6 is an exploratory candidate, not accepted by either1e-3 gate.
The directional screen passes14/96 nodes instead of4/96; some new passes are
near the OUTER edge, where compact-stress moment tails are still unresolved.
No continuous interval of admissibility or pulse residence is certified.

Use `WidthField(6)` / width_screen/candidate.json for the next bounded trial.
Do NOT blindly reuse old AngularModes or PoloidalModes: their normalized
radius and derivative formulas assume width=ri. Generalize them to width=5ri,
then recompute actual moments and transport residence with the new background.
The fixed heat amplitude, original core and original ratio2 implementation
remain available for matched comparisons. Larger width alone does not solve
axial closure, finite global energy, pulse realization or recursive control.


## Width-aware coupled correction basis

`wide_modes.py` generalizes the existing24 swirl, streamfunction and pressure
modes using y=(r-ri)/((ratio-1)*ri). Streamfunction radial derivatives divide
by the ACTUAL width, and y_z=-(1+(ratio-1)*y)*q_z/(2*q*(ratio-1)) uses physical
q_z. Endpoint zeros and analytic divergence are preserved. Original narrow
modules remain unchanged. Nonzero-coefficient ratio2 compatibility differs
by at most3.1e-16 in velocity and zero in pressure on the recorded four points.

`wide_collocation.py` fits the generalized basis on WidthField(6), starting
from zero coefficients rather than transferring the narrow-layer optimum.
Full nonlinear momentum and actual tangential moments enter the objective;
training k=1,4 / eta=-.3,0,.3 and holdouts k=2.5,5.5 / eta=+/-.2 stay distinct.
All residual, moment and volume quadrature radii now cover ri<r<6ri; using
old ri<r<2ri sample helpers would incorrectly omit most of the new transition.
Artifacts: wide_collocation/report.json and wide_collocation/compatibility.json.


Wide coupled outcome: at late k5.5 on MATCHED12x8 volume nodes, sampled max
4629.00 ->4261.34 (7.94%) and physical-volume L2 5.78633 ->4.89297 (15.44%).
At k2.5, max205.326 ->189.200 and L2 1.21799 ->1.02850. Earlier16x6 width
screen maxima differ due to node placement; compare paired values, not maxima
from different grids. Disjoint axial slice checks show about2% full maximum
reduction and42% angular terminal-stress reduction, with modest axial-moment
improvement. Neither moment is zero, and neither1e-3 global gate is passed.
This is a useful wide-background candidate, not an accepted global NS field.
Next recompute the directional screen and actual transport residence for THIS
corrected wider field; earlier narrow-field cone or trajectory results do not
transfer. Keep explicit stress tails, axial closure and finite-energy tasks.


## Corrected wide-background pulse geometry

`wide_pulse_screen.py` recomputes full residual prefix stresses, local shear and
the diagnostic direction conditions on the fitted Width6 field. It uses20
radial nodes at k=2.5,5.5 and eta=-.3,0,.3; the selected late passing point
maximizes physical distance to either radial interface. No old narrow-field
stress direction or trajectory is transferred. Actual particle transport uses
the changing background and moving boundaries, stopping at a radial boundary,
|eta|=.48 probe, or0.1*tau. The initial directional condition is NOT asserted
along the trajectory. A separate frozen full-gradient amplitude calculation
uses the observed residence duration, with27 explicitly recorded initial
wavevectors. This remains preparatory geometry/dynamics, not a supported
non-axisymmetric field or a full momentum cancellation.
Artifact: experiments/root_st073/wide_pulse_screen/report.json.


Wide pulse screen result:14/120 nodes pass the approximate direction test.
The selected point lies at eta=0,y=.956117 (near the outer boundary), with
physical distance7.0616e-4 to that boundary. Its actual particle reaches the
moving outer interface after0.0718045*tau, about20 times the previously tested
narrow-collar residence, although these are different selected positions.
All27 sampled frozen-affine perturbations have peak gain<=1 up to rounding.
This does not certify absence of growth along the changing trajectory.
The target radial stress remains nonzero near the outer boundary, so moment
tails and compact support are still unresolved even where direction passes.
Next evolve phase/amplitude against J(t) along the actual trajectory and
re-evaluate stress-direction compatibility there before constructing a
supported vector-potential wave. Do not treat longer residence as PDE progress
or as a substitute for reducing the remaining full residual and global energy.


## Evolving phase along the corrected wide trajectory

`evolving_pulse.py` evaluates the actual Cartesian velocity gradient at25 saved
particle positions, using fourth-order spatial differences, then integrates
the Kelvin wavevector and pressure-projected viscous amplitude with a smooth
interpolated J(t). It compares27 initial wavevectors with a fixed-start J on the
SAME observed residence time. In both models the sampled peak gain is1.0;
the best evolving final gain is0.97759. The evolving gradient trace is at most
1.22e-7 and wavevector/amplitude orthogonality errors stay below5.23e-11.
This rules out gain for this explicitly sampled local family and horizon only;
it is not an impossibility statement for all waves or a global pulse.

At the seven late points passing the approximate direction test, a mode-one
azimuthal wave has minimum local damping nu/r^2 larger than the reduced
inviscid growth rate: the ratios range about1.30--12.1. This dimensional
comparison explains the observed damping but is not a spectral bound for the
changing, inhomogeneous field. The positive-direction region still occurs
only in narrow inner/outer collars. The next constructive step should change
mean swirl/shear/pressure in the wide transition to create an interior cone
margin with growth that overcomes viscosity, while retaining the full residual
and volume scores. Then re-solve transported phase and stress tails.
Artifact: experiments/root_st073/evolving_pulse/report.json.


## Interior cone/viscosity fit on width6: rejected

`cone_fit.py` uses the exact nonlinear momentum residual of the24-mode wide
field and radial-prefix stress from it. At interior radial Gauss nodes, soft
penalties favor the approximate stress direction and inviscid growth exceeding
the mode-one local viscosity rate nu/r^2. It keeps the existing residual and
moment terms, training k=1,4,5.5 and eta=-.3,0,.3. The diagnostic geometry is
in PHYSICAL cylindrical coordinates and must not be equated with the paper's
normalized leading-cone theorem.

This fit creates three interior passing nodes at late eta=0, whereas the prior
had zero, but creates none at eta=+/-.3. On held-out late physical-volume nodes,
maximum momentum grows4261.34 ->6218.60 and L2 grows4.89297 ->5.39810;
angular terminal mismatch grows from about.307 to1.60--1.63 at eta=+/-.2.
It is rejected as a full-field candidate. The simple cone penalty cannot
replace a radial/axial profile satisfying BOTH stress direction and moment
identities. The current poloidal/pressure basis has axial modes1 and eta only;
add an even eta^2 mode before assessing whether the side failures are structural.
Artifact: experiments/root_st073/cone_fit/report.json.


## Even axial mode extension for the wide background

`even_modes.py` adds eta^2 terms to both streamfunction and pressure radial
bubbles, giving30 coefficients in all (18 poloidal/pressure plus12 swirl).
The physical z derivative of eta^m includes m*(eta/.3)^(m-1)*eta_z/.3;
this is required for exact streamfunction divergence cancellation. With the
new six coefficients zero, the extended evaluator agrees exactly with the
old24-mode evaluator at three nonzero-coefficient sample points.

`even_cone_fit.py` repeats the bounded full-momentum/moment and approximate
interior direction/growth fit with that expanded parity basis. It starts from
the previously accepted width6 fit, with all new eta^2 coefficients zero.
The aim is to check whether a symmetric axial correction can supply interior
pulse conditions at BOTH eta=+.3 and eta=-.3 without sacrificing full momentum
and physical-volume L2. Fitting penalties remain diagnostics, not acceptance.
Artifacts: even_cone_fit/report.json and even_cone_fit/compatibility.json.


Even-mode result: the new eta^2 terms leave the former field EXACTLY unchanged
when set to zero. The cone-penalized fit increases late midplane interior
passing nodes from0 to4/10, but STILL gives0/10 at eta=+/-.3. Its held-out
late slice maxima improve about5%, yet the late matched-volume maximum grows
4261.34 ->4965.24 and physical-volume L2 grows4.89297 ->5.70224.
Angular terminal mismatch at eta=+/-.2 grows about.307 ->1.72--1.85.
Therefore it is rejected. A low-degree axial polynomial did not solve the
side stress orientation and made the user's actual gates worse.

Next use the paper's radial shear modulation idea as a separate controlled
candidate: alter local shear while keeping mean velocities and radial moments
nearly fixed, then compute the FULL viscous residual and pulse compatibility.
High radial frequency can greatly increase viscosity; test this cost explicitly
before claiming any stress cone progress. Do not replace the width6 best
full-momentum candidate with either cone-penalized field.


## Radial shear loop and coupled mean repair: rejected as full-field candidates

`radial_shear_loop.py` adds a compact two/four-cycle axisymmetric swirl loop
with quartic endpoint zeros, keeping exact divergence freedom and C2 interface
jets. At k5.5, a two-cycle amplitude-.3 creates one approximate direction+mode1
growth passing node at y=.729 on BOTH eta=+/-.3 slices; the existing width6
mean field had none there. The late side sampled maximum grows only about0.2%,
and angular terminal stress falls about22%. A four-cycle -.3 loop creates
passing nodes too, but its momentum maximum grows markedly from viscosity.
These are finite-node diagnostics, not a wave or a strict paper cone certificate.

The two-cycle candidate fails the actual physical-volume gate comparison:
late max4261.34 ->4294.04 and L2 4.89297 ->5.90882; early L2
1.02850 ->1.23678. `loop_repair_fit.py` refits all24 mean coefficients in the
presence of the same loop. It recovers part of the L2 cost (late5.45039) but
remains worse than the width6 baseline, and the two side passing nodes disappear
on direct re-evaluation. Thus shear orientation, moment compatibility and full
PDE residual must be solved together. Neither raw nor refitted loop replaces
the accepted best sampled baseline; no supported non-axisymmetric pulse exists.
Artifacts: radial_shear_loop/report.json, candidate_audit.json, and
loop_repair_fit/report.json.

Next work should follow the paper's linked moment-preserving profile continuation
and oscillatory stress program. A sinusoidal shear loop by itself is too weak a
substitute. Preserve failed coefficient sets and full-volume metrics for
comparison; do not count approximate cone passes as an NS residual reduction.

## Spatially compact solenoidal closure: finite-slab construction only

`compact_potential.py` puts the current width6 field into an axisymmetric
vector potential: A_theta=psi/r and A_z=-integral_0^r u_theta(s,z,t)ds.
The new velocity is curl(chi A), where chi is a smooth physical-r/physical-z
cutoff. On the cutoff plateau it replays the entire width6 velocity and pressure;
outside the finite support it vanishes. Thus it is nonzero, spatially compact
and divergence-free by construction for every registered time, and its spatial
kinetic energy is finite at each such time. This is only defined for
tau in [.5/64,.5]; neither uniform energy control nor continuation to tau=0
is established. The radial joins inherited from the width6 field are C2.

`compact_potential_audit.py` checks two times. Four plateau samples match the
width6 field exactly. A coarse 6x6 physical cylindrical quadrature gives
kinetic-energy estimates 6.64e-5 at k=.4 and 1.18e-5 at k=5.5; these are
neither convergence checks nor uniform bounds. Finite differences at four
cutoff-collar points give momentum-residual norms 343--1674 at k=.4 and
5.73e4--3.22e5 at k=5.5. Sampled divergence is <=2.85e-6, consistent with
the curl identity and numerical differentiation. These collar defects exceed
the 1e-3 target by many orders and are not global max/L2 measurements.
The pressure is explicitly tapered with chi; no cancelling force was added.
The outer closure therefore solves the energy/support property on the finite
slab but worsens the full NS momentum gate. Next solve the cutoff-induced
momentum and pressure balance together with the mean/pulse stress construction,
then audit physical-volume norms and critical-time scaling. Report:
`experiments/root_st073/compact_potential/report.json`.

## Weighted residual moments at the compact closure

`compact_moment_audit.py` integrates the *full physical* unforced momentum
residual at k=5.5 on two fixed-z slices. The radial pieces are split at the
inner interface, width6 outer interface, radial cutoff plateau, and support
edge. These are diagnostics analogous to the weighted radial primitives in
paper Section 8, not that section's normalized auxiliary-mean defects.

At z=0, the 12-node-per-piece moments are integral r^2 R_theta dr =
+0.00143895 and integral r R_z dr = +0.000913628. The outer radial cutoff
alone contributes +0.00158829 to the angular moment; the width6 transition
contributes -0.000149333. The inner core is below 2e-6 pointwise in the
sampled residual, so the outer cutoff dominates this slice's angular defect.

At z=0.00394501, inside the axial cutoff collar, the total axial moment is
-0.525289. The width6 transition contributes -0.553372 and the core adds
+0.0264002. This is a distinct axial closure problem, not just a radial
edge-layer defect. Six-versus-twelve-node moment changes are <=2.4e-5 for
these four totals; this is limited quadrature evidence, not a norm bound.

Section 8 requires the relevant weighted source moments to be repaired before
a compact radial stress primitive can cancel the mean residual. A direct
primitive of the current defects would leak beyond its annulus. Next construct
moment-preserving velocity and pressure corrections, and an actual realizable
non-axisymmetric stress if needed, while keeping full momentum max/L2 as the
acceptance gates. Artifact: `experiments/root_st073/compact_potential/moment_audit.json`.

## Conservative balance of the two weighted defects

`compact_balance_audit.py` evaluates the exact axisymmetric conservative
identities behind the preceding residual moments. For the compact field,
Mtheta=int r^2 utheta dr and Fz=int r uz dr (the latter is analytically zero
because uz=r^-1 d_r(chi psi)). The angular balance is
int r^2 Rtheta dr = d_t Mtheta + d_z int r^2 uz utheta dr
                   - nu d_zz Mtheta.
The axial balance is
int r Rz dr = d_t Fz + d_z int r(uz^2+p)dr - nu d_zz Fz.
All derivatives use physical time t=T-tau and physical z.

At k=5.5, z=0, angular time, axial transport and axial viscosity contribute
+0.00107881, -0.000052400 and +0.000411138. Their sum 0.00143755 differs
from direct integrated angular residual by 1.40e-6. The axial kinetic and
pressure flux derivatives contribute +0.000728821 and +0.000184808; their sum
agrees with the direct axial moment to 2.23e-10.

In the axial cutoff collar z=0.00394501, the axial kinetic-flux derivative is
-0.0867909 and the axial pressure-flux derivative is -0.438498. Their sum is
-0.525289109, within 3.28e-9 of the independently integrated axial residual.
The pressure term supplies ~83.5 percent of this signed moment. Fz is ~1e-18
numerically, so its time/viscous contributions vanish at this resolution.
This identifies the main scalar closure defect as a pressure/kinetic axial flux
imbalance, not merely a large pointwise derivative of the radial cutoff.

A moment-only pressure adjustment could cancel the integrated axial defect,
but may create a large radial pressure-gradient error. The next candidate must
reconstruct pressure jointly with radial momentum and velocity/stress fluxes;
it is admissible only if full pointwise and physical-volume L2 momentum improve.
Artifact: `experiments/root_st073/compact_potential/balance_audit.json`.

## Pressure-only axial-moment repair: constructed and rejected

`pressure_moment_repair.py` adds a compact pressure bubble outside the preserved
inner core. At each fixed (z,tau), it subtracts the full integrated axial flux
H=int r(uz^2+p)dr times a radial beta bubble normalized to int r b dr=1.
At strength1 this makes H exactly zero up to quadrature; velocity and analytic
divergence are unchanged. This is an actual correction candidate, not a fitted
force, and it directly tests whether the large axial moment could be the only
remaining obstruction.

At k=5.5, axial collar z=0.00394501, the repaired H is 6.78e-20 and direct
int r Rz dr is -3.76e-9, down from -0.525289. Yet on the 5x5 physical-volume
sample grid over the compact support, full momentum maximum moves
189000.81 -> 189001.45 and sampled L2 moves 209.05091 -> 209.10133.
A least-squares strength -3.05459 lowers that grid L2 only to 208.98472 and
moves opposite to the moment correction. Neither is an accepted field. These
coarse grid numbers are candidate screens, not full-domain maxima or certified
volume norms. The baseline worst residual and the angular component remain
orders above1e-3.

Thus cancelling one weighted pressure/kinetic flux identity is insufficient.
The next correction must jointly handle radial and angular momentum, the
axial pressure profile and velocity/stress covariance, with the full physical
residual retained after each step as in the paper's Sections8-9. Artifact:
`experiments/root_st073/compact_potential/pressure_repair.json`.

## Joint axial-collar velocity/pressure candidate and scale screen

`joint_collar_fit.py` adds six compact axisymmetric modes in the axial collar:
two streamfunction modes (hence exactly solenoidal velocity), two regular
swirl modes, and two pressure modes, each with even/odd axial parity. All vanish
on the preserved z plateau and outside the compact support. A linearized full
Cartesian momentum fit at k=5.5/Gauss5 supplies coefficients; the subsequent
screen retains the exact quadratic velocity term. Strength1 overfits: training
L2 falls209.05 ->141.29, while independent Gauss6 L2 rises314.77 ->1432.30.
At strength.1, training L2 is197.95 and independent Gauss6 L2 is215.36,
with independent max181939 ->91679. This is a bounded candidate, not a new
accepted NS field. Another independent Gauss8 late grid gives max321680 ->
210291 and L2 383.97 ->288.38. The strong grid-order dependence prevents any
claim of converged physical-volume L2.

The constant coefficient set fails at earlier k=.4 (Gauss6 max946 ->2436,
L2 22.56 ->58.82). `joint_collar_scale.py` therefore tests the SAME fitted
coefficients with amplitude multiplier (tau_ref/tau)^1.5 and strength.1. Full
time derivatives of that multiplier are included; k6 uses fourth-order forward
time differences at the registered lower-tau endpoint. On independent Gauss6
samples, baseline -> candidate (max; L2) is:
  k=.4: 946.16 ->932.34; 22.56 ->22.29
  k=4: 38740 ->31923; 144.97 ->125.12
  k=5.5: 181939 ->91715; 314.77 ->215.33
  k=6: 304685 ->127480; 407.60 ->262.11.
This is genuine sampled improvement across the registered finite window, but
still misses both1e-3 momentum gates by many orders and establishes nothing
as tau tends to zero. The angular maximum remains unchanged at these grid
nodes; axial maxima sometimes increase. The sampled quadrature is not a norm
certificate. `ScaledJointCollarField` is a callable full velocity/pressure
candidate. An independent fourth-order Cartesian check at three late collar
points found momentum norms124396 ->122894, 95717 ->95622, and55512 ->55501;
divergence remained at the ~3e-6 finite-difference level.

Next fit across several times and radial/axial quadrature orders, target the
remaining angular channel and interface/cap stresses, then evaluate full
physical-volume gates. Paper Sections7-9 still require a realizable
non-axisymmetric pulse/covariance and iterative mean repair; these six slow
axisymmetric modes do not implement that program. Artifacts:
`experiments/root_st073/compact_potential/joint_collar_fit.json`,
`joint_collar_scale.json`, and `joint_collar_direct.json`.

## Radial-collar angular correction added to the joint field

The remaining angular peak at k6/Gauss8 lies at r=0.02154 and |z|=0.002264,
inside the radial cutoff collar but below the axial cutoff onset. Thus the
previous axial-collar modes are exactly zero there. `radial_swirl_fit.py` adds
two C4 compact pure-swirl radial bubbles, with the same (tau_ref/tau)^1.5 time
weight and the existing smooth axial taper. Their velocity is analytically
solenoidal. Fitting their linearized FULL momentum response at k6/Gauss8 and
retaining the quadratic self-transport gives coefficients7.12829 and6.07674.
`RadialSwirlRepairedField` exposes the complete sum as a callable field.

At fitted k6/Gauss8, angular max drops179634 ->29815 and sampled physical-
volume L2 drops342.08 ->194.37. Total maximum remains276850 because a radial
peak outside these bubbles dominates. On an independent k6/Gauss6 grid, total
max falls127480 ->94315 and L2 262.11 ->191.48. The same coefficients also
lower sampled L2 at k5.5/Gauss6 from215.33 to161.13, at k4 from125.12 to
115.55, and at k=.4 from22.29 to22.17; total maxima at these earlier grids
are unchanged because other components control them. These are finite samples
with strong grid-order dependence, not converged maximum or volume norms.

`radial_swirl_direct.py` independently evaluates the full fourth-order
Cartesian momentum at six physical collar points at k5.5. At the midplane
collar center, the norm falls111121 ->8486; the largest post-correction norm
among these six points is11773, versus111927 before. Sampled divergence stays
around1e-10. This establishes a real local angular improvement, while leaving
large radial/axial errors elsewhere and the 1e-3 target far unmet.

The next work is a coupled correction for the *radial* maximum and interface
stress, together with non-axisymmetric realizability/moment constraints from
the paper. Do not promote this finite-slab screen to an accepted NS solution.
Artifacts: `experiments/root_st073/compact_potential/radial_swirl_fit.json`
and `radial_swirl_direct.json`.

## Late radial-peak incremental joint repair

After the radial-swirl correction, the largest k6/Gauss8 point is at
(r,z)=(0.005689,0.003433), inside the axial cutoff collar. Its residual is
mostly radial (+272322), with angular -29815 and axial -39970. This explains
why the radial-collar swirl modes cannot change the total maximum there.
`radial_repair_refit.py` reuses the six compact streamfunction/swirl/pressure
modes as an *increment* with a steeper (tau_ref/tau)^2 weight. The fit uses
full Cartesian momentum at k6/Gauss8 and retains nonlinear self-transport;
the fitted strength1 overfits the independent k6/Gauss6 grid (max94315 ->
107841; sampled L2 191.48 ->218.15). It is rejected.

Strength.1 is the conservative screen: at k6/Gauss8, total max276850 ->266160
and sampled L2 194.37 ->188.34; independent k6/Gauss6 max94315 ->91544,
L2 191.48 ->190.13. At k5.5/Gauss6 max91715 ->88740 and L2 161.13 ->
158.02. At k4 and k=.4 both sampled measures also decrease slightly.
`IncrementalRadialRepairField` exposes the combined compact field with this
increment; analytic divergence freedom follows from its streamfunction and
swirl construction. `radial_repair_direct.py` independently checks six late
peak-neighborhood points using fourth-order Cartesian derivatives and a
one-sided time stencil at the registered endpoint: the worst of these norms
falls276850 ->266160, with finite-difference divergence around1e-6.

The candidate remains massively above1e-3, the quadrature is not converged,
and no critical-time extension or non-axisymmetric stress is supplied. More
small slow-mode fits will not substitute for the paper's realizable oscillatory
stress and exact moment correction cycle in Sections7-9. Preserve this
candidate as a numerical benchmark, not an accepted Navier-Stokes solution.
Artifacts: `experiments/root_st073/compact_potential/radial_repair_refit.json`
and `radial_repair_direct.json`.

## Exact-curl two-harmonic wave: local covariance passes, full momentum fails

At the current k6 radial peak (r,z)=(0.00568908,0.00343263),
`radial_peak_cone.py` computes the actual corrected velocity gradient and the
physical local stress primitive needed for the tangential residual:
(sigma_rtheta,sigma_rz)=(45.5950,138.9583). This primitive is integrated only
from the axis to this radius; nonzero global weighted moments prevent its
unmodified compact extension. The physical analogue of paper Eq(7.1) gives
lambda^2=3919.85, target dot N=-145.02, and cone ratio0.0673<1. Frozen
Cartesian Kelvin rays give two positive time-averaged covariance weights.
These facts are LOCAL diagnostics, not the paper's normalized leading-cone
hypotheses or supported amplitude solution.

`curl_wave_prototype.py` periodicizes the two selected modes to distinct
integer angular harmonics m=1,2, projects their polarizations transverse to
the corrected normals, solves positive covariance weights23.1621 and75.8107,
and defines an actual compact C4 vector potential. Its analytic curl is
exactly solenoidal. At the selected point, 16-angle sampled covariance of the
exact velocity matches the target within2.3e-15 relative error; different m
values remove cross terms under angular averaging. This completes a concrete
local stress realization, but not the wave amplitude/pressure PDE of paper
Section7 or its mean correction cycle in Sections8-9.

The decisive full momentum test rejects this frozen wave. At target amplitude,
the local 16-angle maximum rises266160 ->42781680 and the angular-mean radial
residual rises262085 ->11826714. `curl_wave_amplitude_screen.py` retains exact
linear and quadratic residual terms: among 0 and101 positive amplitudes from
1e-5 to1, zero has the lowest sampled momentum RMS. At amplitude0.001 the
covariance is only1e-6 of target, yet maximum residual already rises to281524.
The finite-difference divergence of the target-amplitude wave is ~2e-4 while
analytic divergence is zero; this reflects differentiation of a steep compact
wave, not a certified numerical divergence bound.

The axial support halfwidth is0.00075, so its inverse scale is1333, versus
carrier magnitudes239 and419. The dimensional cutoff diffusion scale
nu/h_z^2 is17778, versus the local physical analogue growth rate62.6. These
are only scale diagnostics; the paper's normalized estimates cannot be
replaced by this comparison. They explain why freezing a locally favorable
covariance into a narrow compact pulse has a large omitted curl/viscous cost.
Next derive evolving periodic phases and solve the coupled amplitude-pressure
equation before any wave is added to the benchmark; preserve exact radial
moments and re-evaluate the complete residual. Artifacts:
`experiments/root_st073/compact_potential/radial_peak_cone.json`,
`curl_wave_prototype.json`, and `curl_wave_amplitude.json`.

## Separate the explicit wave's nonzero and mean errors

`curl_wave_amplitude_screen.py` now keeps the exact affine expansion of full
momentum for base + a*wave. At the late 16-angle sample, angular means of the
LINEAR wave residual are below1.2e-5 in all cylindrical components (as
expected from distinct nonzero harmonics). Its pointwise linear maximum is
1.6386e7. The QUADRATIC wave residual has angular mean
(1.15646e7,3.14050e6,-1.92154e4) in (r,theta,z). Thus the enormous radial
mean in the target-amplitude wave is quadratic self-interaction, not merely
a poor mean of the linear phase equation. The nonzero harmonic linear error is
also enormous, so both the paper's amplitude-pressure equation and its mean
correction cycle are independently needed.

`curl_wave_mean_pressure.py` adds a compact axisymmetric pressure bubble with
radial derivative -a^2*1.15646e7 at the selected point. This cancels the
wave-induced mean radial residual LOCALLY: at a=1, its angular mean falls
1.18267e7 ->2.62085e5, back near the original background value. The full
16-angle maximum is still3.14580e7 and RMS1.94002e7 (versus baseline
2.66160e5); at a=.1 the maximum is1.94510e6. The azimuthal mean and nonzero
harmonics remain. This pressure is a local proxy, not the exact compact radial
primitive, and no wave amplitude or harmonic pressure equation was solved.
Do not add either wave to the candidate field. Artifacts:
`experiments/root_st073/compact_potential/curl_wave_amplitude.json` and
`curl_wave_mean_pressure.json`.

## Harmonic pressure and local amplitude time slope: local gain, spatial failure

`curl_wave_harmonic_pressure.py` Fourier-decomposes the exact-curl wave's FULL
linear momentum residual at the selected k6 point. For integer angular modes
m=1,2, only14.9% and11.1% of the respective harmonic norms lie along the
local phase normal. Compact pressure harmonics remove those projected parts.
At target wave amplitude, full local maximum changes42781680 ->42719654;
with the earlier local mean-pressure proxy it changes31457954 ->31342811.
Pressure alone cannot cancel the mostly transverse linear error.

`curl_wave_taylor_amplitude.py` adds an exact-curl time-slope vector potential
whose local Fourier curl amplitudes replay the transverse targets to relative
errors below3.2e-16. At the registered endpoint, its added velocity vanishes
but its physical-time derivative cancels the targeted local linear harmonic.
With wave amplitude.1 and both pressure corrections, the 16-angle maximum
falls1963009 ->571391 and RMS1066886 ->299568. Yet the unperturbed field's
local maximum and RMS are both266160. A dense amplitude scan over[-.3,.3]
finds its best maximum at ZERO amplitude. Best RMS occurs at amplitude-.0205,
266160 ->266094, but the maximum worsens to278588 and covariance is only
0.00042025 of the stress target. The polynomial residual reconstruction agrees
with direct full-field evaluation at that amplitude within3.9e-6.

The local Taylor correction is especially fragile in space.
`curl_wave_taylor_neighbor.py` evaluates complete momentum at eight angles on
five meridional positions. At amplitude.1 the center maximum is571391; moving
radially by +/-0.0005 gives about3.58e6/3.61e6, and moving axially by
+/-0.0002 gives about1.29e7/1.35e7. Each exceeds its no-wave baseline by
large factors. These are sample maxima, not global norms, but they reject the
local Taylor wave decisively. Exact curl and a single-point amplitude derivative
are insufficient. A supported solution of the coupled phase/amplitude/pressure
PDE across the whole patch, followed by the paper's mean moment correction,
is required before revisiting the full-field acceptance gates. Artifacts:
`experiments/root_st073/compact_potential/curl_wave_harmonic_pressure.json`,
`curl_wave_taylor_amplitude.json`, and `curl_wave_taylor_neighbor.json`.
