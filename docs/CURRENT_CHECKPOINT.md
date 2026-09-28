# Current checkpoint - constrained velocity-field integration

## Latest continuation checkpoint - 2026-09-27

The stricter `localized_parent_drift_fit.json` is now also independently
replayed: maximum/L2 8.21304e10 / 1.65408e6 at reference and 8.11335e10 /
1.65466e6 at the short endpoint. Its directly evaluated profile drift is
1999.619777783, restoring the original parent's level. It trades some
momentum reduction for this stricter profile cap; retain both candidates
below as documented alternatives. `global_two_patch_candidate.load(source)`
can load this stricter source explicitly. A velocity-only solenoidal scale
map is also available; see [scale transport](ST073_SOLENOIDAL_SCALE_TRANSPORT.md).
That map is kinematic, with no pressure/time dynamics or NS step accepted.

Latest preferred momentum/shape continuation: the drift-capped two-patch
candidate in `localized_two_patch_constrained.json`. Independent actual
18,720-point replay gives maximum/L2 7.49634e10 / 1.63977e6 at reference
and 7.40549e10 / 1.64031e6 at k0+1e-6, improving both metrics over the
one-patch baseline at both times. `global_two_patch_candidate.load()` is
the assembled global field interface. Direct profile drift is about
2122.104, matching the imposed cap but still above the unpatched parent's
1999.620. A stricter parent-drift fit is the next optimization, not an
accepted result. Full-support residual, critical-time behavior and scale
recursion remain unproved; zero recursion steps are accepted.

The preceding continuation baseline was the constrained localized patch
in `experiments/root_st073/localized_constrained_tangent.json`, with its
fixed enriched 264-control parent. Independent actual Cartesian momentum
replay on 18,720 points confirms simultaneous maximum/L2 improvement at
reference and k0+1e-6: corrected maxima 8.54481e10 / 8.54896e10 and volume
L2 1.81103e6 / 1.81146e6 (about 27.7% / 6.1% below the respective parent).
The physical time interval is only 1.6919e-10; no NS step or recursion is
accepted, and the residual target 1e-3 remains unmet by many orders.

`global_localized_candidate.load()` supplies the same inner candidate with
the global compact mean and disjoint exterior collar. Assembly probes
preserve inner/exterior values and vanish beyond the support union at both
tested times. This is a global diagnostic field, not a full-support residual
certificate or a critical-time extension. The older paragraph below about
missing axial localization describes the pre-extension base, not this field.

Prior acceleration fitting had a weighted-objective gradient bug, now fixed.
The corrected candidate improves independent L2 but increases the holdout
peak, so it is not adopted as a joint improvement. Prior unconstrained
acceleration also failed the aspect gate. Prioritize joint optimization of
the two compact spatial patches, retaining moments, positive cone margins,
endpoint shape constraints and a sampled peak cap. See
[current experiment evidence and next solve](ST073_LOCALIZED_CONSTRAINED_NEXT.md).

## Earlier ST073 research update - 2026-09-27
A globally localized diagnostic now combines the balanced compact wave,
inner acceleration correction and exterior collar tangent through one
callable field. Sampled assembly checks at reference and endpoint preserve
the constituent inner/outer values exactly and return finite axis values.
This does not establish whole-domain momentum or recursion. See
[unified global candidate](ST073_GLOBAL_ACCELERATION_CANDIDATE.md).

The inner acceleration diagnostic reduces sampled endpoint L2 by 4.78%; its
actual shape replay fails relative axial elongation, so constrained acceleration fitting is required before adoption. Independently, exterior collar
L2 falls 7.97% on a separate replay grid. These percentages must not be added.

Global-scope audit confirms a separate required gap: the current base is
registered only on abs(eta)<=0.5, and its radial heat exterior has no axial
localization. The candidate is not yet a globally defined finite-energy
field. A divergence-preserving axial exterior/pressure construction is
being investigated alongside local residual correction. See
[global exterior gap](ST073_GLOBAL_EXTERIOR_GAP.md).

Degree-3 mean enrichment with direct endpoint geometry now reduces
independent actual momentum L2 to 1.92894e6 and maximum to 1.18244e11
(about 27% and 30% below the prior 236-control candidate). All three
actual sampled forward geometric signs pass. Enlarged-mean independent
replay gives moment maximum 2.1084e-7 and 27/27 cones passing, minimum
margin 9.94824e-5. Higher-quadrature refitting now lowers same-grid L2 from 1.92798e6 to 1.89886e6, but increases the peak from 1.13931e11 to 1.15611e11. The peak-capped solve improves independent L2 to 1.89999e6 but worsens its maximum to 1.20729e11. It is not a joint improvement; peak constraints must cover both grids before a separate disjoint check. See [refined momentum](ST073_REFINED_MOMENTUM.md).
No NS step or scale recursion is accepted. See
[enriched mean endpoint](ST073_ENRICHED_MEAN_ENDPOINT.md).

Degree-3 mode-2 enrichment now combines with shape constraints: independent
actual-field momentum L2 is 2.63329e6 and max 1.68320e11, while all three
sampled forward geometric directions remain favorable over delta-k=1e-6.
Moments/cones pass the assembled solve. Mean harmonic now contributes
80.33% of squared residual and is the next enrichment target. This remains
an affine candidate, not a momentum-accepted trajectory or recursion. See
[enriched shape tangent](ST073_ENRICHED_SHAPE_TANGENT.md).

Larger shape-direction margin now gives all three requested sampled forward
changes over delta-k=1e-6: radial RMS decreases, aspect increases, and
weighted angular speed increases. Peak swirl still decreases. This is an
affine candidate increment over physical delta-t about 1.692e-10, not an
accepted NS step or scale recursion. Independent reference-time replay now
gives momentum L2 2.84751e6, max 2.10261e11, moment max 2.1176e-7 and all
27/27 cones passing. Endpoint momentum/matching remain unaccepted. See the updated
[shape-constrained result](ST073_SHAPE_CONSTRAINED_TANGENT.md).

Shape-constrained tangent now preserves independent moments and 27/27
cones while giving all three requested central derivative signs. However,
the forward delta-k=1e-6 interval still expands radius and lowers aspect:
curvature overwhelms the small direction margin. A stronger-margin
candidate is generated and under actual-field replay. No geometric step
or scale recursion is accepted. See
[shape-constrained tangent](ST073_SHAPE_CONSTRAINED_TANGENT.md).

Separately, spatial degree-3 mode-2 enrichment reduces independent-grid
L2 from 2.8472e6 to 2.6329e6; degree 4 is slightly worse on that grid.
This remains to be combined with shape constraints. See
[higher harmonic tangent](ST073_HIGHER_HARMONIC_TANGENT.md).

Target-direction check now rejects advancing the cone-compatible tangent
as successful recursion: fixed-cylinder radial RMS increases and axial/
radial aspect decreases, although weighted angular speed increases.
These are sampled affine-tangent trends, not an NS trajectory. Shape-rate
constraints must join the next tangent solve; residual and matching alone
are insufficient. See [wave shape direction](ST073_WAVE_SHAPE_DIRECTION.md).

Joint moment/cone correction now passes independent actual-field cone
replay at all 27/27 locations (81/81 inequalities, minimum margin 9.9759e-5).
Order-96 integral moment maximum is 2.1219e-7. Independent full momentum
L2 is 2.8472e6 and maximum 2.1019e11: restoring cones costs 6.98% in L2
against the moment-only candidate, retaining about 61.4% reduction versus
the old dense tangent baseline. No recursive step is accepted. The next
construction varies wave shape with both moments and cones in the inner
solve. See [joint moment/cone result](ST073_MOMENT_CONE_WAVE.md).

Moment-constrained wave correction now restores actual order-96 integrated
moments to max 2.54e-8, versus 583600.74 for the earlier candidate, while
independent full momentum L2 changes only from 2.6453e6 to 2.6614e6.
Cone replay improves to 12/27 locations passing but remains a failure.
Higher-grid wave growth is +997.0958. No recursive step is accepted.
See [moment-constrained wave result](ST073_MOMENT_CONSTRAINED_WAVE.md).

Momentum-aware wave-shape optimization now reduces independent actual
patch volume L2 from 7.3766e6 to 2.6453e6 and maximum from 5.6409e11 to
2.0098e11, about 64% for both. Five-node flux constraints remain matched.
The growth sign discrepancy is traced to coarse quadrature of nearly
cancelling production/dissipation terms. Higher z20/r11 quadrature gives
lambda +997.096 for the new wave; positive growth has this finite-grid
support, but the 1000 floor and continuum convergence are not verified.
Renewed mean compatibility remains required. See
[momentum-aware wave construction](ST073_MOMENTUM_AWARE_WAVE.md).

Denser harmonic fitting now improves the new order-13 independent physical
replay: frozen-wave patch L2 9.9155e6 becomes 7.3766e6, with maximum
5.6409e11. The earlier severe fitting instability is reduced, but the
absolute residual remains unacceptable. The next construction optimizes
wave shape against full momentum while retaining local flux/growth targets.
See the denser correction in [joint wave/mean status](ST073_JOINT_WAVE_MEAN_STATUS.md).

The old growing-wave joint 45-control mean fit has completed spatial replay:
annulus volume L2 worsens from 193257.9992 to 193600.2673, despite a small
patch improvement. Do not extend this candidate as a successful recursive
step. See [joint mean/wave spatial result](ST073_JOINT_WAVE_MEAN_STATUS.md).
Full harmonic tangent replay on the locked co-designed candidate now fails
the spatial holdout: volume L2 is 1.1175e8 versus frozen-wave 1.0036e7 and
mean-only 1.2146e5. The independent harmonic budget identifies the second
harmonic as 98.52% of the fitted squared residual. Denser spatial fitting
is next; no time step or scale recursion is accepted.

Joint polarization/growth fitting now matches the required radial flux at
five nodes with relative error 4.69e-11. Higher-quadrature growth lambda is
997.36, positive but below the imposed 1000 floor. Wave RMS is 4.347 times
mean RMS, so this remains a large nonlinear candidate requiring complete
residual correction. See [stress and growth co-design](ST073_STRESS_GROWTH_CODESIGN.md).
No accepted physical trajectory or scale recursion follows from this fit.

The actual physical mean-plus-wave replay is now available in
`mean_wave_replay.py/.json`. The old scalar-amplitude wave candidate fails
to improve the complete residual: independent patch volume L2 increases
from 120498.06 to 120955.20 after fitting oscillatory time/pressure controls.
This rejects a conclusion based on its small angular-mean improvement.
See [physical mean-wave replay](ST073_PHYSICAL_MEAN_WAVE_REPLAY.md).

The constrained short step retains a positive mode-1 instantaneous energy
rate on the fixed patch: split moderate/higher quadrature gives 13161 and
14501, versus initial 19027 and 20438. The 9.24% endpoint sensitivity is
not convergence evidence. The two unconstrained delta-k=0.01 paths instead
lose growth in all eight tested modes. See
[the evolved growth screen](ST073_EVOLVED_WAVE_GROWTH.md). Joint stress and
growth co-design is now the next construction; energy positivity alone is
not covariance matching or actual wave integration.

The first dynamically constrained state step is complete at delta k=0.0001.
Recomputed endpoint compatibility passes 27/27 sampled cones, with integrated
moment maximum 8.88285e-8; endpoint momentum maximum/L2 are 4.43477e9 and
188215.12. This is initial/endpoint local-tangent compatibility, not a
continuous-time or recursive certificate. Fixed-cylinder spin does not
increase. See `meridional_constrained_evolution.json` and the
[state evolution record](ST073_STATE_EVOLUTION.md). Next investigate actual
wave growth/coupling rather than extending mean relaxation alone.

State-dependent unconstrained evolution now improves on frozen slopes over
delta k=0.01: the independent endpoint maximum is 2.27427e9 versus
6.00640e10, and volume L2 is 84464.34 versus 1473368.61. This uses a
165-point fit and two explicit Euler steps, with the full quadratic
nonlinearity. It is not a matched or converged trajectory. A coarse
44-point fit gave misleadingly poor evolution because its axial sampling
aliased basis directions; the dense unregularized matrix is full rank.
See [state evolution and its limitations](ST073_STATE_EVOLUTION.md).

Important target check: on a fixed physical cylinder, that refreshed
delta-k=0.01 path broadens the enstrophy radial RMS from 0.0005595 to
0.0016775, lowers the axial/radial ratio from 0.9337 to 0.3904, and lowers
weighted angular speed from 2.125e6 to 0.968e6. The residual improvement is
not evidence of the requested vortex amplification. The diagnostic is
`vortex_state_observables.json`; it describes the sampled cylinder, not an
identified global core. Do not extend this unconstrained mean-relaxation
path as if it were successful scale recursion.

The constrained 44-direction replay is now complete. On the independent
176-point spatial holdout, momentum maximum decreases from 7.01009e9 to
4.50644e9 and physical-volume L2 from 237497.46 to 188297.75. All 27 sampled
stress cones pass at radial quadrature order 64; the four integrated moments
have maximum absolute value 1.47389e-5 at order 96. The selected solver point
is feasible and improves the objective, but SLSQP returned status 8, not
convergence. See `broad_meridional_constrained_replay.json`. This supersedes
the pending compatibility statement below, not the full PDE/recursion gates.

An opt-in grouped evaluator preserves the original fields and tested Cartesian
jets while reducing the cold 2000-point field benchmark from 5.2885 s to
0.1337 s (39.54 times). The full constrained replay took 121.88 s. This is
computational acceleration, not a reduction of the governing equations.
The next experiment evolves the velocity coefficients and recomputes their
slopes at each state; frozen initial slopes are already rejected below.

For that experiment, `integrated_state_moments.py` now assembles current-state
moment equations directly from fixed-radius conservation integrals. At a
nonzero velocity state, its four predicted moments agree with independent
actual-field integrated replay to 3.59e-8 (order 96). This removes the need
to use the less accurate order-24 pointwise-residual integral rows as the
evolution constraints. It does not establish continuum derivative accuracy.

The first full-vector meridional/pressure correction now reduces both metrics
on the same independent spatial holdout: maximum 7.01009e9 -> 4.38510e9,
physical-volume L2 237497.46 -> 186820.70. Its 44 directions preserve the
instantaneous broad-shear velocity and fit its time derivative plus pressure.
This is an unconstrained instantaneous experiment: the earlier moment/cone
passes do not transfer, and finite-difference divergence remains 0.154809.
The next solve must retain the new full-vector objective while reinstating
moment/cone compatibility. No scale recursion or PDE acceptance is established.
See [the meridional correction record](ST073_MERIDIONAL_MOMENTUM.md).

Off-time replay of that unconstrained affine-k candidate preserves its local
benefit through delta k=0.001, but direct extension to delta k=0.1 and 1
fails badly (momentum maxima 1.52546e12 and 3.82489e14). These are callable
evaluations, not time integration. The combined broad-amplitude slope crosses
zero at delta k=0.034086. Do not use frozen initial slopes as a recursive
trajectory; the next evolution needs state-dependent updates and full nonlinear
terms. Full records are `broad_meridional_time_audit.json` and
`broad_meridional_scale_audit.json`. A separate constrained 44-direction solve
is in progress; no compatibility result is claimed until its replay completes.

The next state-dependent evolution now has a factorized nonlinear residual
engine, `meridional_state_cache.py`. It supports 25 velocity values, their
25 time slopes and 19 pressure coefficients, retaining all quadratic
advection. Two nonzero-state five-point comparisons agree relatively to
2.25e-11, but absolute discrepancies reach 0.02481; this is not final-gate
accuracy or an integrated trajectory. See
[the state-cache record](ST073_MERIDIONAL_STATE_CACHE.md).

A new broad annular r^-2 swirl direction admits positive instantaneous
wave energy growth: at amplitude 44.8774 the mode-1 rate is 1.0217e4
on the refined split quadrature, with all 22 sampled centrifugal growth
conditions positive. This changes the earlier all-decaying energy result
by changing the mean, not by changing viscosity or numerical time steps.
Its added local swirl RMS is about 50 times the old mean velocity RMS.
It is not an accepted low-residual replacement, and modes 2 through 8
still decay. The mean compatibility result is recorded in
`broad_shear_dynamic_control.json`: adding da/dk=-227.669 fixes the wave
region's stress sign. Independent replay gives moment maximum 1.45e-5,
6/6 outer cone passes and 5/5 wave-region cone passes. The instantaneous
growth matrix is preserved, but full sampled momentum remains 7.46e9.
The earlier 18-control seed and its 0/5 local cone result are historical;
see [the original co-design record](ST073_SHEAR_CODESIGN.md).

Radial pressure primitives are now implemented. At the old centrifugal
peak an inner-datum repair lowers the norm from 7.48e9 to 1.55e9, but
creates exterior axial imbalance. A compact primitive restores exterior
pressure while retaining a 4.71e9 collar defect. The numerical radial
integral budget shows that pure swirl plus compact pressure cannot remove
this added radial imbalance. The next construction needs a divergence-free
meridional time correction and/or actual wave stresses, coupled to the
pressure source. See [dynamic and pressure matching](ST073_BROAD_DYNAMIC_MATCHING.md).

Full spatial potential evolution now includes the angular mean and modes
1 through 8, cutoff derivatives, nonlinear transport and viscosity.
An implicit weak-diffusion BDF solve fixes the explicit-step instability.
Three Cartesian time replays have momentum maxima 7.80e9, 3.79e9 and
1.91e9, still far worse than the background-only field. The oscillatory
energy proxy falls to 1.44% of its initial value; initial center covariance
is not retained. This is numerical integration progress, not successful
scale recursion or recursive amplification. See
[the full spatial evolution and decay audit](ST073_SPATIAL_FOURIER_EVOLUTION.md).
The initial-time energy screen now finds no positive growth direction
in any retained mode at support multipliers 1, 2, 4 or 8. At multiplier
8 the best combined amplitude rate is still -2.63e5, while the largest
strain-only rate is 1.42e4. See [the energy budget](ST073_WAVE_ENERGY_BUDGET.md).
The broad-shear seed above passes the initial energy-growth screen; the
remaining blocker is making that growth compatible with the full mean
and wave dynamics and the required two-component stress.
The paper requires growth followed by decay; decay at a pulse tail alone
is not a failure criterion.

A moving-normal principal amplitude/pressure inverse is implemented and
passes manufactured checks, but its current center-path spatial wave is
rejected. Full momentum on spatial holdouts grows from 5.91e9 for the
moving wave to 8.86e9-1.91e10 after its forced correction. The background
alone is about 4.50e5 on the same local grid. This trial does not replace
the feedback mean trajectory. See [the rejected wave experiment](ST073_MOVING_NORMAL_WAVE.md).

The new finite-dimensional spatial solve addresses these omitted operators
but does not provide the paper's supported inverse or recursive estimates.
Do not repeat center-only inverse sweeps or interpret a principal ODE check
as a full-field residual reduction.

State-dependent pressure/swirl-slope integration now completes the short
interval k in [11,11.001]. Three direct time replays have integrated moment
maxima 9.87e-4, 4.99e-4 and 9.57e-4, with 6/6 outer cone samples passing
each time. At k=11.0005, the matched fixed-slope comparator is 1.658,
so feedback reduces that moment drift by about 3,322 times. See
[the completed trajectory diagnostic](ST073_STATE_DEPENDENT_EVOLUTION.md).

This is sampled short-time compatibility, not scale recursion. The worst
moment replay is close to 1e-3 without a verified error margin, spatial
neighborhood holdouts have not been replayed along this trajectory, and
full sampled momentum remains 1.42e6. Full momentum max/volume-L2, finite
energy, forcing/domain requirements and recursive contraction remain
unestablished. No candidate is accepted.

The preceding spatial seed passes eight neighborhood holdouts and nine
inner nodes at k=11. Shared radial integration makes the 22-node benchmark
2.53 times faster. The code and all completed trajectory/comparator
artifacts are in experiments/root_st073.

The snapshot below is historical, not current ST073 acceptance.

Snapshot date: **2026-09-20**.

This integration-branch snapshot predates the 2026-09-22 `main` research
consolidation. On `main`, ST006 is the historical runnable API baseline,
while ST061-D/P are newer residual-oriented controls and ST063-G2R is the
latest documented geometry experiment. The paired ST061/ST063 samples
report full-vector maxima about `0.0225`â€“`0.0272` and volume L2 about
`0.033`â€“`0.034`; their seeds and candidate artifacts differ from the
ST006 protocol below. None meets both `1e-3` momentum gates. See
`main:docs/FINAL_RESEARCH_SNAPSHOT_2026-09-22.md` and
`main:docs/RESEARCH_STATUS.md` for the original study records and
candidate availability. Do not rank them against ST073 annulus screens
without a shared complete-field validation protocol.

Machine-readable authority: [`project_status.json`](../project_status.json).
Active integration branch at snapshot: `codex/cr001-constraints@9491103442e53190b7d90d6f0fc5f3165aacc697`.

This file is intentionally a **current** checkpoint. Historical optimization experiments, rejected basis studies, exact-reconstruction notes, and earlier route-local diagnostics remain available in Git history and artifacts, but they are not merge blockers for the current delivery chain unless they affect the frozen candidate being shipped.

## Final delivery target

The nine lanes converge on a directly usable time-dependent Cartesian velocity field

`velocity(x, y, z, t) -> [u, v, w]`

that can be called from Python, saved/reloaded, exported on reproducible grids for MATLAB/Python, and inspected with streamline/vorticity diagnostics. Its publicly observable geometry and time evolution should be made progressively closer to the public OpenAI velocity-field visualization.

This target is **not** a paper-exact reconstruction target and **not** a complete blow-up proof target.

The three project states are independent:

- `velocity_export_ready`
- `visualization_ready`
- `pde_validated`

CI success, rendering success, optimizer convergence, or visual resemblance does not promote the other states automatically.

## Repository-wide PDE baseline

The retained same-protocol numerical baseline is ST006 on `main`:

- candidate SHA-256: `6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3`
- held-out seed: `9172801`
- held-out Cartesian points: `4096`
- validation times: six
- finest spatial step: `0.005`
- momentum sampled max: `0.1082289305112118`
- momentum volume-L2: `0.10758432876230622`
- registered momentum target: `1e-3`
- `pde_validated=false`

Only directly comparable complete-residual protocols may claim improvement over this baseline. Scoped Kokuno terms, morphology diagnostics, sampled derivative audits, or render metrics are not directly comparable unless they reproduce the same residual contract.

## Canonical constrained delivery

The current canonical constrained delivery remains `eq45_supported_velocity_candidate_v1`.

- public evaluator: `openai_ns_reconstruction.eq45_supported_delivery:velocity`
- save/load: `Eq45SupportedDeliveryField.save_candidate(...)` / `Eq45SupportedDeliveryField.load_candidate(...)`
- grid evaluator: `default_field().grid`
- delivery capsule: `artifacts/constrained/eq45_supported_delivery_capsule.json`
- supported child SHA-256: `2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d`
- parent SHA-256: `48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7`
- registered box: `[-2,2]^3`
- registered time interval: `[0.25,0.75]`
- physical support connection: `r < 2`, `|z| < 2`

Current canonical truth state:

```text
velocity_export_ready = true
visualization_ready   = false
pde_validated         = false
```

`velocity_export_ready=true` means the canonical Eq45 field is callable/saveable/exportable. It does **not** mean visual correspondence, exact OpenAI-field identity, PDE acceptance, paper exactness, or blow-up.

## ST052-M visualization-candidate delivery state

The frozen ST052-M child is no longer awaiting materialization. The following delivery pieces are already live on constrained ancestry:

1. whole-child exact-source-backed runtime and save/load identity;
2. checksum-bound `33^3 x 5 times` sampled-grid export;
3. NPZ and MATLAB-v5 MAT outputs;
4. GNU Octave MAT-load/render smoke;
5. live delivery-identity reconciliation.

Key merged integration PRs:

- #665 â€” ST052 grid export materialization;
- #680 â€” GNU Octave consumer/load/render smoke;
- #706 â€” live ST052-M grid delivery identity reconciliation;
- #777 â€” machine-state synchronization removing stale rematerialization routing.

The #706 repaired exact head `62ed0677675cbba3ff409a6ba124bc157bc5b504` completed repository `tests` workflow run `35475995041` with conclusion `success` before merge.

ST052-M still has a narrower delivery gap than the canonical Eq45 field: the repository has not closed or explicitly accepted the **standalone-package parent-runtime/dependency identity seam** as sufficient for a standalone ST052 export claim. Native MATLAB scientific execution has also not been performed for the ST052 grid; the current evidence is GNU Octave consumption.

Therefore ST052-specific state remains:

```text
source_runtime_backed_callable_save_load_ready = true
exact_source_runtime_identity_closed           = true
grid_export_materialized                       = true
octave_mat_load_render_smoke_passed            = true
standalone_package_parent_runtime_ready         = false
velocity_export_ready                           = false
visualization_ready                             = false
visual_correspondence_verified                  = false
pde_validated                                   = false
```

Do **not** open another ST052 whole-child materialization, NPZ/MAT exporter, or Octave replay lane. Those seams are already closed.

## Current CI truth

After #777 merged, the live constrained head is

`9491103442e53190b7d90d6f0fc5f3165aacc697`.

Its repository `tests` push workflow is run **35507408812**. At this snapshot it is still **queued**, so the live merge head must not be described as newly CI-green yet.

Runner backlog is not scientific evidence. A queued workflow is neither PASS nor FAIL.

## Active sibling delivery assets

### Agent 9 â€” source-observable and delivery diagnostics

`main` already contains ST054 delivery plumbing for:

- continuous Python callable;
- Python streamline/vorticity rendering;
- deterministic NPZ/MAT handoff;
- VTK/ParaView handoff.

Current open Agent-9 work includes:

- #825: forward-port the official-public qualitative observable contract to current `main`;
- #834: renderer-independent cylindrical morphology fingerprint for ST054.

At this snapshot their current exact-head workflows are queued. These are useful sibling assets and reusable diagnostics, but they do not silently replace the constrained Eq45/ST052 candidate identity.

### Agent 7 â€” morphology/capacity lane

Open Agent-7 work studies outer-reservoir morphology, temporal curvature, representation equivalence, and a compact toroidal swirl preflight. These are candidate-capacity or morphology experiments. They should not be promoted into canonical delivery merely because a local capacity gate passes.

New basis growth is not the shortest constrained-delivery path while the frozen ST052 runtime/diagnostic chain is still incomplete.

### Kokuno Agents 1â€“5

The Kokuno lanes are advancing an executable **strict-inner** PA.10 contraction-center stack: velocity, time derivative, spatial derivatives, self-advection, oscillatory composition, cylindrical mean projections, and independent derivative/advection audits.

They remain inner-only and do not yet provide all of:

- outer/global leading join;
- corrected/global candidate velocity;
- matched pressure;
- preregistered restricted forcing;
- complete `velocity/pressure/forcing` candidate API;
- same-protocol complete NS residual;
- final divergence-L2 / momentum max/L2 acceptance.

These lanes are valuable upstream scientific work but are not blockers for exporting a clearly labeled constrained visualization candidate.

## Next integration task â€” shortest delivery chain

The next constrained integration task is **not** another basis experiment.

For the same frozen ST052-M candidate:

1. close **or explicitly accept** the standalone-package parent-runtime/dependency seam;
2. preserve one unified `velocity(x,y,z,t)` identity through save/load;
3. replay the already-live `33^3 x 5` NPZ/MAT export;
4. run fixed-seed/fixed-camera streamline and vorticity diagnostics;
5. run source-observable diagnostics without inventing public numerical targets;
6. produce one Python/MATLAB-facing integration report.

If that candidate becomes a stable visualization candidate, it may be exported and labeled as such even while `pde_validated=false`.

Only if this exact frozen candidate is selected for PDE work should the integration route rebuild compatible pressure/restricted forcing and run a fresh independent 4096-point momentum/divergence validation before any `pde_validated` promotion.

## Duplicate-work firewall

Do not duplicate:

- ST052 whole-child materialization;
- ST052 NPZ/MAT grid export;
- ST052 GNU Octave MAT consumer smoke;
- ST054 Python render / NPZ-MAT / VTK plumbing;
- Agent-9 public-observable contract work;
- Kokuno strict-inner derivative/advection subterms already owned by Agents 1â€“5.

New work should close a currently open seam in the shortest `[u,v,w]` delivery chain or provide a genuinely independent validator for a selected candidate.

## Truth boundary

The following implications are forbidden:

```text
CI green              != PDE-valid
optimizer converged   != PDE-valid
export works          != visual correspondence
visual resemblance    != PDE-valid
visual resemblance    != exact OpenAI field
sampled local audit   != whole-domain validation
PDE-valid             != blow-up proof
```

The project should prefer a reproducible, clearly labeled callable visualization candidate over blocking all export on unresolved historical proof goals, while keeping PDE and exact-source claims fail-closed.
