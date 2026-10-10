# Current original Rp correlated physical source delivery (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. [Source c8808a2a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c8808a2aa3dfd4b6bb1ea95a7936db059a40e089) introduces a source-bound original-scale physical log-input family, proves its exact inverse using the original implicit equation and uniqueness, routes the resulting coordinates into actual current source callbacks, and maps them through unchanged Cartesian spatial/time operators. Six actual original-scale physical inputs now deliver source packets and 840 spatial plus 24 fixed-x time rows in exact signed scale groups. These remain scaled enclosures. Numerical physical u/v/w/p, unrestricted inputs, global/axis coverage, stress, actual coefficient recursion and oscillatory correction remain open.

## Callable current interface

```python
owner = CurrentOriginalRpCorrelatedPhysicalSource(accepted_native_box_source)
point = owner.physical_input('heat_exterior', '9/2', '1/2', '0', '7/10')
inverse = owner.invert(point)
located = owner.locate(inverse)
delivery = owner.source(inverse)
rows = owner.physical_rows(delivery)
rows['physical_velocity_pressure']['ux'].groups()
rows['Cartesian_spatial_rows']['ux']['x1_y0_z0']
rows['fixed_x_time_rows']['ux']
```

Implementation and receipt quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_correlated_physical_source.py`, `.json.gz`, `_check.py`, `_check.json`. Reuse accepted current runtime objects. Original source selection, five-moment repair, future/heat construction and old 33-chart registry are not replayed. The current exact-rational physical source-map guard is unchanged; the new correlated caller has its own live source guard.

## Supported physical input family

Choose an admitted original source chart, exact rational native coordinate v, axial similarity coordinate Z with -1<Z<1, finite exact log-time offset s and angle theta. Let L_R be the **actual original source logR** from that chart. The exact physical log inputs are

```
log_tau = -L_R + s
q_forward = (log_tau - log(1-Z^2))/2
log_r = q_forward + (log(2)+L_R)/2
log_abs_z = log(abs(Z)) + (1-delta)*q_forward   # Z != 0
physical_z_sign = sign(Z)
```

The exact physical time is t=1-exp(log_tau), with 0<tau<1 on this admitted near-critical family. The viscosity is the original nu=1. For Z=0, physical z is exactly zero, log_abs_z is omitted and q=log_tau/2. Original delta, original N, actual Rp, analytic P0 and same selected repair are retained.

These are genuine physical log inputs, with their original enormous time/axial scales retained as exact functions. They are a **restricted correlated family** parameterized through current source coordinates; they are not arbitrary independent x/y/z/t inputs. Physical log_r reduces exactly to `(s+log(2)-log(1-Z^2))/2`, so its finite information survives despite huge original L_R. No convenient replacement radius, independent interval oracle, rounded t=1 or axial z=0 approximation is used.

## Original inverse identity and exact cancellation

For nonzero physical z, the original implicit inverse is retained as its existing exact root node:

```
2q = logaddexp(log_tau, 2log_abs_z + 2delta*q)
F_prime >= 2*(1-delta) > 0
```

The source-correct forward exponents give

```
exp(log_tau) = (1-Z^2)*exp(2q_forward)
exp(2log_abs_z + 2delta*q_forward) = Z^2*exp(2q_forward)
(1-Z^2) + Z^2 = 1
```

Thus q_forward solves the original equation; original strict monotonicity proves it is the same unique root. The original root node remains in the result, together with a separate equality certificate. A live, registered equality alias is used **before** rational graph normalization and directed numeric arithmetic. It is not a replacement root chosen from a coarse interval. Exact Z and the source logR are recovered with this proof, retaining nonzero/negative/zero axial cases.

Source positive quotients are normalized as numerator times a shared exact reciprocal. This identifies the original `v/mu` and `v*(1/mu)` forms without numerical approximation. Bound original parameter functions, particularly integral-defined waiting W, remain atomic; their defining integrals are not expanded into unsupported numerical oracles or replaced by endpoints.

The native affine inverse is separately proved: `logR-base=J*v`, J>0. Request geometry and inverse-token geometry therefore have the same exact logR. This identity remains attached to the source delivery and must accompany future physical scale evaluation. The source native enclosure is refined by intersecting two enclosures of that **proved same function**, one from directed native location and one from the exact source-coordinate identity. This reduces dependency/rounding error without selecting a midpoint or endpoint as a field value.

## Actual source and physical operator delivery

Inputs, inverses and source deliveries each require unchanged live objects owned by this exact caller. Fingerprints bind physical input function nodes, geometry, inverse mapping, equality aliases, signed source coefficient rows and scale metadata. Imported/copy/foreign-owner/mutated objects are rejected. All possible native chart candidates remain in the locator result; source dispatch requires complete native and axial admission on the original request chart. Separate seam/crossing protocols remain open.

The actual current whole native/Z enclosures reach the same pulse, postpulse and closed-heat algorithms. The scoped physical proxy then invokes the **unchanged** `CurrentOriginalRpPhysicalSourceMap.map_source_view`, including its original moving basis, viscosity-one spatial operators and fixed-physical-x time derivatives. Coefficient enclosures retain the source and coordinate errors. Exact point Z/native metadata comes from the proved inverse identity, never a numerical midpoint. The time rows are the original fixed-x partial derivatives evaluated at these physical points, not derivatives along the correlated input family.

## Focused evidence

- Six original-scale input/source/operator routes: active pulse exit, flatten, outer angular correction, waiting, finite heat collar and infinite heat exterior. They return actual source packets, 840 spatial rows and 24 fixed-x time rows.
- Eight exact physical inverse identities include negative axial and zero axial cases. Independent symbolic checks recover source logR, signed Z, finite physical log_r, exact native coordinates and request-versus-box geometry; 62 current graph identities pass.
- Independent algebra verifies the original root factorization and monotonicity margin. Three finite-log original numerical root fixtures use the **current original delta**, solely as consistency evidence; they do not define original scales or source values.
- Source/scaled-row contexts, exact original units, same-graph factors and complete ordinary derivative grids pass. Physical value materialization and physical accuracy remain explicitly false.
- Thirteen inadmissible-input checks cover copied/foreign inputs, inverse and source packets, axial endpoint, old/unknown chart, imported native intervals, graph log oracle, another source chart and altered live geometry/Z.
- Accepted constructor and a fresh off-axis heat-exterior Z=0 public source/physical-row call pass. The sole read-only reviewer remains GPT-5.6 Luna / max and found no material mathematical/provenance defect. Dependency bytes are matched against the Git index before upload.

**Known amplitude-binding gap:** the current interval inlet carries the same-family legacy `CompliantPowerInletC4.constants['U']` box (in Utheta/Pstar units), copied into current selected pulse owners. The current exact `U0` is separately defined by substituting Z=0 into `actual_identified_Rp_native_frame_C3['u']`. Existing checks identify the symbolic function/family and use the legacy box for interval rows, but do not yet export a direct numeric inclusion proof from that exact current function to the U box, or a source-bound directed logU0. This is not evidence that the legacy U values are wrong; it is missing evidence needed for certified physical amplitude/error. The bounded callback/operator gates do not claim that missing physical certification. It is the first next implementation task.

## Completed bounded tasks

- [x] **CURRENT-RP-SOURCE-BOUND-CORRELATED-INPUT-FAMILY** Original log-time/log-radius/log-axial correlations with exact current graph provenance and original scales.
- [x] **CURRENT-RP-CORRELATED-ORIGINAL-INVERSE-IDENTITY** Original root retained; separate source-correct equality proof and uniqueness; nonzero, negative and zero axial identities.
- [x] **CURRENT-RP-CORRELATED-NATIVE-IDENTITY** Exact affine inverse and request/box logR equality with positive J; identity-based enclosure refinement.
- [x] **CURRENT-RP-CORRELATED-ACTUAL-SOURCE-DELIVERY** Six genuine original-scale physical inputs consume live current interval source callbacks.
- [x] **CURRENT-RP-CORRELATED-SCALED-PHYSICAL-ROWS** Current physical operators consume these live source packets and retain signed exact scale groups for velocity/pressure/spatial/time rows.

## Next tasks and acceptance

Prioritize actual value/error delivery next. Mark only achieved bounded tasks complete and preserve the full reconstruction objective. Reuse accepted runtime state instead of repeating source/repair/future solves for a receipt.

- [ ] **CURRENT-RP-TRUSTED-AMPLITUDE-SCALE-BINDINGS** Implement a typed direct evaluation/inclusion proof for raw `U0` / `logU0`, the actual original inlet function at Z=0. Prove that the live U box encloses `evaluate(actual_identified_Rp_native_frame_C3['u'], Z=0)` in Utheta/Pstar units, or rebuild the interval amplitude from that exact source with rigorous error. Trace `current_original_Rp_native_constants.py:75-85`, `current_original_Rp_interval_inlet.py:65-71` and `power_inlet_C4.py:46-49,70-80`. Same-family/source/datum equality and symbolic constant equality alone are insufficient for this numeric inclusion. Bind every other physical amplitude atom with the same source/repair; preserve independent P0 from normalized analytic datum and the Pstar^2 scale. Caps, old fixtures, selected endpoints and archived point rows cannot define these amplitudes.
- [ ] **CURRENT-RP-EXACT-GROUPED-SCALE-EVALUATOR** Evaluate complete source + physical_radial + physical_lambda log-scale sums using the stored original root and native equality proofs. Preserve exact shared L_R, logRp, pulse terms, inlet amplitude, delta and time factors before directed arithmetic. Do not sum each enormous factor numerically first. Support signed coefficients and exact integer powers/positive denominators.
- [ ] **CURRENT-RP-SIGNED-LOG-VALUE-DELIVERY** Return each physical component as signed directed coefficient times exact positive scale groups; add groups only after establishing common exact units. Where feasible, expose directed log-magnitude/mantissa representations without huge exponent materialization. Handle coefficients straddling zero honestly. Preserve component sign and cancellation, particularly ux/uy moving-basis terms and independent pressure baseline.
- [ ] **CURRENT-RP-PHYSICAL-ROW-MATERIALIZER** Materialize finite u/v/w/p and needed derivatives when their full grouped log-scale/coordinate/source bounds support the requested absolute/relative accuracy. If the original scale or uncertainty is outside feasible materialization, return explicit factored/log representations and an unresolved accuracy reason. Do not label scale groups as already produced ordinary Cartesian numbers.
- [ ] **CURRENT-RP-PHYSICAL-POINT-ERROR** Report actual phase, quadrature, source coefficient, Gamma/tail, coordinate and scale errors in physical units. Refine source intervals and true parameter bounds when needed; finite-width q enclosures cannot be naively exponentiated beside astronomical original scales. Use exact correlation and variation bounds, not a chosen coefficient or root.
- [ ] **CURRENT-RP-UNRESTRICTED-CORRELATED-INPUT** Extend beyond the current source-parameterized input family where a true inverse is needed. Accept supported exact physical log graph expressions only through explicit provenance and a correct inverse/coordinate enclosure proof. The current fixed family does not satisfy arbitrary independent x/y/z/t admission.
- [ ] **CURRENT-RP-INVERSE-BOX-SOURCE-SUCCESS** Finish the broader unified physical input caller, including ordinary inverse/refinement when it reaches outer charts and correlated families beyond the current one. The present six original-scale family routes are complete; do not redo their root identities. Preserve every uncertain chart and the distinction between native cover and conditional clipped overlap.
- [ ] **CURRENT-RP-PHYSICAL-BOX-MAP** Extend physical mapping from proved exact points to nondegenerate correlated physical coordinate boxes with same-graph Z/time functions, physical operators, error propagation and sign-fixed/nonzero/zero branches. Current independent similarity `source_cell` boxes are not physical coordinate boxes.
- [ ] **CURRENT-RP-BOUNDARY-AND-SUPPORT-PARTITION** Split chart/support crossings at exact boundaries, call accepted fourteen seam protocols and convert to common units before hulling. Source dispatch must not pick a possible chart or discard a conditional overlap.
- [ ] **CURRENT-RP-ANGLE-EVALUATOR** Bind actual physical atan2/sin/cos and moving-basis boxes, including branch/support changes, for the general point API. Current theta is an exact supplied angle.
- [ ] **CURRENT-RP-PHYSICAL-POINT-API** Combine current inverse, complete source admission, true scale/source accuracy and materialization into callable u/v/w/p. Publish supported domains and actual errors. The present family + factored operator rows is progress toward this API, not its completion.
- [ ] **CURRENT-RP-INDEPENDENT-PHYSICAL-VALUES** Independently corroborate actual materialized/factored values and fixed-x/spatial derivatives with true errors and original scales. Moderate root fixtures alone cannot certify original-scale physical values.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Add actual core, transition, Rh/O2/O3, compact/quiet and leading-to-Rp source rows and seams on the same current physical graph. Existing ordinary Cartesian fixtures route before Rp and still need these providers.
- [ ] **CURRENT-PHYSICAL-AXIS-LIMIT** Construct the actual current core parity, r=0 Cartesian limits/derivatives, pressure and divergence. The new zero-axial case is still off radial axis and does not complete this task.
- [ ] **CURRENT-RP-UNIFORM-SOURCE-AND-SEAM-BOUNDS** Prove full-domain source derivatives/support and quantitative seam bounds beyond finite source cells/physical points. Propagating coordinate intervals does not itself give a uniform chart theorem.
- [ ] **CURRENT-PHYSICAL-ENERGY** Integrate real physical-volume kinetic energy with the full radial tail and lambda-dependent volume Jacobian; give controlled analytic tail and coordinate/source error.
- [ ] **CURRENT-STRESS-DERIVATIVE-ORDER-LEDGER** Identify and supply all ordinary physical source derivatives needed by signed stress and flat remainder, across current core-to-exterior layers.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Implement the original divergence-form stress and separate flat remainder with cone margins, max bounds and physical-volume L2 bounds. Do not mislabel background stress as a failed corrected residual.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement actual first-order recovery and independent five-moment repair on the shared inner interval, with divergence preserved by curl/streamfunction construction.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n-dependent recovery, independent repairs, controlled truncation, finite-order remainder and smooth summation. Coordinate rescaling is not coefficient recursion.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Implement both original families and mean corrections, admissible stress realization, averaged quadratic cancellation and independent remainder.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure actual contraction, axial/radial aspect ratio, amplitude/vorticity growth, material winding, recursive scale relations and energy on supported time-dependent physical queries.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently bound the full corrected forced NS residual and physical-volume norms after correction/recursion/source accuracy are available; retain the requested full target.
- [ ] **CURRENT-FAST-RUNTIME-BOOTSTRAP** Restore current immutable source parameter/datum state with hash and semantic family/branch validation so new machines/agents can rebuild genuine analytic source functions efficiently. Do not import saved point rows as a field, silently set acceptance flags or replay source selection/future just to create a new wrapper owner. Keep exact original N/delta/selected logC and current context.

Predecessor: [whole-box source callbacks and full backlog](CURRENT_ORIGINAL_RP_NATIVE_BOX_SOURCE_2026_10_10.md). Full goal remains ACTIVE. Next concrete delivery: source-bound original amplitude/log-scale evaluation and signed physical value/error delivery, then general/global/axis coverage.
