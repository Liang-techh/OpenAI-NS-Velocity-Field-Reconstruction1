# Update: current native affine locator installed (2026-10-10)

[Source ca93c2e7](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ca93c2e7ea3a7e4cb7ed67b34664330d6f4dbb0a) completes the bounded fifteen-chart affine radial locator and source-bound logRp evaluator. Read the [new handoff](CURRENT_ORIGINAL_RP_NATIVE_CHART_LOCATOR_2026_10_10.md). Native interval source delivery and seam calls remain open, so the broader caller tasks below are not marked complete. Numerical u/v/w/p and global/axis coverage remain open.

---

# Current original Rp directed physical inverse (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. [Source 210c14c9](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/210c14c91e8a40b51499e18d61037a773e209bf5) adds the physical-coordinate inverse on the accepted current original Rp function graph. Exact Cartesian x/y/z/t or logarithmic cylindrical inputs return exact inverse function nodes and directed log-lambda, Z, logR and angle enclosures. The original delta, finite correction N, analytic P0 and complete source history remain unchanged. This completes the bounded coordinate inversion step. Native chart admission, source point-error delivery and numerical u/v/w/p remain open.

## Current public interface

```python
owner = CurrentOriginalRpPhysicalInverse(accepted_physical_source_map)
view = owner.cartesian('0.8', '-0.3', '0.2', '0.95')
view['coordinate_functions']['log_lambda']  # exact unique-root FunctionRef
view['coordinate_functions']['Z']
view['source_logR_enclosure']
view['directed_inverse_mapping']['actual_log_lambda']
view['theta_enclosure']
view['inverse_log_lambda_absolute_width']
view['strict_numerical_Z_interior_resolved']

near = owner.log_cylindrical(log_r='0', log_abs_z=None,
    sign_z=0, log_tau='-1e1000', theta='7/10')
```

Reuse the accepted current runtime and source map; do not replay unchanged selection, five-moment repair, postpulse future or heat solves. The new class is `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_physical_inverse.py`. Its compressed candidate, focused checker and check receipt have the same stem. No old tensor registry or 33-region locator owner is constructed.

Cartesian scalar strings, integers and Fractions are exact rational physical inputs. A finite binary float denotes its exact binary value: 0.1 is not silently replaced with decimal 1/10. Intervals and FunctionRefs are not accepted as scalar physical point inputs. The log interface takes exact finite log_r, log_tau, angle and, for nonzero z, log_abs_z with explicit sign_z=+1/-1. Zero z requires sign_z=0 and log_abs_z=None. It never constructs a numerical exp(log_abs_z) or absolute similarity radius.

The current inverse domain is **r>0, t<1, nu=1, actual original 0<delta<1**. All finite axial physical coordinates have a unique inverse with exact |Z|<1. Axis points require the current core limit and are rejected here. The unique inverse theorem does not establish that logR belongs to the Rp-to-heat source region. `current_original_Rp_native_chart_location_installed` and `current_original_Rp_inverse_intervals_consumed_by_point_source` remain false.

## Exact inverse and directed solver

For q=log(lambda), the current forward source is

```
lambda^2*(1-Z^2)=tau, tau=1-t>0
z=Z*lambda^(1-delta), r=lambda*sqrt(2R)
2q=logaddexp(log_tau, 2log_abs_z+2delta*q)
Z=sign(z)*exp(log_abs_z-(1-delta)*q)
logR=2log_r-log(2)-2q
```

At z=0, q=log_tau/2 and Z=0 exactly. In the nonzero branch the implicit function has derivative between 2(1-delta) and 2, giving a strictly increasing unique root. The exact node retains the actual graph delta=exp(-4logPstar-30), rather than a numerical enclosure endpoint or a replacement zero. The numerical solver consumes the directed actual delta enclosure separately.

The Cartesian path calls the previously accepted pure directed implicit solver. The log-z adapter copies that solver's original nonzero root iteration, brackets, derivative enclosure, interval Newton contraction and Z-complement calculation. Only the explicit log_abs_z input and exact nu=1 prefix are added. AST bindings and dependency hashes bind this copy to the current original source. The current physical forward coordinates method is separately AST-bound to its nu=1 radial/axial assignments and full method hash.

Numerical exp bounds for very negative arguments remain enclosure bounds. They do not define the exact coordinate function. Directed stopping tolerance is relative to the scale of q, and the actual absolute q-width is reported separately. A very large log input can therefore leave a wide absolute uncertainty. Near |Z|=1, numeric bounds may touch an endpoint even though the exact finite inverse is strictly interior. `strict_numerical_Z_interior_resolved=False` keeps this unresolved numerical admission visible; no endpoint source evaluation is performed.

The assembled result requires a private live solve witness from the same inverse owner. Mapping/ref dictionary identities and a recomputed packed fingerprint bind the result to its actual inputs, source graph and directed solve. An unwitnessed `assemble()` call is rejected. This protects the current call path; it is not a general hostile-process import API. The exact source map's downstream interval/source admission still has to be implemented.

## Focused evidence

- Nine actual current input observations: five Cartesian quadrants/axial signs including zero z and the negative-x angle cut, plus four extreme log-time/log-axial cases with magnitude 1e1000.
- Seven exact nonzero inverse function nodes and 89 exact coordinate function identities, with exact delta and unit viscosity retained.
- Ten independent reference roots at the actual delta enclosure corners for the five Cartesian cases. Directed logR, Z and angle contain the reference results; the current root widths are below 1e-50 for those cases. The corner values check enclosure containment only and never select the defining source delta.
- Three independent moderate log-z fixtures agree with the original direct-z root program. These are explicitly fixtures, not substitutes for current coefficients.
- Twelve invalid-source/domain rejections, including axis/terminal time, nonfinite/interval inputs, inconsistent axial logs/signs, foreign graph, defining-delta replacement and an unwitnessed solve import. Binary float semantics are checked exactly.
- An accepted constructor and fresh Cartesian and extreme zero-z log queries pass on the current runtime. The existing sole read-only worker is **GPT-5.6 Luna / max**; it reports no material false acceptance in the witness/source-binding changes. Dependency bytes are compared with the staged Git blobs before upload.

## Completed bounded tasks

- [x] **CURRENT-RP-DIRECTED-PHYSICAL-INVERSE** Actual Cartesian/log cylindrical inverse with exact function definitions and directed coordinate enclosures on the current accepted graph.
- [x] **CURRENT-RP-LOG-AXIAL-INVERSE** Preserve the original pure root loop without numerical materialization of extreme physical z, lambda or R.
- [x] **CURRENT-RP-INVERSE-SOURCE-BINDING** Same graph/actual delta, nu=1 current forward AST binding, original program AST binding, current receipt hashes and live solver witness.
- [x] **CURRENT-RP-INVERSE-INPUT-SEMANTICS** Exact rational/binary input distinction, zero/sign/axis/time domain and atan2 branch.

## Next implementation tasks and acceptance

These tasks extend the accepted operator; each new agent must update its checkbox only when its stated acceptance is achieved. Inherited physical spatial4/time1 rows are already complete at the chart-parametric scope described in the predecessor handoff.

- [ ] **CURRENT-RP-PHYSICAL-INVERSE-CONTRACT** Finish the broader inverse-plus-source-domain contract. The directed inverse subtask above is complete, but source-domain admission and interval consumption are not. Admit only native chart locations established from exact current radius functions and propagated inverse errors. Keep outside-Rp, uncertain boundary and endpoint-touching Z results explicit.
- [ ] **CURRENT-RP-NATIVE-CHART-LOCATION** Partial: all fifteen exact affine inverse maps, directed radial chart candidates, current source-bound logRp and exact common-factor cancellation are installed. The broader task remains open for native interval source callbacks and accepted seam provider calls. See the latest native locator handoff; do not repeat the completed affine implementation.
- [ ] **CURRENT-RP-NATIVE-INVERSE-MAPS** For pulse s, stage a/b/c coordinates, integral-defined W, waiting, heat collar and heat exterior, invert the actual native logR expression with the correct Lrel-4, Ts or W Jacobian. Do not return an arbitrary rational representative. Preserve unknown/overlapping admission intervals for refinement.
- [ ] **CURRENT-RP-NATIVE-BOUNDARY-CALLERS** Route an exact inherited boundary through the accepted fourteen same-function seam callers. Convert left/right units explicitly, including full mixed ordinary derivative conventions. A directed interval crossing a seam requires both source enclosures or further narrowing; it cannot be rounded into one chart.
- [ ] **CURRENT-RP-INTERVAL-SOURCE-INPUT** Extend actual current source callbacks beyond their present exact rational Z/native inputs. Accept directed coordinate boxes tied to exact inverse FunctionRefs; enclose their dependence with actual derivative/phase bounds. Do not substitute interval midpoint, archived receipt rows or independently fabricated jets.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Deliver a real error budget in each amplitude unit for phase inversion, signed quadrature, narrow beta supports, coefficient recovery and Gamma tail. Include coordinate uncertainty and report absolute/relative error of each physical scale group. Existing directed jets alone do not establish sufficient point accuracy.
- [ ] **CURRENT-RP-EXACT-POINT-GRAPH-EVALUATOR** Add admitted evaluation of exact_original_physical_log_lambda_inverse and exact_physical_atan2 nodes, together with existing sin/cos and signed integer-power nodes. Bind numerical results to the current exact input/source graph. Combine correlated log factors before exponentiation and report unresolved magnitude/precision honestly.
- [ ] **CURRENT-RP-PHYSICAL-ROW-MATERIALIZER** Consume admitted source/error rows and exact grouped scales to produce directed finite physical values. Retain signed cancellations, pressure baseline/radial unit split and independent P0. No practical replacement frequency or smaller N is allowed.
- [ ] **CURRENT-RP-PHYSICAL-POINT-API** Return u/v/w/p plus errors/provenance for requested supported (x,y,z,t), combining inverse, chart location, interval callback and materializer. Exercise pulse, flatten/angular, waiting/collar and exterior points, with genuine source accuracy and inherited interfaces.
- [ ] **CURRENT-RP-PHYSICAL-POINT-INDEPENDENT-CHECK** Independently compare completed point queries and original fixed-x time/spatial derivatives at bounded representative points. Include angle cut, axial signs and correlated near-critical time. Coordinate fixtures alone do not validate source velocity.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Add actual current core, transition, Rh/O2/O3, compact/quiet and leading-to-Rp source rows/seams to this same physical graph. The fifteen outer charts do not establish global field coverage.
- [ ] **CURRENT-PHYSICAL-AXIS-LIMIT** Construct actual current core parity/limits and Cartesian derivatives at r=0; check pressure/divergence independently.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver the full callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] after point/error and global/axis admission. Publish supported physical domain and actual errors.
- [ ] **CURRENT-RP-POSTPULSE-QUANTITATIVE-SEAMS** Establish full-domain mixed derivative and seam bounds, beyond exact function joins and finite diagnostics.
- [ ] **CURRENT-RP-PULSE-UNIFORM-BOUNDS** Establish actual active-pulse support/integral bounds with original amplitudes and physical units.
- [ ] **CURRENT-RP-STRESS-DERIVATIVE-ORDER-LEDGER** Identify the actual mixed space/time orders needed for stress/divergence/flat remainder and compute only missing necessary orders.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Construct signed original stress and independent flat remainder with full-domain cone margins and max/physical-volume L2 estimates.
- [ ] **CURRENT-PHYSICAL-ENERGY** Integrate actual physical-volume kinetic energy and full tail with the lambda-dependent Jacobian; establish required finite-energy bounds.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement the actual n=1 coefficient recovery equations and independent five-moment repair, then map those coefficients physically.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n>=2 recovery/repair, curl-preserving truncation, finite-order error and controlled smooth summation.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Construct both original pulse families and mean corrections from admitted stress/coefficients; establish averaged quadratic cancellation and remainder bounds.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure actual radial/axial widths, aspect ratio, growth, material winding, coefficient recurrence and energy from the admitted physical field. Coordinate scaling is not coefficient recurrence.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently validate the full corrected forced NS max and physical-volume L2 residual after stress/recursion/oscillatory corrections.

Predecessor: [current original physical Cartesian/time source map](CURRENT_ORIGINAL_RP_PHYSICAL_SOURCE_MAP_2026_10_10.md). Next agent implements native chart location and interval/error source delivery. Keep the full reconstruction ACTIVE.
