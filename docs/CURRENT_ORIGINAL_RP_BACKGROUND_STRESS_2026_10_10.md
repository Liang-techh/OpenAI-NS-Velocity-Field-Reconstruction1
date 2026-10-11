# Original full background stress and separate remainder — 2026-10-10

The actual original pulse_exit source now supplies a completed symmetric stress tensor, its divergence, a separate remainder and the full background momentum residual in cylindrical and Cartesian components. All 30 output signs/zeros are bounded: 5 exact zero outputs and 25 strictly signed nonzero outputs. Independent Cartesian NS source assembly agrees with the decomposition. Full goal **ACTIVE / INCOMPLETE**.

Source commit: [0567e214](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/0567e2147677251ef96a6acd48574b42fa3be29e). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_background_stress.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual delivery and result

The checked delivery remains pulse_exit, xi=21/2, Z=371/1000, finite log-time offset=1/5, theta=7/10. Original viscosity is 1; the existing original similarity/physical map and actual current source parameters are retained. The incoming cumulative moments, nonlinear cross histories, independent P0 and remaining absolute pressure are all consumed. No surrogate velocity, pressure or parameter family is introduced.

| Output | Signs in component order |
| --- | --- |
| cylindrical_stress | r_theta: +, r_z: +, theta_theta: +, rr: 0, theta_z: 0, zz: 0 |
| cylindrical_divergence | radial: 0, theta: +, axial: + |
| cylindrical_remainder | radial: +, theta: +, axial: + |
| cylindrical_momentum_residual | radial: +, theta: −, axial: − |
| Cartesian_stress | xx: +, xy: −, xz: +, yy: +, yz: +, zz: 0 |
| Cartesian_divergence | x: −, y: +, z: + |
| Cartesian_remainder | x: +, y: +, z: + |
| Cartesian_momentum_residual | x: +, y: +, z: − |

The exact radial stress divergence is zero because the completed diagonal is retained. The other exact zeros are T_rr, T_theta_z, T_zz and Cartesian T_zz. The residual is nonzero, as expected at the background stage. Its signs do not establish a corrected NS solution or a residual smaller than 1e-3.

The focused independent checker included 377 collected shared-source product coefficients and 25 nonzero output bounds with directed 800-digit arithmetic. No output has an unresolved scale ratio. **None of the 25 nonzero outputs yet meets the total physical 0.1% relative-width target.** Exact signs and source-scale enclosures are available; ordinary nonzero absolute materialization remains false. The inherited 144 derivative rows and their 26 total-width passes are unchanged.

## Original decomposition and source binding

```text
R_B = partial_t u_B + (u_B dot grad)u_B + grad p_B - Laplacian u_B
    = -div(T_B) + E_B

T_rtheta = lambda**(-2-delta) * (Itheta + Stheta)
T_rz     = lambda**(-2-delta) * (Iz + Sz)
T_thetatheta = r * partial_z(T_rz)
T_rr = T_thetaz = T_zz = 0; T is symmetric
```

The pure `original_full_stress` AST in `pulse_end_physical_C2.py` is replayed with the actual full current-Rp source functions and histories. Source sectors carry their exact lambda and R powers. Shared primitive monomials and exact cumulative FTC identities are collected before interval multiplication; the already admitted Utheta_y=-(1/2+mu)Utheta law preserves original mu. Pressure differentiation replaces only positive radial orders through P_y=Utheta^2/2; independent P0 and total pressure at radial order zero remain intact.

The runtime radial velocity is explicitly bound to `CompliantPulseHighJets.radial`, the original normalized recovery numerator and denominator, and its physical product rows in `transport_mixed`. The paper 3.9 numerator, 3.8 log-radial derivative and physical axial product-rule identities are required. The actual refined pulse/dispatch/raw/mixed callback chain is inherited from the accepted axial owner and checked again. Each delivery must carry its original forward source packet and physical-radial-product evidence. Source products, exact operators and their scale DAGs are sealed against mutation.

The independent checker constructs div(T) directly in Cartesian coordinates, then assembles the three NS components from the actual typed derivative rows: time derivatives, quadratic velocity-gradient products, pressure gradients and velocity Laplacians. It proves agreement after the bound original recovery/FTC identities; it also separately checks the signed product and final output intervals. Copied packets, invalid targets, substituted radial callbacks, missing actual packets, changed defining DAGs and a removed tensor diagonal are rejected. The accepted constructor and public evaluate/report calls pass. The existing read-only reviewer is reused: **GPT-5.6 Luna / max**; no additional worker is spawned.

## Executable use

Create the live accepted derivative and axial owners using [the axial checkpoint](CURRENT_ORIGINAL_RP_AXIAL_VORTICITY_2026_10_10.md#executable-use). Append:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_background_stress import CurrentOriginalRpBackgroundStress as BackgroundStress

background = BackgroundStress(axial)
packet = background.evaluate(field, '1/1000')
report = background.report(packet)
assert report['cylindrical_divergence']['radial']['exact_zero_enclosure']
assert report['cylindrical_stress']['r_theta']['signed_log_value']['sign'] == 1
assert report['convention'] == 'R_B=-div(T_B)+E_B'
```

JSON is a report, not a hydrated owner. Use the fresh-process typed chain referenced by the earlier checkpoints. Prepare distinct original deliveries before freezing derived views; never reset failed source guards. The adapter admits only the original pulse_exit source chart, and this receipt independently exercises only the supported point.

## Completed bounded tasks

- [x] **CURRENT-RP-SIGNED-BACKGROUND-SOURCE-OPERATORS** Replay the full original stress and remainder bodies with actual incoming histories, independent P0, total pressure, source viscosity/scales and physical derivative rows. Independently verify the completed tensor and all three Cartesian NS source identities at the supported point.
- [x] **CURRENT-RP-RADIAL-RECOVERY-SOURCE-BINDING** Connect the live radial recovery callback, its exact AST and paper derivative identities to the actual current physical product rows used in the remainder.
- [x] **CURRENT-RP-SIGNED-BACKGROUND-POINT-OUTPUTS** Return all 30 cylindrical/Cartesian outputs with 377 independently checked collected source products, 5 exact zeros and 25 resolved nonzero signs. Keep physical accuracy, admissibility, flatness and recursion gates separate.

## Next executable tasks

- [ ] **CURRENT-RP-SIGNED-CONE-MARGINS-POINT** Locate the original admissible stress inequalities for this pulse sector. Build the actual signed relative combinations from T_rtheta, T_rz and the completed diagonal using common primitive products/scales before enclosure. Include incoming history and pressure contributions. Report each lower margin, its source identity, the applicable sector and whether its inequality passes. Individual tensor signs do not satisfy this task.
- [ ] **CURRENT-RP-BACKGROUND-RELATIVE-WIDTH** Diagnose coefficient width and centered log-scale width separately for the 25 nonzero outputs. Refine only the dominant original source dependencies, preserving exact mu/delta/logC and incoming histories. Acceptance: independently included narrower enclosures and explicit 0.1% pass/fail counts; never reset a width flag or replace an original source constant.
- [ ] **CURRENT-RP-FLAT-REMAINDER-SCALE-CANCELLATION** Partition E_r, E_theta and E_z into their original exact lambda/R/source-scale sectors. Compare them to the paper reference decay scales with common-source cancellation, retaining time offsets and Jacobians. Distinguish a point relative remainder bound, finite-order decay and uniform flatness. Acceptance: actual controlled relative decay evidence; one small or positive remainder value is insufficient.
- [ ] **CURRENT-RP-OFFCENTER-AND-TIME-BACKGROUND** Prepare distinct nonzero Z/pulse sections and exact log-time offsets before freezing owners. Reuse the original inverse and source algorithms, evaluate derivative/axial/stress packets and reject conflicting cross-query bindings. Acceptance: measured radial/axial width ratios and velocity/vorticity/stress/remainder scale exponents at multiple deliveries.
- [ ] **CURRENT-RP-PRESSURE-END-AND-ALL-CHART-FIELD** Resolve the xi=13 late-tail lower-bound guard, then install source-bound postpulse, flatten, steep, waiting, collar and exact-heat derivative/stress packets. Include actual source joins and pressure datum compatibility. Mark sampled point, seam and uniform bounds separately.
- [ ] **CURRENT-SIGNED-STRESS-CONE-FLAT-REGIONS** Extend accepted point margin/remainder machinery over inner exit, matching annuli, pulse/end, flatten and collar. Bound the appropriate signed cone margins and physical stress/remainder maxima and volume-L2 norms with actual coordinate Jacobians. Every claimed region needs controlled coverage and errors.
- [ ] **CURRENT-PHYSICAL-MATERIAL-WINDING-ENERGY** Compute physical trajectories and cumulative material winding at several times separately from instantaneous streamlines. Integrate energy over the actual physical domain, including certified radial tails. Acceptance: bounded time-dependent observables and explicit sampled/uniform scope.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t, solve the original inverse, preserve all support joins and pressure compatibility, and implement axis parity/Cartesian limits. The current off-axis point adapter does not establish a globally available velocity/stress field.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** After leading moments and required stress conditions hold, implement distinct n=1 and n>=2 radial recovery equations on a common core domain. Keep independent moment repair, finite-order remainder and smooth summation. Truncate streamfunction/vector potential before curl. Acceptance requires actual order-dependent coefficients and identities; coordinate scaling is insufficient.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-FULL-NS** Implement both oscillatory pulse families and mean corrections, verify realizable averaged quadratic stress cancellation, then evaluate fixed-forcing corrected Cartesian residual maxima and physical L2 norms against 1e-3. Do not choose forcing as the measured residual.

Next implementation priority: source-correlated admissible-cone margins and relative remainder scales, followed by controlled multitime/all-chart coverage. Cone admissibility, flat remainder, uniform/global/axis claims, actual n-dependent recursion and oscillatory correction remain open.
