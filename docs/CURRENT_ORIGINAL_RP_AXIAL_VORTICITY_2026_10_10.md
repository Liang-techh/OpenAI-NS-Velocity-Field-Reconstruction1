# Original nonzero axial vorticity with exact mu — 2026-10-10

The original pulse_exit axial vorticity now has a negative, proved nonzero source-scale value and meets the requested 0.1% total physical relative-width bound. Its formerly unresolved subtraction is reduced using the actual original theta derivative law, retaining mu before enclosure. Full goal **ACTIVE / INCOMPLETE**.

Source commit: [b214209b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b214209bbe5d1b49e7767783ee602a8aa5997eff). Source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_axial_vorticity.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual result and source identity

The supported delivery remains pulse_exit, xi=21/2, Z=371/1000, finite log-time offset=1/5, theta=7/10. The preceding derivative interface exposes all 144 spatial4/time1 values and 13 linear fields. Its numeric omega_z enclosure crossed zero because Utheta_y and Utheta/2 were bounded separately.

The actual source packet uses the radial law below. This includes the current refined inlet coefficient, incoming histories and selected pulse source; radial recovery contributes no surviving term to axial curl.

```text
y = log(R)
Utheta_y = -(1/2 + mu) * Utheta
omega_z = sqrt(2/R) * lambda**(-2-delta) * (Utheta_y + Utheta/2)
        = -mu * sqrt(2/R) * lambda**(-2-delta) * Utheta
mu = exp(-4*(exp(40)+11))/1000
```

The new value is exported as `exp(L_theta_curl + logmu) * (-sqrt(2)*theta_coefficient)`. mu is an unchanged exact source function in its defining log scale. It is neither dropped nor recovered by subtracting two rounded nearly equal numbers. The physical R and lambda factors already supplied by the original Cartesian map are retained; no extra radial-velocity sqrt(2) is introduced.

- omega_z is strictly negative and nonzero at this actual delivery.
- Independent 800-digit coefficient arithmetic is included in the directed result. The recovered mu expression is included in the original broader unresolved curl enclosure after exact common-scale cancellation.
- The exact pure inverse-mu coefficient in its source log partition is zero. Its total physical width is still computed from the original centered residual scale and signed coefficient bounds; it passes 0.1% without promoting any inherited u/v width flag.
- omega_x remains negative and omega_y positive as established by the preceding linear interface. All three vorticity component signs are now known at this point.
- Ordinary absolute nonzero numeric materialization remains false. The inherited 144-row width count is still 26; the derived omega_z result is reported separately and is not added to that count.
- Time growth, material winding, energy, uniform charts, arbitrary xyz/t, stress admissibility, flat remainder, n-dependent recursion and oscillatory corrections remain open.

## Bound actual source path

The adapter requires the accepted derivative owner and its owner-issued pulse_exit field. It follows SourceProduct → RemainingPressureTail → RefinedPulseCoefficients, including the actual copied refined pulse, native raw pulse callback and mixed velocity callback. It checks the actual refined main/data/high-packet methods, dispatch owner, transport callable and original mu/delta enclosures. The defining AST binds `bp=.5+mu`, `[u*(-bp)**k for k in range(5)]`, `u` from `self.data(Z)`, and the original inlet coefficient `r*k['U']`.

From the actual typed Cartesian rows, exactly two Ur primitive operators cancel first. Only the base and first-log-radial Utheta primitives may survive; their source units, exact powers and complete log scales must agree. Symbolic Cartesian operator algebra then proves the -mu*sqrt(2) factor. Reports and opaque values cannot hydrate a new owner, change source parameters or substitute a different callback.

The focused checker independently derives Cartesian-to-cylindrical axial curl and differentiates the physical pure-power swirl, verifies the actual source coefficient/sign/scales/units and relative-width result, checks inclusion in the older unresolved interval, and rejects copied values/fields, a substituted refined source method, changed original mu and modified curl proof. The accepted constructor and public value/report calls pass. One existing read-only reviewer is reused: **GPT-5.6 Luna / max**; no additional or Astra worker is spawned.

## Executable use

Construct the live accepted `arithmetic`, `differential` and `field` from [the derivative checkpoint](CURRENT_ORIGINAL_RP_DIFFERENTIAL_FIELD_2026_10_10.md#executable-use). Append:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_axial_vorticity import CurrentOriginalRpAxialVorticity as AxialVorticity

axial = AxialVorticity(differential)
omega_z = axial.evaluate(field)
omega_z_report = axial.report(omega_z)
assert omega_z_report['value']['sign'] == -1
assert omega_z_report['value']['original_absolute_width_target_satisfied']
omega_z_over_itself = arithmetic.ratio(omega_z, omega_z)
```

The source law is algebraically valid in its original off-axis pulse domain; this adapter admits only current pulse_exit deliveries and this receipt independently exercises only the supported point. Additional queries require their own actual inverse/source rows and checks. JSON is a report, not a serialized live owner. Prepare fresh original queries before freezing derived cache views; never reset failed source guards.

## Completed bounded tasks

- [x] **CURRENT-RP-EXACT-AXIAL-VORTICITY-LAW-POINT** Bind the actual refined theta radial derivative source, cancel original Ur operators and preserve exact mu plus original R/lambda factors before enclosure.
- [x] **CURRENT-RP-NONZERO-AXIAL-CURL-POINT** Restore negative nonzero omega_z at the supported delivery, independently check its signed interval/source partition and 0.1% total width, and expose a live source-scale value.

## Next executable tasks

- [ ] **CURRENT-RP-SIGNED-BACKGROUND-SOURCE-OPERATORS** Feed actual current-Rp raw similarity rows and shared five moments/P0 into the original pure stress/remainder operator bodies in `pulse_end_physical_C2.py` and `current_pre_pulse_stress_operator.py`. Reuse `current_pulse_main_exit_background_tensor.py` for pulse-exit incoming/cross-term structure. Bind all AST assignments and original viscosity/scales. Acceptance: independently checked signed momentum decomposition into divergence-form stress and separate remainder at this point; preserve paired source products before interval evaluation. Cone/flatness remain false until bounded.
- [ ] **CURRENT-RP-OFFCENTER-AND-TIME-VORTICITY** Prepare distinct nonzero Z/pulse sections and several exact log-time offsets before frozen views. Check complete derivative fields and exact axial-curl law per delivery; implement compatible cross-query scale ratios with conflicting source bindings rejected. Acceptance: measured radial/axial scale ratios and vorticity growth exponents, not only known signs at one time.
- [ ] **CURRENT-PHYSICAL-MATERIAL-WINDING-ENERGY** Compute physical trajectories and cumulative material winding separately from instant streamline geometry; preserve coordinate Jacobians, time scales and tails in energy/volume integration. Acceptance: actual bounded observables at several times, finite-energy tail control and explicit sampled/uniform scope.
- [ ] **CURRENT-RP-PRESSURE-END-AND-ALL-CHART-FIELD** Resolve the xi=13 late-tail lower-bound guard, then install controlled postpulse/flatten/steep/waiting/collar/exact-heat derivative fields with actual original source definitions. Keep point, seam, axis and uniform claims separate.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t, solve the original inverse, retain support joins, pressure compatibility and axis parity/Cartesian limits. Current source-law and point APIs do not satisfy this task.
- [ ] **CURRENT-SIGNED-STRESS-CONE-FLAT** After source operator admission, bound signed regional cone margins and separate flat-remainder maxima/physical L2 norms at inner exit, annulus, pulse/end, flatten and collar. A numerical decomposition alone is insufficient.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement the different n=1 and n>=2 equations on a common original inner domain, independent moment repair, finite-order remainder and smooth sum. Truncate streamfunction/vector potential before curl. Acceptance requires actual higher-order coefficients and n-dependent identities; self-similar coordinate rescaling alone is insufficient.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-FULL-NS** Implement both pulse families and mean corrections, independently verify averaged quadratic stress cancellation, then check fixed-forcing corrected Cartesian NS residual/max/L2 norms.

Previous checkpoint: [original derivative and correlated linear fields](CURRENT_ORIGINAL_RP_DIFFERENTIAL_FIELD_2026_10_10.md). The next implementation priority is the signed original background stress/remainder source adapter, alongside controlled time/query coverage. Mark only source-demonstrated scope complete.
