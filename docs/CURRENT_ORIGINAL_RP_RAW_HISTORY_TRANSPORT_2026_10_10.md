# Current factorized raw postpulse histories and radial transport (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [80e83353](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/80e83353c21e2713aaabe1016b3bdb0490e29073) connects the accepted absolute source geometry to actual five-history outputs on all nine postpulse charts. It supplies callable factorized raw source enclosures and first native radial derivatives. Unrestricted numerical absolute point values, absolute heat terminal closure, all-chart mixed contracts, Cartesian u/v/w and temporal recursion remain open.

## Result

`CurrentOriginalRpRawHistoryTransport` consumes an accepted `CurrentOriginalRpSegmentedRadius`. It keeps the same current postpulse providers, selected pulse/flatten, exact original inlet functions, independent P0, complete future/Gamma tails and unique selected repair. Only the radius expression view is copied. No old prefix constructor or second repair solve is used when an existing accepted owner is supplied.

The current exact inlet u at Z=0 defines U0. The positive formal terminal velocity scale is

```
Ev0 = Pstar * U0 * exp(-13/(2mu)-13)
Utheta = Ev0 * theta
```

`theta` is the actual packet's `theta_over_Ev0_Taylor` from flatten, power/angular, steep entry/power/exit, waiting, collar or full Gamma exterior. The numerical Ev2 cap never defines Ev0.

| Output | Absolute scale | Live directed axial Taylor coefficient |
|---|---|---|
| Mz | R Ev0 | theta times inherited normalized Mz |
| Mtheta | R^(3/2) Ev0 | sqrt(2) theta times normalized angular history |
| Mtheta_z | R^(3/2) Ev0^2 | sqrt(2) theta^2 times inherited normalized mixed history |
| Mztheta | R Ev0^2 | theta^2 times complete energy history |
| Mp | Pstar^2 | original P-P0 normalized by Pstar^2 |
| P0 | Pstar^2 | independent analytic axis datum |
| pressure | Pstar^2 | actual inherited total pressure |

Mz and Mtheta_z vanish only by the admitted current compact terminal source theorem. Angular, energy and pressure memory are not reset at seams. The actual heat angular defect and pressure-infinity constants remain present; no overlap containing zero is promoted to a terminal identity.

## Exact large-scale cancellation

A `FactorizedSourceTaylor` represents a true positive scale function times six directed axial coefficients, ordered as dZ^j/j!. For units R^r Ev0^e Pstar^p, the scale retains separate exact parts:

```
r logRp + r local_offset
+ 13(r-e/2)/mu
+ (e+p) logPstar + e logU0 - 13e
```

The R Ev0^2 energy scale cancels the common inverse-mu pulse term exactly before arithmetic. The R^(3/2) Ev0 angular scale retains 13/mu. A short offset is never added numerically to the enormous origin. Signed coefficient enclosures and directed width bounds remain in these explicit units. Neither a midpoint nor an interval endpoint supplies a point value.

The scales are Z independent, so each axial coefficient receives the same factor. These are finite derivative enclosures of source functions, not polynomial fits or unrestricted numeric physical points.

## Actual first native radial transport

With J=dlogR/d(native coordinate), the original cumulative equations after Uz=0 give

```
d_native Mz = 0
d_native Mtheta = J sqrt(2) R^(3/2) Ev0 theta
d_native Mtheta_z = 0
d_native Mztheta = -J R Ev0^2 theta^2 / 2
d_native Mp = J Ev0^2 theta^2 / 2
d_native P0 = 0
d_native pressure = d_native Mp
```

J comes from the accepted true coordinate function and its directed bound. It is applied once. Derivative output rows are conservative source enclosures, not exact numerical Taylor point coefficients. In particular, Mp/pressure derivative units are Ev0^2, while their history units are Pstar^2; consumers must not silently discard this difference.

## API and evidence

Reuse the accepted owner from the ongoing construction:

```python
owner = CurrentOriginalRpRawHistoryTransport(accepted_radius_owner)
view = owner.evaluate('waiting', '.449', '1/2')
raw = view['raw_histories']['Mtheta']
raw.report()  # exact scale function nodes, directed coefficients and widths
view['first_native_radial_derivatives']['Mp'].report()
```

- Twenty actual current source calls include active flatten/angular regions and every postpulse seam; all nine charts are covered.
- 1440 live coefficient identities and 960 directed scaled width checks.
- 72 independent scale-function identities, including nine exact energy-scale pulse cancellations.
- Five signed radial densities read from the original general moment implementation; derivatives independently derived from cumulative integrals with retained nonzero seed memory.
- 336 axial seam coefficient overlaps plus exact scale identities across eight seams. Overlaps remain diagnostics, not high mixed or terminal closure proofs.
- An accepted runtime and a fresh rational waiting call succeed. All 1290 dependency bytes match the Git index.
- Scoped static review: **GPT-5.6 Luna / max — PASS**. The native-Jacobian enclosure distinction is retained in output and documented here.

Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_raw_history_transport.py`, `.json.gz`, `_check.py`, `_check.json`.

## Detailed next tasks

- [x] **CURRENT-RP-FACTORIZED-POSTPULSE-FIVE-HISTORIES** Connect all nine actual postpulse charts to exact R/Ev0/Pstar scale functions, retain five histories plus P0/pressure, and expose directed scaled coefficients/error widths.
- [x] **CURRENT-RP-FIRST-NATIVE-HISTORY-TRANSPORT** Apply the actual coordinate Jacobian once to the five cumulative equations, with correct signed energy/pressure half-factors and derivative units.
- [ ] **CURRENT-RP-SAME-REPAIR-PRESSURE-CLOSURE** Inject the existing `owner.post.before.seed.exact.repair` into downstream angular/balance/raw/pressure closure. Rebind current native Xp/U/H/Pin, literal waiting/radius/amplitude functions, roots, weights and callbacks. Do not solve a second repair branch. Retain actual terminal constants until functional identities prove their elimination.
- [ ] **CURRENT-RP-ABSOLUTE-HEAT-CLOSURE** Match current absolute heat Rtail/amplitude and complete collar/Gamma tails to independent P0. Prove Dtheta/Cp as functions, separately from relative repair equations and interval overlaps.
- [ ] **CURRENT-RP-POSTPULSE-FIVE-HISTORIES** Finish numerical absolute/raw physical point evaluation and validated error bounds using the new factorized transport. Keep this broad task open until actual point values are available; function handles and scaled enclosures alone do not finish it.
- [ ] **CURRENT-RP-RAW-PULSE-HISTORIES** Extend the common factorization to the six active pulse charts without inheriting post-Rv zero theorems. Preserve nonzero Mz/Mtheta_z, all incoming memory and selected phase dependence.
- [ ] **CURRENT-RP-RAW-HIGH-MIXED-TRANSPORT** Extend actual radial/axial derivatives to the required joint orders and every seam; differentiate amplitude and history factors by product rule. Use the same physical units before comparing seam rows.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Carry signed error/remainder bounds through the lazy absolute factors and pressure cancellation. Define error contracts usable by the numeric point oracle without choosing cap endpoints.
- [ ] **CURRENT-RP-UNIFORM-SELECTED-SEAMS** Complete mixed-y/Z contracts across six pulse charts, entire flatten and all postpulse charts. Current geometry identities and coefficient overlap diagnostics do not close this task.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Join the actual core, transition, Rh, O2/O3, compact repair, quiet and current exterior under one derivative contract.
- [ ] **CURRENT-FREQUENCY-REPRESENTATION** Make exact selected N=2^(2^J), J=1358356628656378313 phase executable lazily; never substitute a smaller frequency.
- [ ] **CURRENT-INVERSE-PICARD-ORACLE** Implement same-source monotone inverse, flat branches, Picard controls and signed integrals with quantified error.
- [ ] **CURRENT-NUMERIC-POINT-ORACLE** Evaluate defining functions on valid native cells with phase limits and error. Caps/recipes remain bounds and never become point values.
- [ ] **CURRENT-PHYSICAL-TIME-MAPPING** Complete original similarity-to-physical coordinates/time and axis limits from the same coefficients.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] and pressure with derivative/error contracts.
- [ ] **CURRENT-STRESS-DECOMPOSITION** Restore signed divergence-form stress and separate flat remainder with all pressure and inertial terms.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Quantify global cone margins and flat decay through exits/end supports/collar.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement n=1 recovery on the common inner domain with an independent five-moment repair.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n>=2 recovery, order-dependent repair, curl-preserving truncation and controlled summation.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Build both original oscillatory families and mean corrections from admitted stress/recursive coefficients.
- [ ] **CURRENT-AVERAGED-STRESS-CANCELLATION** Prove averaged quadratic flux cancellation and bound remaining oscillatory/flat terms.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure contraction, aspect ratio, velocity/vorticity growth, material winding, scale recurrence and finite energy from computed fields.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently evaluate final corrected forced NS residual and L-infinity/volume-L2 targets.

Agents: implement the next unmet dependency, reuse accepted unchanged mathematics, publish evidence for changed boundaries, and mark only scoped achieved tasks complete with their source commit. Preserve the full ACTIVE goal and every unproved physical/global gate.

Predecessor: [absolute segmented source radii](CURRENT_ORIGINAL_RP_SEGMENTED_RADIUS_2026_10_10.md).
