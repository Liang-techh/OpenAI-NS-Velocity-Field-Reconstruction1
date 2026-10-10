# Current fifteen-chart raw mixed transport (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [69692731](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/6969273120a8fcbca390bb15d3bf3f0338fe1908) supplies ordinary mixed log similarity R / axial Z derivatives through total order four for the current velocities, pressure and all five raw cumulative histories on all fifteen actual pulse/postpulse charts. Native-coordinate derivatives carry exact Jacobian powers as separate lazy scale functions. This completes the bounded all-chart mixed source interface; uniform/global seams, numerical physical u/v/w, stress/cone/remainder, true temporal recursion and oscillatory cancellation remain open.

## Actual current interface

```python
owner = CurrentOriginalRpMixedTransport(accepted_Rv_common_unit_owner)
view = owner.evaluate('pulse_main', '.521', 1)
view['log_radius_mixed_rows']['Mztheta']['y2_Z1'].report()
view['native_coordinate_mixed_rows']['Utheta']['n1_Z3'].report()
heat = owner.evaluate('heat_exterior', '.549', '9/2')
heat['log_radius_mixed_rows']['pressure']['y2_Z1'].report()
```

The same accepted pulse/closed-heat caller, selected C5 source, beta kernels, repair, complete future, analytic P0 and expression graph are reused. No source selection or repair is replayed. Six pulse and nine postpulse routes retain ten outputs (five histories, P0, full pressure and three similarity velocity components), each with fifteen derivative indices K+J<=4. There are 2250 logR/Z rows and 2250 native/Z rows.

Every row consists of a directed signed ordinary derivative enclosure and an explicit lazy positive scale. A scalar enclosure is not a selected numerical point or a polynomial approximation. Original factorized values and their full admitted source jets remain accessible separately.

## Velocity sources and unit guards

- Six pulse routes consume actual `physical_velocity_and_pressure_y_derivative_Taylor` velocity rows. Those rows already differentiate the physical radial prefactors, including the sqrt(R) factor in Ur; the adapter does not differentiate them a second time.
- Flatten, power/angular, steep and waiting consume the original actual Ev0-normalized velocity ODE rows. The retained legacy labels contain `over_Pstar`; mandatory producer metadata establishes the real Ev0 units. An intentionally Pstar-normalized packet is rejected.
- Heat collar and exterior use the same closed owner's full Gamma K derivative rows and current `theta_base*exp(-bh*t)` amplitude. The amplitude source is checked against the pressure closure's published `exact_pressure_scale` object. The closed packet does not publish the entire terminal object.
- Published Ur follows the preceding raw velocity convention: the provider's `Ur_over_sqrt_R_over_2_*` row is divided by sqrt(2), with sqrt(R) retained in the exact scale. Generic recovery's normalized Q=Ur/(S*sqrt(R/2)) uses the un-divided coefficient. Downstream generic field recovery must apply that conversion explicitly.

The full infinite Gamma source is retained, rather than defined by a finite radius series. Closed base coefficients and forward subtraction diagnostics retain distinct meanings.

## General five-history mixed recovery

For true actual swirl E and axial velocity V, the cumulative equations are

```
d_y Mz       = R V
d_y Mtheta   = sqrt(2) R^(3/2) E
d_y Mtheta_z = sqrt(2) R^(3/2) E V
d_y Mztheta  = R (V^2-E^2/2)
d_y Mp       = E^2/2
d_y P0       = 0
d_y pressure = E^2/2, with pressure=P0+Mp.
```

The actual radial velocity derivative rows already include their amplitude derivatives. Products E*V, V^2 and E^2 use the ordinary binomial product rule. Remaining R or R^(3/2) powers are then differentiated once. This recovers all positive radial orders through four while preserving incoming histories at radial order zero.

The main and active end pulse retain nonzero axial contributions. The inactive gap has zero current axial velocity while its incoming linear histories remain nonzero. Post-Rv linear zeros are used only on the admitted terminal/postpulse branch. P0 is independent, and nonzero angular, quadratic and pressure memories are not reset.

Pressure base values remain in Pstar^2 units. Positive radial derivatives use exact pulse Fpulse^2*Pstar^2 units or postpulse Ev0^2 units. Old cap-backed pressure derivative grids are not used as derivative definitions. Directed caps can enclose a source but do not replace its exact positive amplitude.

## Exact native-coordinate derivative scales

All current radial maps are affine in their native coordinate and have Z-independent positive Jacobian J=dlogR/dnative. Thus d_native^K dZ^J equals the corresponding logR/Z row multiplied by the exact source J^K. The adapter appends K*log(J) to the scale graph, retaining the same coefficient object, instead of multiplying coefficients by a Jacobian interval endpoint or bound.

| Charts | Exact J |
|---|---|
| pulse main, exit, gap | 1/mu |
| outer power phase | Lrel-4 |
| steep power phase | Ts |
| waiting phase | true original integral-defined W |
| entrance, gap-end, end, flatten, angular, steep entry/exit, heat | 1 |

Pulse row units are `(R,Fpulse,Pstar,native_Jacobian)`, postpulse row units are `(R,Ev0,Pstar,native_Jacobian)`. Their common Rv conversion remains the accepted predecessor interface. Huge absolute radii, selected frequencies and tiny absolute amplitudes are not materialized.

## Focused evidence

- 60 independent arbitrary-source cumulative/product-rule identities, including nonzero V, signed quadratic energy and independent P0.
- All fifteen actual radius maps independently checked for affine native dependence, the correct Jacobian and Z independence.
- Fifteen actual current routes at fresh Z=.521 compared to recorded source outputs; all 2250 logR/Z and 2250 exact-native rows checked.
- 675 actual velocity-row, 420 first-density, 240 original pulse primitive, 30 independent closed-Gamma log-rate and 300 accepted Rv mixed-seam overlap diagnostics. Function and scale proofs are separate from overlap.
- All 29700 final scale parts belong to the same current expression graph.
- Invalid domains, missing admitted Rv source and wrong Pstar/Ev0 producer metadata are rejected.
- The heat metadata guard was corrected using the actual published pressure-scale field. Thirteen unchanged typed source observations were reused with matched mixed-core bytecodes and refreshed unit/output checks; the two heat routes were then called and checked. Duplicate expensive original pulse quadrature was removed from the check path.
- Accepted constructor and a new heat exterior call at Z=.549, t=9/2 succeed.
- All 1305 defining input byte dependencies match the staged Git index.
- Scoped static source/unit review: **GPT-5.6 Luna / max — PASS**, including the corrected actual heat return field.

Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_mixed_transport.py`, `.json.gz`, `_check.py`, `_check.json`. A standard fresh check works from an accepted owner. The explicit refresh/resume options only reuse same-source typed observations after unchanged mixed-core methods are verified; they do not replace them with receipt values or relax dependencies.

## Completed bounded scopes

- [x] **CURRENT-RP-RAW-HIGH-MIXED-TRANSPORT** All fifteen actual routes now expose exact-scale mixed source rows through total logR/Z order four in [69692731](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/6969273120a8fcbca390bb15d3bf3f0338fe1908). Higher required stress/remainder orders and uniform/global admission remain separate.
- [x] **CURRENT-RP-MIXED-VELOCITY-PRESSURE** Actual current velocity and density-derived pressure derivative interfaces are callable through total order four.
- [x] **CURRENT-RP-MIXED-NATIVE-JACOBIAN** Native derivatives retain exact J^K source factors for all fifteen affine Z-independent maps.

## Detailed next tasks

- [ ] **CURRENT-RP-PULSE-ENTRANCE-MAIN-SEAM** Join t=.02/mu and xi=.02 as defining functions, preserving the lazy boundary and conversion between native derivatives. Use the new ordinary logR grid as the shared comparison coordinate.
- [ ] **CURRENT-RP-PULSE-MAIN-EXIT-SEAM** At xi=10, prove primitive, velocity and pressure total-order-four joins using the actual flat source and incoming memories.
- [ ] **CURRENT-RP-PULSE-EXIT-GAP-SEAM** At xi=11, prove flat support termination while retaining nonzero cumulative memories and their scaled derivatives.
- [ ] **CURRENT-RP-PULSE-GAP-GAPEND-SEAM** Join xi=12 and the formal s=-1/mu boundary, keeping defining geometry separate from outward legal coverage. Do not substitute a rounded reciprocal bound.
- [ ] **CURRENT-RP-PULSE-GAPEND-END-SEAM** At s=-4, join actual backward source functions and compact support derivatives, retaining future energy and pressure memory.
- [ ] **CURRENT-RP-RAW-PULSE-DOMAIN-SEAMS** Aggregate the preceding five source/derivative/coverage joins. Rv is already admitted; all pulse-domain coverage is not.
- [ ] **CURRENT-RP-CLOSED-HEAT-MIXED-SEAMS** Use the same closed Gamma branch to establish waiting/collar and collar/exterior quantitative mixed joins, including pressure and energy.
- [ ] **CURRENT-RP-POSTPULSE-QUANTITATIVE-SEAMS** Connect flatten/power/angular/steep/waiting using explicit source-scale units, derivative errors and positive widths. Reuse existing exact defining identities.
- [ ] **CURRENT-RP-UNIFORM-SELECTED-SEAMS** Extend local derivative interfaces to whole-domain quantitative contracts, including narrow bump supports and flatten. Finite source rows alone do not close this gate.
- [ ] **CURRENT-RP-STRESS-DERIVATIVE-ORDER-LEDGER** Determine the derivative orders consumed by physical stress and flat-remainder recovery. Add genuine higher source rows where needed; do not pad Ur or use Taylor values as missing derivative proofs.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Propagate signed integration, inverse, Taylor/derivative and phase errors through exact scales; define usable point/error contracts separately from finite directed ranges.
- [ ] **CURRENT-FREQUENCY-REPRESENTATION** Keep N=2^(2^J), J=1358356628656378313 lazy, with exact phase operations and declared coordinate limits; never substitute a smaller frequency.
- [ ] **CURRENT-INVERSE-PICARD-ORACLE** Evaluate actual same-source monotone inverse, flat branches, Picard controls and signed quadrature with documented error.
- [ ] **CURRENT-NUMERIC-POINT-ORACLE** Evaluate genuine source values on legal cells under phase/error contracts; a Taylor jet or enclosure is not an unrestricted point.
- [ ] **CURRENT-RP-RAW-PULSE-POINT-HISTORIES** Deliver numeric pulse velocities and five histories with errors while preserving the selected source and P0.
- [ ] **CURRENT-RP-POSTPULSE-FIVE-HISTORIES** Deliver usable numeric postpulse values under true absolute scales; mixed source rows alone do not finish this broader task.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Join core/transition/Rh/O2/O3/compact/quiet to the current fifteen-chart caller under one quantitative derivative contract.
- [ ] **CURRENT-PHYSICAL-TIME-MAPPING** Install original similarity/physical coordinate and time maps, keep huge origin and finite/time offsets separate, include axis limits and apply chain rules to the new derivative rows.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] and pressure from actual current point/error contracts, with cylindrical-to-Cartesian component conversion.
- [ ] **CURRENT-STRESS-DECOMPOSITION** Recover signed divergence-form stress and separate flat remainder from actual physical inertial, pressure and viscous terms.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Prove global margins and decay with max/physical-volume L2 ledgers through supports, exits and annuli.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement actual n=1 recovery and independent moment repair on the common inner domain.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n>=2 order-dependent recovery, independent repairs, curl-preserving truncation and controlled smooth summation. Spatial derivative transport is not temporal coefficient recursion.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Recover both original oscillatory families and mean corrections from admitted stress and recursive coefficients.
- [ ] **CURRENT-AVERAGED-STRESS-CANCELLATION** Establish averaged quadratic momentum-flux cancellation and bound remaining oscillatory/flat terms.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure contraction, relative axial length, velocity/vorticity growth, true material winding, scale recurrence and finite energy from computed fields.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently evaluate corrected forced NS max/volume-L2 targets after the correction layer is installed.

Next agent: implement bounded mixed source seams using the new grid and advance physical/time mapping and point-error contracts. Prioritize delivery and changed-boundary evidence; reuse unchanged accepted proofs. Mark achieved task scopes with their commit and keep the full goal ACTIVE. Predecessor: [Rv common-unit function join](CURRENT_ORIGINAL_RP_COMMON_UNIT_SEAM_2026_10_10.md).
