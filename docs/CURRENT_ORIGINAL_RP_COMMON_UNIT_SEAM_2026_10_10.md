# Current Rv common-unit function join and mixed logR/Z rows (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [e445ee9b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e445ee9b23f5195fe9ef26a8ff729472663ffe3c) joins the actual current pulse endpoint `pulse_end(s=0)` to `flatten(t=0)` in one explicit `(R,Ev0,Pstar)` scale schema. The same-source function join retains all five histories, independent P0, pressure, and the original future/repair. It also makes endpoint mixed derivatives through total order four callable. Uniform/global admission, unrestricted numerical physical u/v/w, stress cone, temporal coefficient recursion and oscillatory corrections remain open.

## Implemented interface

```python
owner = CurrentOriginalRpCommonUnitSeam(accepted_pulse_closed_caller)
view = owner.evaluate('.537')
view['reduced_pulse_terminal']['Mtheta'].report()
view['actual_flatten_terminal']['pressure'].report()
view['pulse_log_radius_mixed_rows']['pressure']['y2_Z1'].report()
view['flatten_log_radius_mixed_rows']['Utheta']['y1_Z3'].report()
```

The constructor reuses the accepted fifteen-chart caller and existing selected pulse, flatten, future, repair, context and datum. It does not construct a second repair or replay incoming solves. The future-energy facade is bound to the actual C5 provider; facade identity is not confused with the underlying future owner.

Every final exact scale expression belongs to the current caller's graph. Values are directed source enclosures times lazy exact positive scales, never an enormous materialized R or a selected amplitude midpoint. The three endpoint views deliberately retain separate meanings:

- `converted_pulse_terminal_enclosures`: the literal directed division of the actual pulse rows by the positive U0 enclosure.
- `reduced_pulse_terminal`: the exact admitted u/U0=q^-1 function reduction, applied before interval arithmetic.
- `actual_flatten_terminal`: the independently called original flatten provider at t=0.

## Exact common-unit conversion

```
Rv = Rp * exp(13/mu)
Fv = exp(-13/(2mu)-13)
Ev0 = Pstar * U0 * Fv
u(Z) = U0/(1+Z^2)

pulse units (R,Fpulse,Pstar):(r,e,p)
   -> post units (R,Ev0,Pstar):(r,e,p-e)
      coefficient -> coefficient/U0^e, only at Rv
```

All ten published output power schemas obey the exact log-scale identity. This conversion is rejected at an active end-pulse point s=-3. Both endpoint coordinate Jacobians are one; the main pulse's 1/mu Jacobian is not transferred here. The formal gap boundary s=-1/mu remains lazy and is outside this adapter's scope.

The exact u=U0/q identity is independently checked against the accepted complete native-frame source formula. U0 is the Z-independent defining function used by the current raw graph. Its directed box encloses the function and is not an equality proof or a numerical choice.

## Function join and pressure memory

The current selected endpoint binding checks unchanged native assignments and same-object providers, and rebinds the accepted original endpoint theorem: 155 functional, 34 inlet/datum and 24 endpoint identities. The new constructor explicitly verifies that theorem's family, byte ledger and scope. Directed overlap remains diagnostic.

At Rv, theta=1/q, X=Xv, and energy is the same complete future half. Linear histories, Uz and Ur are zero by empty future support at this actual endpoint. Nonzero angular, quadratic and pressure memories remain present. Independent pressure is

```
P(Rv,Z) = P0(Z) + Pin/q^2
          + (U0/q)^2 * (1-exp(-13/mu-26))/(2*(1+2mu)).
```

The pressure decay cap is an enclosure, not the defining exponential. The adapter does not reset P0, discard terminal energy, infer zeros from a box containing zero, or replace flatten's Rp reference data with Rv histories.

## Callable mixed derivatives

Both sides expose ten outputs times fifteen ordinary derivative indices `yK_ZJ`, K+J<=4. y=log similarity R; these are not cylindrical-r, physical-z, time or Cartesian derivatives. Axial factorials are applied once. Base Ur retains order four and other current source histories/velocities retain their admitted order five.

Let bp=1/2+mu and theta=1/q. For K>=1, the normalized radial derivative rows are:

| Output | Exact source factor | Directed axial coefficient |
|---|---|---|
| Utheta | Ev0 | theta*(-bp)^K |
| Mtheta | R^(3/2)*Ev0 | sqrt(2)*theta*(1-mu)^(K-1) |
| Mztheta | R*Ev0^2 | -theta^2*(-2mu)^(K-1)/2 |
| Mp, pressure | Ev0^2 | theta^2*(-(1+2mu))^(K-1)/2 |
| Mz, Mtheta_z, Uz, Ur | their explicit common units | source zero at Rv |
| P0 | Pstar^2 | zero radial derivative |

The R powers contribute to these derivatives. The pressure base remains in Pstar^2 units while radial derivatives use exact Ev0^2 units; no tiny cap is substituted. Independent exact ODE germs verify the product rules and axial derivatives. The original providers' actual mixed velocity rows are also retained as comparison evidence.

## Focused evidence

- 10 independent exact common-unit identities; exact Rv geometry and both unit Jacobians.
- 40 independent radial/axial source product-rule identities plus analytic inlet and absolute pressure identities.
- Fresh Z=.521 and a directed cell Z in [.520,.522] invoke both actual providers.
- 236 endpoint coefficient, 300 paired mixed-row and 180 original-provider velocity-row overlap diagnostics; none is used as the defining function proof.
- All 3960 final scale parts belong to one current source graph.
- Unclosed bypass, illegal Z domains and conversion at an active pulse point are rejected.
- Accepted constructor and a fresh Z=.537 call succeed with ten outputs and fifteen derivative indices each.
- All 1301 input byte dependencies match the staged Git index.
- Scoped static review: **GPT-5.6 Luna / max — PASS** for the local common-unit and derivative source adapter.

Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_common_unit_seam.py`, `.json.gz`, `_check.py`, `_check.json`.

## Task handoff

- [x] **CURRENT-RP-COMMON-UNIT-SEAMS** Done for the only pulse-to-post amplitude-unit conversion at Rv. Exact conversion and same-current function join are installed in [e445ee9b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e445ee9b23f5195fe9ef26a8ff729472663ffe3c). This completion does not include other pulse-domain seams or a global mixed certificate.
- [x] **CURRENT-RP-RV-MIXED-LOGR-Z4** Rv endpoint rows through total order four are callable in explicit source units in [e445ee9b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e445ee9b23f5195fe9ef26a8ff729472663ffe3c).
- [ ] **CURRENT-RP-PULSE-ENTRANCE-MAIN-SEAM** Use the exact t=.02/mu / xi=.02 source boundary, retaining lazy coordinate and correct Jacobian. Join full cumulative histories and velocity derivative rows.
- [ ] **CURRENT-RP-PULSE-MAIN-EXIT-SEAM** At xi=10, use the actual flat original main source; prove same histories and total-order-four logR/Z rows.
- [ ] **CURRENT-RP-PULSE-EXIT-GAP-SEAM** At xi=11, prove B support termination, retain nonzero history memory, and transport actual pressure/velocity rows.
- [ ] **CURRENT-RP-PULSE-GAP-GAPEND-SEAM** Keep the exact formal s=-1/mu / xi=12 identity and separate legal outward coverage. Never choose an interval endpoint as the defining seam.
- [ ] **CURRENT-RP-PULSE-GAPEND-END-SEAM** At s=-4, join actual backward histories, compact bump support and derivative rows without deleting future memory.
- [ ] **CURRENT-RP-RAW-HIGH-MIXED-TRANSPORT** Extend exact-scale mixed history/velocity transport throughout all fifteen routes using original actual source rows, radial product rules and Jacobians; preserve Ur C4 recovery/error limits.
- [ ] **CURRENT-RP-RAW-PULSE-DOMAIN-SEAMS** Aggregate the five pulse boundary tasks only after their source, derivative and coverage evidence is installed. Local Rv completion does not close this task.
- [ ] **CURRENT-RP-CLOSED-HEAT-MIXED-SEAMS** Use the same closed Gamma branch for waiting/collar and collar/exterior derivative joins, including pressure and energy in correct units.
- [ ] **CURRENT-RP-POSTPULSE-QUANTITATIVE-SEAMS** Connect flatten/power/angular/steep/waiting with actual error and derivative contracts. Reuse accepted defining identities and avoid broad unchanged validation.
- [ ] **CURRENT-RP-UNIFORM-SELECTED-SEAMS** Establish quantitative whole-domain contracts for all fifteen routes and narrow supports. Separate bounded enclosures, function identity and numerical usability.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Propagate signed integration, inverse and derivative errors through exact scales. Specify error contracts before extracting point values.
- [ ] **CURRENT-FREQUENCY-REPRESENTATION** Preserve N=2^(2^J), J=1358356628656378313 using lazy phase operations and explicit admissible-coordinate limits; do not reduce N.
- [ ] **CURRENT-INVERSE-PICARD-ORACLE** Implement same-source monotone inverse, flat branches, Picard controls and signed quadrature with bounded errors.
- [ ] **CURRENT-NUMERIC-POINT-ORACLE** Evaluate genuine defining functions on legal cells under a phase/error contract. Source enclosures and Taylor jets alone are not unrestricted points.
- [ ] **CURRENT-RP-RAW-PULSE-POINT-HISTORIES** Deliver numerical pulse velocity/history points with errors while preserving the selected phase, all five cumulative memories and independent P0.
- [ ] **CURRENT-RP-POSTPULSE-FIVE-HISTORIES** Deliver usable numerical postpulse values/errors under exact absolute scales; factorized source completion alone does not close this broader task.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Connect core/transition/Rh/O2/O3/compact/quiet to this current pulse/exterior owner under one derivative contract.
- [ ] **CURRENT-PHYSICAL-TIME-MAPPING** Preserve huge origin, finite radial offset and time offset separately; implement original physical coordinate/time mapping and axis limits.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Supply callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] and pressure with actual current-source point/error/derivative contracts.
- [ ] **CURRENT-STRESS-DECOMPOSITION** Recover signed divergence-form stress and a separate flat remainder from actual physical inertial, pressure and viscous terms.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Establish global signs/margins and decay, with max and physical-volume L2 ledgers through all relevant supports and annuli.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement the actual n=1 recovery and independent moment repair on the common inner domain.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n>=2 recovery, n-dependent repair, curl-preserving truncation and controlled smooth summation. Radial derivative recurrences are not temporal coefficient recursion.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Construct both original oscillatory families and mean corrections from admitted stress and recursion coefficients.
- [ ] **CURRENT-AVERAGED-STRESS-CANCELLATION** Verify averaged quadratic flux cancellation and remaining oscillatory/flat errors.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure contraction, relative axial elongation, velocity/vorticity growth, true material winding, recurrence and finite energy from computed fields.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently evaluate the corrected forced NS max/volume-L2 targets after the correction layer is installed.

Next agent: implement **CURRENT-RP-RAW-HIGH-MIXED-TRANSPORT** and the associated bounded pulse seams, reusing accepted source providers. Publish changed-boundary evidence and the commit, mark only achieved scopes, and retain full goal ACTIVE. Predecessor: [fifteen-chart pulse/closed exterior](CURRENT_ORIGINAL_RP_RAW_PULSE_TRANSPORT_2026_10_10.md).
