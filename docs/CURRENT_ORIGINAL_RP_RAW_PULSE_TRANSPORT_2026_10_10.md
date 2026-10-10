# Current six-pulse raw transport and common closed-heat caller (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. Source [b7ac3996](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b7ac39966ab599c96f806fb687e97780af3cea82) extends executable factorized raw histories to all six active pulse charts and connects them to the nine current postpulse charts, using the already closed collar/exterior. The caller exposes similarity velocity components Ur, Utheta and Uz, all five histories, independent P0, full pressure and first native history derivatives. Unrestricted numerical physical u/v/w, uniform/global mixed admission, cone and true temporal recursion remain open.

## Concrete implementation

`CurrentOriginalRpRawPulseTransport(accepted_closed_heat)` requires the accepted existing-repair heat closure. It keeps the same selected pulse, C5/energy/angle providers, exact repair/future branch, beta source and independent P0. A cheap independent expression graph holds all published scales. No source selection, repair or old incoming construction is replayed when an accepted owner is supplied.

All 15 actual routes are callable: entrance, main, exit, gap, gap-end, end, flatten, outer power/angular, steep entry/power/exit, waiting, collar and Gamma exterior. Heat calls retain the closed source coefficients and rebase only their exact scale functions into the common caller graph. Closed coefficient objects are unchanged; forward constants remain diagnostic. Every final history/velocity/derivative scale belongs to `owner.expression_owner(chart).graph`, which is the same current raw expression graph on all 15 charts.

## Pulse units and exact cancellations

Define the exact Z-independent source

```
Y = log(R/Rp)
Fpulse = exp(-(.5+mu)*Y)
Utheta = Pstar * Fpulse * u(Z)
```

The inlet u(Z) is the same current C5 source function. Pulse history packets are normalized by local Utheta; base inlet packets use Pstar units. The new adapter distinguishes these schemas.

| Output | Exact factor | Directed axial coefficient |
|---|---|---|
| Mz | R Fpulse Pstar | u*m |
| Mtheta | R^(3/2) Fpulse Pstar | sqrt(2)*u*X |
| Mtheta_z | R^(3/2) Fpulse^2 Pstar^2 | sqrt(2)*u^2*H |
| Mztheta | R Fpulse^2 Pstar^2 | u^2*e |
| Mp, P0, pressure | Pstar^2 | original packet pressure rows |
| Utheta | Fpulse Pstar | u |
| Uz | Fpulse Pstar | u*B |
| Ur | sqrt(R) Fpulse Pstar | u*A/sqrt(2) |

Here B=Uz/Utheta, A=Ur/(sqrt(R/2)*Utheta), and m/X/H/e are actual transported pulse primitives. Histories, Utheta and Uz retain true axial order 5; Ur retains order 4 after radial recovery. No C5 radial velocity is invented.

For R^r Fpulse^e Pstar^p, source logs are reduced before arithmetic:

```
entrance: r*logRp + p*logPstar + (r-e/2-e*mu)*t
xi charts: r*logRp + p*logPstar + (r-e/2)*xi/mu - e*xi
end charts: r*logRp + p*logPstar + 13*(r-e/2)/mu
            -13*e + (r-e/2-e*mu)*s
```

Thus the energy R*Fpulse^2 and radial velocity sqrt(R)*Fpulse have exact inverse-mu pulse cancellation. A finite offset is never recovered by subtracting enormous rounded logarithms. Caps used by native packets remain directed coefficient bounds, never the defining Fpulse or an unrestricted point value.

Pulse factors use explicit `(R,Fpulse,Pstar)` units. Postpulse/heat factors use `(R,Ev0,Pstar)`, where `Ev0=Pstar*U0*exp(-13/(2mu)-13)`. Consumers must read `scale_units`; the different pulse and postpulse coefficient normalizations cannot be compared directly. At Rv, convert Fpulse*Pstar to Ev0/U0 using the actual U0 defining function.

## Nonzero axial source and all five density equations

The main and active end sources retain nonzero axial velocity. The inactive gap has Uz=0 while its cumulative histories retain incoming/end-tail memory. Only the actual Rv endpoint has source-proved zero Mz/Mtheta_z; no terminal zero theorem is applied inside active pulse charts.

With J=dlogR/d(native coordinate):

```
d_native Mz       = J R Uz
d_native Mtheta   = J sqrt(2) R^(3/2) Utheta
d_native Mtheta_z = J sqrt(2) R^(3/2) Utheta Uz
d_native Mztheta  = J R (Uz^2 - Utheta^2/2)
d_native Mp       = J Utheta^2/2
d_native P0       = 0
d_native pressure = d_native Mp
```

J is 1/mu on main/exit/gap xi charts, and 1 on entrance/gap-end/end charts. It is applied once; the existing native mixed packets already contain ordinary logR derivatives. Pressure history units are Pstar^2, while its derivative uses Fpulse^2 Pstar^2. Full pressure remains P0+Mp, with no axis reset or loss of incoming memory.

## Usage and evidence

```python
owner = CurrentOriginalRpRawPulseTransport(accepted_closed_heat)
main = owner.evaluate('pulse_main', '.491', '5/4')
main['similarity_velocity']['Uz'].report()
main['raw_histories']['Mztheta'].report()
main['first_native_radial_derivatives']['Mtheta_z'].report()
heat = owner.evaluate('heat_exterior', '.479', 4)
heat['raw_histories']['pressure'].report()
graph = owner.expression_owner('heat_exterior').graph
```

- Fifteen actual current calls at Z=.479 cover all routes, including active main/end and inactive gap memory.
- 60 independent absolute scale identities and 12 exact energy/radial inverse-mu cancellations.
- Independent general cumulative integrals with nonzero Uz establish all five signed density equations, the single Jacobian and independent P0.
- 642 history/density coefficient checks and 90 overlap diagnostics against original physical velocity mixed rows. Overlap diagnostics do not certify uniform high mixed interfaces.
- 1572 final scale parts belong to one common current expression graph.
- Main/end axial sources remain present; gap memory is not zeroed; actual Rv linear terminal zeros are preserved.
- Invalid chart domains and an unclosed source bypass are rejected.
- Accepted constructor and a fresh main call at Z=.491, xi=5/4 succeed. Radial velocity retains C4 and axial velocity C5.
- All 1297 dependency bytes match the staged Git index.
- Scoped static review: **GPT-5.6 Luna / max — PASS**, including the closed coefficient rebase into the common expression graph.

Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_raw_pulse_transport.py`, `.json.gz`, `_check.py`, `_check.json`.

## Completed at source/enclosure scope

- [x] **CURRENT-RP-RAW-PULSE-SCALE-FUNCTIONS** Install the six actual pulse radius/decay scale functions with symbolic cancellation.
- [x] **CURRENT-RP-RAW-PULSE-FIVE-HISTORIES** Restore all five active pulse histories, independent P0 and full pressure without deleting nonzero terms.
- [x] **CURRENT-RP-RAW-PULSE-NATIVE-DERIVATIVES** Restore the general five density equations with one true native Jacobian.
- [x] **CURRENT-RP-CLOSED-HEAT-CALLER-INTEGRATION** Supply one 15-chart current source caller selecting closed collar/exterior and retaining the same repair/providers.
- [x] **CURRENT-RP-RAW-PULSE-HISTORIES** The factorized source portion of the older aggregate task is done; numerical raw point delivery remains open below.

## Next bounded implementation work

- [ ] **CURRENT-RP-COMMON-UNIT-SEAMS** Convert pulse `(R,Fpulse,Pstar)` and postpulse `(R,Ev0,Pstar)` using the exact U0/Rv identity before comparing derivative/history rows. Prove defining-function continuity at pulse-end/flatten; retain diagnostic enclosures without promoting overlap to a high mixed proof.
- [ ] **CURRENT-RP-RAW-PULSE-DOMAIN-SEAMS** Keep original six guards and supplemental reciprocal gap coverage. Complete source/physical derivative seams for entrance/main/exit/gap/gap-end/end; do not replace s=-1/mu with an interval endpoint as a defining boundary.
- [ ] **CURRENT-RP-RAW-HIGH-MIXED-TRANSPORT** Factor actual native y/Z derivative rows through required total orders for pulse and postpulse. Differentiate every radial amplitude and coordinate scale. Ur's C4 recovery needs its own error/remainder contract.
- [ ] **CURRENT-RP-CLOSED-HEAT-MIXED-SEAMS** Use the same closed Gamma branch to prove waiting/collar and collar/exterior derivative joins, including pressure and energy. Source Cp/Dtheta identities alone do not close quantitative mixed joins.
- [ ] **CURRENT-RP-UNIFORM-SELECTED-SEAMS** Finish whole-domain mixed contracts over all 15 routes, including narrow end supports and flatten. Current factored cell bounds do not finish uniform admission.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Carry signed quadrature/inverse/derivative errors through each exact positive source scale and cancellation. Specify usable point-value contracts without midpoint or cap choices.
- [ ] **CURRENT-RP-RAW-PULSE-POINT-HISTORIES** Deliver actual numeric pulse histories/velocities under a documented admissible coordinate/phase and error contract, beyond directed factorized enclosures.
- [ ] **CURRENT-RP-POSTPULSE-FIVE-HISTORIES** Deliver actual numeric raw postpulse point values with error; keep this broader task open after source transport.
- [ ] **CURRENT-FREQUENCY-REPRESENTATION** Preserve selected N=2^(2^J), J=1358356628656378313 in lazy phase operations; do not substitute a smaller frequency.
- [ ] **CURRENT-INVERSE-PICARD-ORACLE** Evaluate same-source monotone inverse, flat branches, Picard controls and signed integrals with bounded error.
- [ ] **CURRENT-NUMERIC-POINT-ORACLE** Evaluate true defining functions over legal cells with phase limits/error, separately from range covers, source handles and Taylor jets.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Join core/transition/Rh/O2/O3/compact/quiet to this current pulse/exterior caller under one admitted derivative contract.
- [ ] **CURRENT-PHYSICAL-TIME-MAPPING** Map the current lazy source geometry and vector components into original physical coordinates/time, keeping enormous origin and finite/time offsets separate and retaining axis limits.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] and pressure with current coefficients and point/error/derivative contracts.
- [ ] **CURRENT-STRESS-DECOMPOSITION** Recover full signed divergence-form stress and separate flat remainder from current physical inertial, pressure and viscous terms.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Prove global margins/flat decay and report maximum/volume norms separately through exits, supports, matching and heat collar.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement actual n=1 coefficient recovery and independent moment repair on the common inner domain.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement n>=2 recovery, order-dependent repair, curl-preserving truncation and controlled smooth summation.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Build both original oscillatory families and mean corrections from admitted stress/recursion coefficients.
- [ ] **CURRENT-AVERAGED-STRESS-CANCELLATION** Establish averaged quadratic momentum flux cancellation and bound residual oscillatory/flat terms.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure core contraction, relative axial length, velocity/vorticity growth, true material winding, scale recurrence and finite energy from computed fields.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently evaluate corrected forced NS L-infinity/volume-L2 targets after the correction layer is installed.

Agents: implement the next unmet dependency, reuse unchanged accepted proofs, publish only changed-boundary evidence, and mark achieved scopes with the source commit. The full goal remains ACTIVE. Predecessor: [existing-repair closed heat](CURRENT_ORIGINAL_RP_SAME_REPAIR_HEAT_CLOSURE_2026_10_10.md).
