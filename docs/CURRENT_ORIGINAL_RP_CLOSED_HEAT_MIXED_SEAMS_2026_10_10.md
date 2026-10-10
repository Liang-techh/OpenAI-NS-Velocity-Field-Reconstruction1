# Update: remaining six raw postpulse function joins installed (2026-10-10)

[Source d00208c7](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d00208c721b24f3765012803bd30f16302741dc9) completes the four bounded task IDs below covering six interfaces, and aggregates all eight postpulse function joins. The broader quantitative/uniform task remains open. Read the [new handoff](CURRENT_ORIGINAL_RP_POSTPULSE_MIXED_SEAMS_2026_10_10.md); next priority is physical/time and point/error delivery.

---

# Current original raw closed-heat seams (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. [Source 274d5518](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/274d5518b2844375aaae24b035c7ed832649f6d0) installs the two same-current raw function joins: waiting phase=1 to closed collar t=0, and closed collar t=3 to full Gamma exterior t=3. Both sides retain the original analytic P0, full energy/pressure histories, unique repair, complete infinite Gamma future and exact absolute radius/amplitude functions. Ordinary logR/Z derivatives have total order at most four. This is a source-function interface with directed consistency bounds; numerical physical u/v/w, uniform/global norms, stress/remainder, genuine temporal coefficient recursion and oscillatory cancellation remain open.

## Callable interface

```python
owner = CurrentOriginalRpClosedHeatMixedSeams(accepted_pulse_mixed_seam_owner)
view = owner.evaluate('waiting_collar', '.521')
view['left']['log_radius_mixed_rows']['Mztheta']['y2_Z1'].report()
view['right']['native_coordinate_mixed_rows']['pressure']['n1_Z3'].report()
exterior_join = owner.evaluate('collar_exterior', '.521')
```

| Seam | Left actual source | Right actual source | Left / right native Jacobian |
|---|---|---|---|
| waiting_collar | waiting phase=1 | closed collar t=0 | true waiting W / 1 |
| collar_exterior | closed collar t=3 | full Gamma exterior t=3 | 1 / 1 |

Each side exposes five complete cumulative histories Mz/Mtheta/Mtheta_z/Mztheta/Mp, independent P0, absolute pressure and three velocity components. Every output has fifteen ordinary mixed derivative indices k+j<=4. Source units are (R,Ev0,Pstar); positive radial pressure derivatives retain Ev0 squared, while base pressure/P0 retain Pstar squared. Native rows multiply the exact lazy positive J^k scale; phase-to-y conversion is W^(-k). No rounded waiting-root endpoint becomes the defining Jacobian or radius.

## Function identities consumed

The original source-bound collar/Gamma theorem is recomputed and required to equal both its accepted receipt and the live current heat provider's theorem. It identifies theta, angular X, energy and pressure at each seam as whole-Z functions, including their ordinary mixed derivatives through total order four and base axial derivatives through order five. The current source proof separately binds the shared exact infinite Gamma integrands, both epsilon energy atoms and their normalization:

```
H_waiting(Z) = E_collar(0,Z)/(1-epsilon)^2
energy_waiting_exit = H_waiting/2
theta_waiting_exit = theta_base*(1-epsilon)
```

The same-current unique angular repair identifies Dtheta as the zero function. Original analytic pressure identification proves P0 plus the complete native pressure integral is zero, so Cp is the zero function. The closed collar/exterior pressure is the original forward pressure function after this identity; its diagnostic forward Cp enclosure is retained separately and is never added back as a pressure patch. At t=3 the flat cutoff joins the original full Gamma source, retaining the complete tails rather than a finite S polynomial.

AST bindings target the actual closed provider's assignments and raw outputs. Raw transfer uses Mtheta=sqrt(2)*R^(3/2)*Ev0*theta*X, Mztheta=R*Ev0^2*theta^2*energy and Mp=P-P0. The inherited admitted selected terminal identities propagate zero meridional histories. Positive radial orders use the existing general cumulative density equations, preserving their exact radius and true amplitude scales.

## Focused evidence

- 120 source-bound primitive mixed identities on the two current joins; canonical base axial-five identities retained.
- Two exact lazy absolute-radius equalities, 300 exact common-scale derivative identities and 600 exact native Jacobian-power identities.
- 300 actual same-source signed mixed derivative difference/overlap diagnostics. Function equality is established separately before these comparisons.
- Unknown seams, missing full-Gamma function certificate, changed complete-energy normalization and wrong waiting boundary are rejected.
- Accepted constructor succeeds on the same live upstream graph. The actual source observations remain typed source outputs; no derivative scalar is reconstructed from an archived receipt.
- All 1314 defining dependency byte hashes match the Git index. The existing read-only **GPT-5.6 Luna / max** worker's scoped review found no material blocker.

Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_closed_heat_mixed_seams.py`, `.json.gz`, `_check.py`, `_check.json`. No upstream selection, repair, inverse-radius or provider solve is replayed for this milestone. Previously accepted defining source files are unchanged.

## Completed bounded tasks

- [x] **CURRENT-RP-WAITING-CLOSED-COLLAR-SEAM** Actual current waiting phase=1 / closed collar t=0, original pressure/energy memories, raw mixed4 rows and exact native W distinction.
- [x] **CURRENT-RP-CLOSED-COLLAR-EXTERIOR-SEAM** Actual current t=3 branches, complete Gamma tails, exact common units and raw mixed4 histories/velocity/pressure.
- [x] **CURRENT-RP-CLOSED-HEAT-MIXED-SEAMS** Aggregate those two source-function joins with directed signed diagnostics. Uniform and global quantitative norms remain a separate task.

## Next tasks and acceptance

- [x] **CURRENT-RP-FLATTEN-POWER-SEAM** Join flatten t=100 and outer-power phase=0 on the accepted mixed caller. Instantiate the existing full-future flatten theorem, preserve nonzero angular/quadratic histories, prove exact radii/Ev0 scales and raw mixed4 transfer. Publish actual signed differences.
- [x] **CURRENT-RP-POWER-ANGULAR-SEAM** Join power phase=1 and angular t=-4. Use the current repaired angular source and full future. Require exact local radius and amplitude equality before comparing fifteen derivative rows per output.
- [x] **CURRENT-RP-ANGULAR-STEEP-ENTRY-SEAM** Join angular t=0 and steep-entry phase=0. Consume the original arbitrary-function theorem on the current repair and retained analytic P0; keep actual history memory intact.
- [x] **CURRENT-RP-STEEP-PHASE-SEAMS** Join entry phase=1 / steep-power phase=0, steep-power phase=1 / exit phase=0, and exit phase=1 / waiting phase=0. Keep exact Ts/W Jacobian powers; compare ordinary y rows rather than unequal native coordinates.
- [ ] **CURRENT-RP-POSTPULSE-QUANTITATIVE-SEAMS** Aggregate the six remaining postpulse joins with these two heat joins. Establish function identities and exact unit maps plus signed bounds. State separately which bounds are finite observations and which are uniform over full domains.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Connect current core/transition/Rh/O2/O3/compact/quiet/Rp owners to the admitted mixed caller. Preserve one original defining source family; supply quantitative norm contracts across the additional joins.
- [ ] **CURRENT-PHYSICAL-TIME-MAPPING** Implement the original similarity-to-physical coordinate/time functions on this current original-radius graph. Preserve the huge origin/finite offset separation, true N, moving radial/axial chain rules, fixed-x time derivative and axis limits. Convert published Ur coefficient by sqrt(2) where the generic recovery convention requires it.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Propagate real quadrature/inverse/phase errors through exact scales and provide a defined point/error pair. A directed source jet, source range or upper bound alone does not complete this task.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] from the admitted physical map and point/error interface. Include axis regularity and independent divergence checks.
- [ ] **CURRENT-RP-PULSE-UNIFORM-BOUNDS** Bound actual signed derivative norms throughout all pulse charts and narrow beta supports. Endpoint joins and finite observations do not provide these norms.
- [ ] **CURRENT-RP-STRESS-DERIVATIVE-ORDER-LEDGER** List the derivative orders required by stress/divergence/remainder; compute genuine missing orders and retain Ur's real C4 limit without padding jets.
- [ ] **CURRENT-STRESS-CONE-AND-REMAINDER** Recover signed divergence-form stress, global cone margins, the separate flat remainder and physical-volume max/L2 bounds on the actual field.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** Implement the true n=1 recovery equations and independent moment repair on the common inner domain; coordinate rescaling alone does not satisfy this task.
- [ ] **CURRENT-TEMPORAL-RECURSION-HIGHER** Implement actual n>=2 equations, independent repairs, curl-preserving truncation, finite-order remainder and a controlled smooth sum.
- [ ] **CURRENT-OSCILLATORY-FAMILIES** Recover both oscillatory families and mean corrections from admitted stress/recursive coefficients; demonstrate averaged quadratic stress cancellation.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** Measure contraction, axial relative length, velocity/vorticity growth, material winding, true scale recurrence and finite energy from computed physical fields.
- [ ] **CURRENT-FULL-CARTESIAN-RESIDUAL** Independently evaluate the corrected forced Navier–Stokes max and physical-volume L2 residuals after the correction layer exists.

Next agent: implement the six remaining postpulse raw joins using the existing accepted mixed runtime and generic source theorems, then the original physical/time and point/error interface. Mark bounded completed tasks with their source commit; keep full reconstruction ACTIVE. Predecessor: [five pulse mixed seams](CURRENT_ORIGINAL_RP_PULSE_MIXED_SEAMS_2026_10_10.md).
