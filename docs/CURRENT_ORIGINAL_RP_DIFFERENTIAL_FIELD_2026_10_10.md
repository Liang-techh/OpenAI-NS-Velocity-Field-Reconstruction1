# Original Cartesian derivatives and correlated linear fields — 2026-10-10

Later progress: [CURRENT_ORIGINAL_RP_AXIAL_VORTICITY_2026_10_10.md](CURRENT_ORIGINAL_RP_AXIAL_VORTICITY_2026_10_10.md) restores negative nonzero omega_z and its passing 0.1% total width using the actual original theta radial law. The earlier unresolved direct-sum result remains historical; broader time/global/stress/recursion tasks remain open.

The supported original pulse_exit field now exposes 140 spatial derivatives through total order four and four fixed-physical-position time derivatives as live nonzero source-scale handles. It also computes 13 signed linear differential fields from the actual shared source primitives: divergence, three vorticity components, three velocity Laplacians and six symmetric strain components. Full goal **ACTIVE / INCOMPLETE**.

Source commit: [ff569e26](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ff569e26bacace902e5ec8366f75ec2b7744df9e). New source quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_differential_field.py`, `.json.gz`, `_check.py`, `_check.json`.

## Actual supported point and results

The physical delivery remains pulse_exit, xi=21/2, Z=371/1000, finite log-time offset=1/5, theta=7/10. All 144 derivative handles are nonzero at this delivery. They retain the existing exact physical scales, source pressure datum, signed coefficient enclosures and ordinary derivative units. No second factorial, native Jacobian or amplitude derivative is introduced.

| Result | Supported result at this delivery |
| --- | --- |
| Spatial derivatives | All four components ux/uy/uz/p, total orders 0 through 4: 140 rows |
| Time derivatives | One fixed-physical-position derivative of each component: 4 rows |
| Vorticity | omega_x negative; omega_y positive; omega_z sign unresolved |
| Velocity Laplacians | All three negative |
| Symmetric strain | 00 positive, 01 negative, 02 positive, 11 positive, 12 positive, 22 negative |
| Divergence | Numeric enclosure contains zero; original analytic source divergence identity retained separately |
| Exact primitive cancellation | Two Ur groups cancel in omega_z before interval evaluation; two source groups cancel in divergence |
| Accuracy counts | 59 factored target passes; 26 total physical relative-width passes; 96 nonconstant residual-logC rows |
| Ordinary absolute numeric export | Zero nonzero materializations; unrestricted xyz/t API remains false |

The 0.1% total physical width flags are inherited without promotion. Base u/v widths remain unresolved and base w/p widths pass. Exact source-scale handles make derivative arithmetic executable; they do not establish global accuracy, vorticity growth with time, stress admissibility, flatness or scale recursion.

## Source correlation and admission

Each linear form first combines exact rational weights for identical derivative rows. It then groups actual original source primitives by source label, source derivative, source units, exact scale powers and canonical log-scale polynomial. An exact zero graph operator is cancelled before its source coefficient is enclosed. Surviving coefficients are summed only in matching typed scale groups, then rebased against an existing exact source scale with the inherited directed exponential/tail guard. A shared scale ratio is applied once to the complete coefficient group.

Independent interval overlap never proves an operator identity. In particular, the numeric divergence sum is not promoted to exact zero. The accepted original recovery algebra supplies a separate analytic divergence witness. omega_z still contains different Utheta derivative primitives; its unresolved small difference needs an actual source derivative law, rather than numerical cancellation of two almost equal intervals.

Fields and derivative handles belong to their live issuing owners. The adapter requires the exact 4-component by 35-spatial-plus-1-time grid and ordinary order-zero source rows; native-Jacobian and Taylor-polynomial rows are rejected. Defining DAG nodes, source context/coefficient/order, source units/powers, bindings, aliases and whitelist are sealed. The accepted receipt binds current code/candidate hashes, the Cartesian map definitions and all 13 operator definitions.

The focused checker independently derives curl/strain/divergence/Laplacian orientation with symbolic Cartesian functions, verifies all 144 actual typed rows, independently checks exact cancelled operators, and encloses all 13 shared-scale sums at 800 digits. It rejects copied fields, missing time rows, unavailable components/orders, mixed units, removed handles, changed aliases, substituted primitives and Taylor-order mutation. Accepted construction, evaluation, report and derivative retrieval pass. Only this new adapter's checks and source-index audit are run; inherited accepted suites are reused.

The existing read-only worker was reused: **GPT-5.6 Luna / max**. No new worker or Astra child was spawned.

## Executable use

Use the complete typed original bootstrap in [the pressure checkpoint](CURRENT_ORIGINAL_RP_REMAINING_PRESSURE_TAIL_2026_10_10.md#fresh-process-use) and construct the accepted SourceProduct and ScaleValues owners as in [the value checkpoint](CURRENT_ORIGINAL_RP_SCALE_VALUE_ARITHMETIC_2026_10_10.md#executable-use). Then append:

```python
from lei_ren_part1_paper_compliant_current_original_Rp_differential_field import CurrentOriginalRpDifferentialField as DifferentialField

differential = DifferentialField(arithmetic)
field = differential.evaluate(delivery, '1/1000')
report = differential.report(field)
ux_x = differential.derivative(field, 'ux', (1, 0, 0))
ux_t = differential.derivative(field, 'ux', ('t',))
ux_x_record = arithmetic.export(ux_x)
omega_z = report['actual_linear_differential_fields']['vorticity_z']
laplacian_u = report['actual_linear_differential_fields']['Laplacian_ux']
```

Keep fields and handles with their owners. JSON exports are reports, not hydrated live fields. Prepare independent source queries before constructing frozen derived views. Do not reset failed cache guards or substitute independently enclosed row products into exact nonlinear source identities.

## Completed bounded tasks

- [x] **CURRENT-RP-SCALE-VALUE-ALL-ROWS-POINT** Issue live source-scale handles for the complete original spatial4/time1 grid at the supported point; preserve component, derivative, units, sign and absolute-width flags.
- [x] **CURRENT-RP-SOURCE-CORRELATED-LINEAR-POINT** Combine exact same-source primitive operators before interval evaluation and preserve shared scale correlation for actual divergence, curl, strain and Laplacian.
- [x] **CURRENT-RP-DIFFERENTIAL-POINT-API** Supply accepted derivative retrieval and 13 linear reports with exact-grid, ordinary-row, ownership, source-map and receipt guards.

## Next executable tasks

- [ ] **CURRENT-RP-EXACT-AXIAL-VORTICITY-LAW** Read the actual original pulse-exit Utheta source definition, including incoming/pulse contributions. Prove the radial coefficient dependence and derive Utheta_y + Utheta/2 as one source expression. If the pure power law holds, retain its exact mu factor and physical sqrt(2/R), lambda powers before bounding. Bind the source AST and verify the resulting expression against the actual Cartesian operator. Do not replace an interval difference by zero or assume a radial law from equal coefficient intervals. Acceptance: actual source-issued omega_z signed result with unchanged original parameters; explicit chart/domain and no growth/recursion claim.
- [ ] **CURRENT-RP-SIGNED-BACKGROUND-SOURCE-OPERATORS** Reuse pure operator bodies from `lei_ren_part1_paper_compliant_pulse_end_physical_C2.py` and `lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator.py`; use `current_pulse_main_exit_background_tensor.py` as the closest pulse-exit decomposition reference. Rebind them to actual current-Rp raw source rows, incoming histories, five moments and original P0. Keep paired source products before enclosure. Acceptance: signed stress and separate remainder sectors at the supported point, source-AST replay of the full momentum decomposition, exact original viscosity/physical scaling; cone and flatness stay false until separately bounded.
- [ ] **CURRENT-RP-OFFCENTER-DIFFERENTIAL-COVERAGE** Prepare distinct nonzero Z and pulse main/exit/gap queries before freezing accepted views. Issue complete 144-row fields per delivery and report coefficient/total width independently. Acceptance: actual inverse delivery, all source identities and pressure FTC rows preserved at each query; no inferred uniform chart bound.
- [ ] **CURRENT-RP-TIME-SCALE-DYNAMICS** Introduce compatible cross-query bindings for several original log-time sections at a fixed similarity point. Preserve exact coordinate/time/amplitude functions; compare radial/axial scales and source-signed vorticity with directed ratio operations. Acceptance: independently checked scale ratios/exponents and separate Eulerian/material winding diagnostics. Time rescaling alone does not complete n-dependent recursion.
- [ ] **CURRENT-RP-PRESSURE-END-LAYER** Extend toward xi=13 with stable late-tail/expm1 arithmetic, a positive same-source late C0 lower bound, controlled derivatives and original P0/raw waiting root. Do not weaken the existing endpoint guard without its missing lower-bound proof.
- [ ] **CURRENT-RP-ALL-CHART-DIFFERENTIAL-FIELD** Extend measured-error-driven source evaluation to postpulse, flatten, steep, waiting, collar and exact heat exterior. Keep axis, seams and higher derivatives as separate bounded tasks; collect point and uniform scope separately.
- [ ] **CURRENT-GENERAL-GLOBAL-AXIS-FIELD** Admit independent xyz/t, solve the actual original inverse and handle radial support joins plus axis parity/Cartesian limits. Acceptance: nonzero off-axis evaluations, controlled axis limits and original pressure compatibility, without legacy amplitude/forcing substitutes.
- [ ] **CURRENT-SIGNED-STRESS-CONE-FLAT** After source operator admission, bound signed regional cone margins and separate flat remainder maxima/physical volume L2 norms. Include inner exit, annulus, pulse/end, flatten and heat collar. A stress/remainder decomposition alone does not prove admissibility or flatness.
- [ ] **CURRENT-N-DEPENDENT-RECOVERY** Implement distinct n=1 and n>=2 equations with a common original inner domain, independent moment repair, finite-order remainder and smooth sum. Truncate streamfunction/vector potential before curl. Acceptance: actual computed higher-order coefficients, moment and divergence identities, n-dependent residual/remainder evidence; coordinate scaling cannot satisfy this task.
- [ ] **CURRENT-OSCILLATORY-CORRECTION-FULL-NS** Implement both original pulse families and mean corrections, verify averaged quadratic stress cancellation, then check fixed-forcing full Cartesian NS residual and physical max/L2 norms. Keep the background and corrected fields separate in reports.

Previous checkpoint: [nonzero original scale-value arithmetic](CURRENT_ORIGINAL_RP_SCALE_VALUE_ARITHMETIC_2026_10_10.md). Continue from accepted live source operators and mark only demonstrated scope complete. The reconstruction remains a leading/background point interface; stress admissibility, flat remainder, true recursion and corrected NS completion remain open.
