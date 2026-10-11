# Genuine n=1 radial recovery and corrected common-inner geometry

The first coupled n=1 inner field now retains genuine first Z derivatives and recovers the radial velocity coefficient V1. Its initial integral retains second Z derivatives. A new source-bound API restricts this field to a corrected common inner interval below the original later-shear start. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before the earlier [first coupled step](CURRENT_ORIGINAL_N1_COUPLED_FIRST_STEP_2026_10_11.md). The earlier phase-1/2 computation remains valid for the canonical **unrepaired** source. It cannot serve as the completed-leading common inner interval because its outer radius exceeds the later shear modification start.

## Computed derivative and velocity extension

The actual leading mixed source rows now retain total radial/Z order four. N1theta, N1z, N1p and regular forcing g retain Z order two; B0/B1 retain Z order three. Genuine ordinary F0 derivatives and actual own cumulative-mean derivatives are retained. Missing derivatives are rejected rather than padded with zeros.

The actual initial integral Wbase has Z order two. The computed prefix Wbase+H Wbase, including pressure feedback 2F0F1 and B1 partial_Z Wbase, has Z order one. The same source logarithmic sectors are preserved throughout.

For temporal n=1, the exact recovery is

```text
mean1 = Uz1 + K1
L = 1-delta*Z^2
Q1 = [(1-delta)*Z*Uz1 - (1+delta)*Z*K1
      - (1-Z^2)*partial_Z(Uz1+K1)] / L
V1 = rho*Lambda^-1*Q1.
```

This uses q1=1+delta. Substituting the leading q0 recovery is incorrect. `n1_radial_recovery_prefix` exposes the computed mean1/Q1/V1; `n1_radial_recovery_series_enclosures` includes the omitted-series value and derivative bounds. All axis traces remain exactly zero.

The existing analytic tail bounds the omitted W series by E on Z radius eta/8. Reserving the additional Cauchy gap eta/16 yields a first-derivative bound E/(eta/16) on eta/16. The report exposes these radii separately in `derivative_tail`. The resulting Q1/V1 bounds concern the canonical inner series, with combined source/integration uncertainty. They do not establish useful component relative accuracy or a completed global n=1 field.

## Corrected common inner interval

The checked original strict collar selects a positive scalar sc below 1/4. Its original exact binary endpoint is restored; no new point is selected from a rounded enclosure. The same source family, implicit source, datum and original width ledger are required.

```text
Ra    = 4*epsilon_core
Rkeep = Ra*exp(hb*sc/8)
Rin   = Ra*exp(hb*sc/4)
r_-   = Ra*exp(hb*sc/2)

Ra < Rkeep < Rin < r_-.
```

The microscopic positive offsets remain formal. They are never added to a huge rounded logRa, and hb is never replaced by a numerical cap. A new actual coupled sample at Z=371/1000 and source phase sc/4 evaluates the corrected common-inner exit, including Q1/V1 recovery and its derivative tail.

The original Rm factor has basis power `(0,5,0,0,0)`, whose second logarithmic base is **2 logPstar**. Therefore

```text
log Rm = 10*logCstar + 10*logPstar + log110 - 6
log Ra = log4 - 4*logPstar - 1000
log(Rm/Rin) = 10*logCstar + 14*logPstar
              + log(110/4) + 994 - hb*sc/4 > 0.
```

The accepted original parameter/radius frame and the previously proved hb*sc<log(4.1/4) give the directed positive lower bound. This latter inequality is used only to bound the retained offset, never to replace the source radius. The first original Rm bump starts at (49/40)Rm, so it also lies beyond the selected Rin.

Accepted new gate: `current_original_n1_common_inner_geometry_below_known_modification_starts_selected`. This certifies the selection and source-bound geometry. **`completed_leading_supports_outside_Rin_certified` remains false**: the assembled completed-leading source does not yet provide an executable same-source zero-extension dispatch on this interval. Geometry alone cannot identify the completed-leading field with the canonical field.

## Files and API

The inner operator, initial integral and coupled-first-step quartets have been regenerated with their new derivative budgets. The new quartet in `experiments/root_st073/` is

```text
lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection.py
lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection.json.gz
lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection_check.py
lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection_check.json.
```

```python
from lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection import CurrentOriginalN1CommonInnerSelection

inner = CurrentOriginalN1CommonInnerSelection(accepted_coupled_owner)
packet = inner.evaluate('.371', '1', 'first')
view = inner.report(packet)
# First-chart coordinate is a FRACTION of sc/4, not the old source phase.
field = view['canonical_n1_field']
radial = field['n1_radial_recovery_series_enclosures']
assert field['ordinary_Z_prefix_derivative_order'] == 1
assert not view['completed_leading_supports_outside_Rin_certified']
assert not view['n1_full_inner_interval_solution_certified']
```

An already accepted coupled packet may be reused for the **exact same** Z/chart/source coordinate. Owner and coordinate checks reject a foreign or different query. This avoids repeating source computations solely to wrap a result in the corrected-domain API.

Focused checks cover 306 independent source/operator derivative rows; 36 independent initial-integral rows; 132 normalized-axis/range and signed coupled rows; nine independent mean/Q1/V1 polynomial rows; the extra Cauchy derivative margin, exact zero axis traces, source ownership and formal geometry. Prior ancestor suites were not rebuilt. The one reused read-only worker is **GPT-5.6 Luna / max**; it identified the support conflict and corrected the Rm logarithmic basis. No new child was spawned.

## Next executable tasks, in order

- [x] **CURRENT-RP-N1-Z-DERIVATIVE-BUDGET-EXTENSION** Extend actual mixed leading rows to total order four, forcing/g to Z2, matrices to Z3, Wbase to Z2 and the first coupled prefix to Z1. Preserve actual F0 derivatives and own means.
- [x] **CURRENT-RP-N1-RADIAL-VELOCITY-RECOVERY** Recover mean1/Q1/V1 using the n=1 temporal weight. Add the compatible Cauchy derivative tail and independently check the polynomial recovery and zero axis. Completion applies to the canonical inner coefficient series.
- [x] **CURRENT-RP-N1-CORRECT-COMMON-INNER-SELECTION** Bind original sc/source/width/parameter receipts, select Ra<Rkeep<Rin<r_-, prove Rm>Rin without materializing microscopic offsets, and compute an actual coupled sample at the corrected exit.
- [ ] **CURRENT-RP-N1-COMPLETED-LEADING-SUPPORT-GUARD** Implement the original leading-source branch below the corrected r_- in the assembled dispatch. Bind the later generic-shear branch, original Rm bump, later repair and outer matching supports to the same family/P0/histories. Return the canonical original leading source on the whole selected interval without changing the incoming cumulative moments. Only then promote the completed-leading support gate. Do not reuse the old phase-1/2 Rin.
- [ ] **CURRENT-RP-N1-COMPONENT-ERROR-BUDGET** Split source uncertainty, core integration, actual-collar integration and omitted series error. Preserve amplitude-dependent scales for F1/Uz1/Q1/V1/P1 and their derivatives; the uniform W-norm tail can swamp a smaller swirl component. Refine only the cells or source jets causing material uncertainty. Report widths and dominating terms.
- [ ] **CURRENT-RP-N1-CONTINUITY-AND-MEAN-IDENTITY** Bound the discrepancy between independently assembled radial integration of Uz1 and the recovered mean1; retain the exact regular-system identity rather than comparing unrelated interval midpoints. Include the core/first inlet and corrected Rin, pressure recovery and both signed feedback contributions.
- [ ] **CURRENT-RP-N1-INNER-EQUATION-DIAGNOSTICS** Assemble each original n=1 coefficient equation independently from the returned fields and required derivatives. Report source/integration/tail contributions over whole radial/Z cells. The first iterate itself is not an exact solved field; distinguish its defect from the infinite-series enclosure. Promote solved/full-inner gates only with the original domain, derivatives and support assumptions closed.
- [ ] **CURRENT-RP-N1-WHOLE-Z-AND-DERIVATIVE-COVERAGE** Deliver accepted query coverage for the whole real Z interval, including endpoints, on the same corrected radial interval. Continue ordinary derivative bounds from the common holomorphic tube, with explicit Cauchy losses and no zero padding. Additional point samples alone do not close this task.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** Extend the solved inner correction through the original annular construction, calculate its own five whole-Z defects, solve its independent original bump controls and preserve the same pressure. Truncate streamfunction/vector potential before curl. Verify terminal moment identities, joins/derivatives, axis regularity and finite-energy tail. Leading-order repair is not this repair.
- [ ] **CURRENT-RP-N2-AFTER-REPAIRED-N1** Form the actual n=2 quadratic sources from solved/repaired n=1 and leading data. Implement the true n-dependent recovery and regular equations on the same corrected interval. Solve and repair n=2 independently; do not substitute a radial Taylor degree or rescaling for the temporal order.
- [ ] **CURRENT-RP-RECURSIVE-SUM-AND-PHYSICAL-DIAGNOSTICS** After successive admissible temporal coefficients exist, bound the finite-order/time remainder and smooth summation. Measure radial/axial core widths, vorticity/swirl growth, material winding and energy. Preserve regional cone/flatness and oscillatory averaged-stress cancellation as separate requirements before full NS residual validation.

The two nearest unresolved deliverables are the completed-leading support dispatch and a useful separated field error budget. Their completion is required before claiming a repaired n=1 recursion layer or proceeding to n=2.
