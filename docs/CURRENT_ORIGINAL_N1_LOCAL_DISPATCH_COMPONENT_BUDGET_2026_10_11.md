# Original n=1 local leading dispatch and component error budgets

The corrected common inner interval now has an executable original leading-source dispatch. It preserves the original background velocity, pressure and cumulative moments while returning exact zero later-modification increments. A second API reports 16 source-sector component/derivative budgets for the actual first coupled n=1 prefix at the corrected exit. The full reconstruction goal remains **ACTIVE / INCOMPLETE**.

Read this checkpoint before [the radial-recovery checkpoint](CURRENT_ORIGINAL_N1_RADIAL_RECOVERY_2026_10_11.md). The new work advances source compatibility and exposes the actual numerical bottlenecks; it does not deliver a solved/repaired n=1 layer or n=2.

## Executable original local branch

The selected original geometry is unchanged:

```text
Ra < Rkeep=Ra*exp(hb*sc/8) < Rin=Ra*exp(hb*sc/4)
   < r_-=Ra*exp(hb*sc/2).
```

`CurrentOriginalN1LocalLeadingDispatch` binds the accepted original finite-N support selector and the actual selected positive eta from the same generic-shear source family. It replays `OriginalMicroFiniteN.modification_support` on the selected first chart, using the original sc. The whole selected source phase is at most sc/4, so the original unchanged-left branch applies. No active Section 11 inverse or ancestor source constructor is needed.

The original strict left-collar inequality and selected eta retain the exact flat-loop identity:

```text
Delta >= left_excess >= 2*eta
q=A=B=0
E_N=E, V_N=V.
```

This makes the later velocity/pressure/five-moment increments exactly zero for any positive N on this local branch. It does **not** set total velocity, P0 or incoming cumulative histories to zero. The returned leading source still carries its actual F0/Uz0/Q0/P0, own means, forcing and operator data. The core is a separate original stress-free branch; no Section 11 cone argument is applied at the axis. Directed original Rm/Rin and later Rc/Rin logarithmic gaps keep the repair supports beyond the selected interval.

The accepted new gate is `current_original_n1_selected_inner_local_leading_branch_identified`. Canonical holomorphic bounds transfer through this local identity. The older `completed_leading_supports_outside_Rin_certified` remains false: the complete global source dispatch, finite-N source graph and whole-Z incoming histories are still unfinished. The new local gate must not be substituted for that broader gate.

Four deliveries cover the whole-Z core axis, a core exit, the actual core/first inlet and the corrected actual inner exit. Independent checks exercise the four original support-selector branches, flat-loop velocity and five signed-density identities, nonzero total axis background/P0, and five ownership/state guards.

## Component budgets and current numerical bottlenecks

`CurrentOriginalN1ComponentBudget` reports the six regular state fields at Z orders zero and one, mean1 at both orders, and recovered Q1/V1 at order zero: 16 rows in total. Each row retains its original `(hb, F0_anchor, Pstar_squared, Lambda)` exponent tuple. It normalizes against a source monomial times a directed coefficient norm, not against a selected field value.

Equal source exponents cancel before combining enormous logarithms. Extremely small relative terms receive a directed small-log cover; excessively large positive relative tails remain logarithmic. Neither F0 nor an astronomical physical amplitude is expanded. Signed prefix enclosures and remaining-series tails are reported separately. Source uncertainty and integration uncertainty inside the prefix are still combined.

At Z=371/1000 and the corrected first-chart exit, **0/16 rows meet the 0.1% total relative-width target**. The conservative width `(upper-lower)/min(abs(endpoints))` is available only when the enclosure excludes zero.

| Component | Z orders | Prefix and total relative width, rounded | Current limitation |
|---|---|---:|---|
| Uz1 | 0, 1 | 0.204234 | Combined source/integration prefix width |
| K1 | 0, 1 | 0.397668 | Combined source/integration prefix width |
| P1 | 0, 1 | 0.133333 | Combined source/integration prefix width |
| d_x_Uz1 | 0, 1 | 0.133333 | Combined source/integration prefix width |
| mean1 | 0, 1 | 1.11125 | Source correlation/cancellation needs tighter bounds |
| Q1, V1 | 0 | Not finite in this metric | Signed prefix enclosure does not exclude zero |
| F1, d_x_F1 | 0, 1 | Not finite in this metric | Prefix uncertainty plus unweighted tail swamping the F0 sector |

These are conservative **enclosure widths**, not measured errors of the exact solution. The scalar absolute series tail can be tiny while still failing to resolve an F0-weighted swirl component. Therefore increasing arithmetic precision alone, or calling a small absolute tail “accurate swirl,” does not close the task. The budget makes two separate needs visible: tighter signed/correlated prefix bounds, and a sector-preserving analytic tail.

For real queried Z, the independent radial tail is

```text
mean_error <= 2E
Q_error <= [(abs(1-delta)+abs(1+delta))*abs(Z)*E
            + 2*abs(1-Z^2)*E1] / lower(1-delta*Z^2)
V_error <= rho*Lambda^-1*Q_error.
```

E and E1 use the already accepted common analytic tube and derivative gap. The Lambda exponent -1 remains explicit in the V1 tail. A complex-Z extension must bound abs(1-Z^2) on the admitted complex tube rather than reuse this real-point expression.

A directed-rounding defect was corrected in the new budget implementation: scalar `abs` and unary minus on high-precision endpoint objects could round them through the lower global scalar precision before interval conversion. Magnitude, sign and symmetric-tail operations now stay in the owner interval context. Five independent doubled-context rational prefix/tail containment checks include a negative non-dyadic endpoint; 48 independent radial extreme-error cases also pass. This is a focused new-feature check, not a replay of broad ancestor suites.

## Files and accepted API

Both quartets live under `experiments/root_st073/`:

```text
lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch.py
lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch.json.gz
lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch_check.py
lei_ren_part1_paper_compliant_current_original_n1_local_leading_dispatch_check.json

lei_ren_part1_paper_compliant_current_original_n1_component_budget.py
lei_ren_part1_paper_compliant_current_original_n1_component_budget.json.gz
lei_ren_part1_paper_compliant_current_original_n1_component_budget_check.py
lei_ren_part1_paper_compliant_current_original_n1_component_budget_check.json
```

```python
local = CurrentOriginalN1LocalLeadingDispatch(accepted_common_inner_owner)
local_view = local.report(local.evaluate('.371', '1', 'first'))
assert local_view['current_original_n1_selected_inner_local_leading_branch_identified']

budget = CurrentOriginalN1ComponentBudget(local)
view = budget.report(budget.evaluate(accepted_common_inner_packet))
assert view['current_original_n1_source_sector_component_error_budgets_computed']
assert view['passing_relative_total_width_rows'] == 0
assert view['total_component_derivative_rows'] == 16
assert not view['n1_full_inner_interval_solution_certified']
```

The accepted common-inner packet must come from the exact owner retained by `local`; a foreign packet or mutated live budget is rejected. Constructors require checked ancestors and unchanged source hashes. Raw generation intentionally does not promote acceptance; checked public construction/reporting returns both new gates true while all inherited global/solved/repair gates remain false.

The single reused read-only reviewer is **GPT-5.6 Luna / max**. No new child or GPT-6 Astra child was spawned.

## Ordered tasks for the next agent

- [x] **CURRENT-RP-N1-LOCAL-ORIGINAL-LEADING-DISPATCH** Bind the original finite-N selector, selected eta and strict collar; replay the unchanged-left branch on the corrected interval; preserve total backgrounds/P0/histories; keep core and open-loop domains separate. The completion scope is the new local gate only.
- [x] **CURRENT-RP-N1-SOURCE-SECTOR-BUDGET-DIAGNOSTICS** Publish signed source-sector prefix enclosures, component-tail scales, relative-width metrics and the actual 0/16 result. Keep full-inner/solved/repair gates false. This completes diagnostics, not useful numerical accuracy.
- [ ] **CURRENT-RP-N1-SECTOR-PRESERVING-TAIL** Derive a per-field, per-source-monomial majorant for Wbase and H, retaining exponent tuples through the operator and integral instead of injecting a scalar zero-sector error. Start with F1 and d_x_F1; include coupling from Uz1/K1/P1. Normalize only by original source-sector identities and directed norms. Deliver a checked tail at the corrected exit and explicit common analytic-domain assumptions.
- [ ] **CURRENT-RP-N1-WEIGHTED-AXIAL-GAUGE** If normalized raw Phi is used, retain `dZ(F0*Fhat)=F0*(Fhat_Z-Lambda*Gprime*Fhat)` exactly once. Existing `ordinary_mixed_profile_grids` already contains this gauge; do not add it again. The weighted B1 drift epsilon*Lambda*Gprime can be order Gprime, so the old unweighted small-contraction bound cannot simply be reused. Put the leading drift into the fundamental propagator or prove a compatible weighted majorant before accepting the new tail. Do not shrink the analytic tube to an epsilon scale merely to make an F0 ratio finite.
- [ ] **CURRENT-RP-N1-PREFIX-UNCERTAINTY-ATTRIBUTION** Separate actual core-source jets, core quadrature, first-collar source jets, collar quadrature and series tail. Preserve signed shared-source correlations. Identify which sources/cells create the 13.3%, 20.4%, 39.8% widths and the mean/Q cancellation interval. Refine those terms only; retain original width and amplitudes, not chosen interval midpoints.
- [ ] **CURRENT-RP-N1-COMPLETED-LEADING-SUPPORT-GUARD** Embed the local branch into the complete same-source completed-leading dispatch, including whole-Z histories, admitted finite-N source state, later repair and outer matching supports. Keep the existing broader support gate false until the assembled contract holds. Do not use the old phase-1/2 interval.
- [ ] **CURRENT-RP-N1-CONTINUITY-AND-MEAN-IDENTITY** Independently assemble radial Uz1 integration and compare with mean1 through the exact regular-system identity and correlated bounds. Cover the axis, core/first inlet and corrected exit, both signed feedbacks and the same pressure datum. Do not compare independent interval midpoints as an identity proof.
- [ ] **CURRENT-RP-N1-INNER-EQUATION-DIAGNOSTICS** Assemble the true temporal n=1 equations from returned fields and necessary derivatives. Separate finite-prefix equation defect from infinite-series enclosure; preserve pressure feedback and B1 partial_Z W. Deliver meaningful whole radial/Z-cell bounds before solved/full-inner promotion.
- [ ] **CURRENT-RP-N1-WHOLE-Z-AND-DERIVATIVE-COVERAGE** Continue accepted same-source derivatives on the corrected common interval for all real Z including endpoints. Use the common holomorphic tube and explicit Cauchy losses, reject missing derivatives, and keep actual incoming histories. Point samples alone do not complete this.
- [ ] **CURRENT-RP-N1-FIVE-MOMENT-REPAIR** After useful inner bounds and equations are available, extend n=1 through the original annulus, compute its own five whole-Z defects and independent original bump controls, and truncate a streamfunction/vector potential before curl. Retain the same pressure; check terminal identities, joins, regular axis and finite-energy tail. Leading-order repair does not substitute for this layer.
- [ ] **CURRENT-RP-N2-AFTER-REPAIRED-N1** Build actual n=2 quadratic sources from solved/repaired n=1 and leading data, implement the true n-dependent recovery on the same interval, and solve/repair it independently. Radial Taylor degree or rescaling is not temporal recursion.
- [ ] **CURRENT-RP-RECURSIVE-SUM-AND-PHYSICAL-DIAGNOSTICS** Only after successive admissible coefficients exist, bound finite-order/time remainder and smooth summation; measure core widths/aspect ratio, swirl/vorticity growth, material winding and energy. Regional cone, flat remainder and oscillatory stress cancellation remain separate requirements before full NS residual validation.

The nearest useful advance is a sector-preserving tail together with tighter, attributed signed prefix bounds. Do not claim completed scale recursion from the new local dispatch or budget API.
