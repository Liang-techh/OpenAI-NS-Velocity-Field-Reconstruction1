# Full pre-repair original n1 moment-value enclosures — 2026-10-13

All five original n1 moments now have a common, source-bound value interface on every real Z in [-1,1]. The actual cutoff adapter has been evaluated on four full axial cells. Its three velocity moments remain signed factored physical sectors with their original value/Z1 budgets. Both pressure moments now include the actual inner cutoff terms plus the complete original external integrals, represented by 512 continuous common-grid models. Maximum absolute pressure errors in the fixed positive source units remain 1.7214685424e-5 / 2.3006530659e-5.

These are the complete **pre-repair value enclosures**, not five terminal zeros or selected coefficient functions. Full pressure derivative accuracy, the positive-order repair map, repaired own M/V/P, true n2/temporal recursion, global stress and final NS validation remain open. Full original goal **ACTIVE / INCOMPLETE**.

## Actual inner and external composition

The checked live `CurrentOriginalN1RadialCutoffExtension` is queried on [-1,-.5],[-.5,0],[0,.5],[.5,1], at coordinate "1", chart "first", with four radial cutoff cells. The independent pressure datum, actual selected common-inner phase, signed inner solution and its error enclosures are inherited from the same checked source family. F1/Uz1 are cut first; their own M, V and P are recovered afterward. No leading-background pressure or old unmodified mean is substituted.

The output `physical_cutoff_moment_values_at_Rin` supplies all five physical source-sector values. Each axial cell retains 2 source sectors for m11, 17 for m12 and 45 for m14. The pressure cutoff outputs retain 19 sectors for m13 and 27 for m15. Each record includes its original four source exponents, signed ordinary Z derivatives and combined source logarithm. A zero coefficient is removed only when its interval is exactly zero; no uncertain value or moment is set to zero.

The complete pressure sums are:

```text
m13_full = m13_cut(Rin) + external_m13
m15_full = m15_cut(Rin) + external_m15
external_m13 = integral_Rin^infinity N1p0 dR
external_m15 = -integral_Rin^infinity R*N1p0 dR
```

The external m15 already carries its original minus sign. It is not negated again during composition. The external source is the accepted [complete pressure provider](CURRENT_ORIGINAL_N1_EXTERNAL_PRESSURE_2026_10_13.md); all original radial segments are included once, including exact meridional support beyond Rv.

For each pressure moment, the inner signed sectors are compared with the same positive unit Phi*A*Rv^p, p=1,2. Its logarithm is formed using the original correlated source identities, before enormous physical exponentials are taken. In particular the inverse-mu coefficient is 13*p-15/2, including the -1/mu term in logE. On every whole axial cell, every nonzero inner sector's normalized magnitude logarithm is below -1000. Its original sign is retained with a directed nonzero exp(-1000) allowance. The resulting pressure values are signed interval sums, not chosen zeros or affine midpoint fields.

Each model keeps the external affine coefficients and both inherited Taylor remainders. The inner cutoff contribution is added as a separately identified nonzero source error. These affine coefficients describe the external part only. They do not approximate the full pressure derivatives or define the physical repair functions. Three velocity moments are left in their correct signed physical source units rather than collapsed through impossible floating-point exponentials.

The unchanged parent cutoff source record still states that its external tails were uncomputed; that historical record is preserved. The new composed provider explicitly supplies those tails. Its current scope is the complete radial support for all five pre-repair values.

## Reproduction and completed work

Run `python scripts/run_original_n1_five_defect_values.py`. It restores checked owners once, follows the original constructor type guards, computes the four genuine whole-Z cutoff packets and combines the checked external models. It does not run ancestor producers/checkers or overwrite their reports. The public provider is `CurrentOriginalN1FiveDefectValues(checked_cutoff,checked_external_pressure)`, with `evaluate(Z)` and `report(packet)`. Accepted loading reuses the saved source cells; it does not recompute them.

The new source quartet is `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_five_defect_values`, with `.py`, `.json.gz`, `_check.py`, `_check.json`. The only installed gate is `current_original_n1_full_five_defect_value_enclosures_installed`.

- [x] Restore the checked original cumulative/cutoff owners without ancestor suite reruns.
- [x] Evaluate actual cutoff moments on all four whole axial cells.
- [x] Preserve signed velocity source sectors and ordinary Z1 budgets.
- [x] Bind actual cutoff pressure sectors to the same independent P0/source family.
- [x] Add both complete external pressure tails with the correct m15 sign.
- [x] Prove each retained inner pressure allowance in the original positive units; preserve its sign and nonzero uncertainty.
- [x] Extend all 512 continuous pressure models with the actual inner contribution.
- [x] Independently check units, sector signs, all models, endpoint/seam/interior query sums and scope.
- [x] Confirm the accepted public interface at default precision with source precision restoration.

## Correct positive-order repair route

The current Rc nonlinear equation `B*h+d+Q(h)/N=0` is a leading-background finite-N control mechanism. It is not the positive-order n1 Section 15 repair and must not be reused as that repair without a separate derivation.

The Section 15 linear repair uses a region where the actual leading axial velocity is zero and the leading angular velocity is Astar(Z)*(2R)^(-1/2-mu). For bumps normalized by integral b(R)*dR/sqrt(2R)=1, the axial pair has rows (1,(2R)^(-mu)) and the angular triple has rows ((2R)^(1/2),(2R)^(-3/2-mu),(2R)^(-1/2-mu)). Its equations are:

```text
Bz * alpha = -[m11_full, m14_full/Astar]
Btheta * beta = -[m12_full, m13_full/(2*Astar), -m15_full/Astar]
```

The third angular right-hand-side component therefore becomes positive m15_full/Astar after the outer minus sign. These equations use genuine signed defect functions, not magnitude-majorant ranges or chosen model centers. C1 defects/amplitude/bump functions are the minimum for first derivatives; C3 is additionally required to claim C3 repaired fields.

A source inspection identified a useful exact gap between the existing leading repair bumps. Their log(R/Rc) supports are [7,9]*log(2)/40, [19,21]*log(2)/40 and [31,33]*log(2)/40. The reserved interval [11,17]*log(2)/40 is strictly outside all three supports. Accepted `compact_raw_beta` source nodes and their lazy derivative nodes vanish there, so delta_E=delta_V=0 and all ordinary Z derivatives of those corrections vanish, independently of their finite-N coefficients. The base axial velocity is zero and the original angular power law remains there. Accumulated histories, pressure and density are not zeroed.

Evidence is in `current_original_whole_Z_limit_repair_band.py`, `current_native_repair_band_cone_budget.py`, `current_original_whole_Z_all_N_control_connection.py` and `current_limit_band_physical_recovery.py`. Before installing an n1 repair, bind the physical angular amplitude of the actual original family: if its physical convention is S*A_rc(Z)*(R/Rc)^(-alpha), then Astar=S*A_rc(Z)*(2Rc)^alpha. S and A_rc must come from that actual source, not a free normalization.

The paper's prescribed old centers/shapes do not all fit this single gap. Narrower new n1 supports are an **autonomous support choice**, permitted for the staged reconstruction but requiring fresh bump normalization and determinant proofs. Do not call them the paper's original fixed bump set. The subsequently checked [linear operator](CURRENT_ORIGINAL_N1_LINEAR_OPERATOR_2026_10_13.md) uses log centers 12,13,14,15,16 times log(2)/40 and half-width log(2)/400; all five supports are strictly disjoint inside the reserved gap. Its support and normalization certificates are now installed; actual repaired coefficient functions remain open.

## Next repository tasks

- [ ] Bind complete signed defect function handles and their genuine Z derivatives, retaining the source graph, N, mu, P0 and radius definitions. The current value interface is not a substitute for this function adapter.
- [ ] Supply full pressure derivative enclosures from the original source integrals/endpoint equations; do not differentiate magnitude bounds or affine approximation errors.
- [x] Install a same-family reserved-gap source certificate for the actual leading velocities and their necessary derivatives; retain all nonlocal histories.
- [x] Bind actual physical Astar and Astar_Z, prove positivity on the required domain and preserve its logarithmic scale.
- [x] Choose and record the new n1 bump supports explicitly; verify they are inside the exact gap and mutually disjoint.
- [x] Compute each new positive normalization integral with physical dR/sqrt(2R), retaining its Rc powers and nonzero quadrature errors.
- [x] Construct the original axial two-by-two and angular three-by-three matrices from those actual normalized bumps.
- [x] Prove determinant/invertibility bounds from ordered positive supports and the distinct exponential powers; do not assume the old leading matrices apply.
- [ ] Evaluate the axial divided-difference row stably when mu is tiny. Preserve the signed right-hand-side correlation before dividing by mu; do not subtract independently boxed nearly equal rows.
- [ ] Construct alpha(Z) and beta(Z) as linear operators on the genuine signed full defects. Interval coefficient ranges are diagnostics, not a choice of fields.
- [x] Check the paper moment-map signs with a separate nonunit finite-scale fixture, especially m13's factor two and m15's sign.
- [ ] Add the new corrections to the correct F1/Uz1 family and recover its own M/V/P on the complete physical support.
- [ ] Verify all five terminal identities as full axial functions; keep repair closure false until these identities and source bindings are established.
- [ ] Check axis regularity, high-order joins, exact support and analytic preheat/heat compatibility of the repaired family.
- [ ] Install n1 repaired fields and their measured errors, explicitly separating any remaining full derivative limitation.
- [ ] Implement genuine n2 recovery on the same original inner interval, followed by its own cutoff and five-moment repair.
- [ ] Continue higher-order/temporal coefficient recursion, finite-order remainders and controlled smooth summation.
- [ ] Establish admissible stress and flat-remainder behavior across inner exit, matching, pulse/end, flatten and heat regions.
- [ ] Add both oscillatory pulse families and mean corrections; verify averaged quadratic momentum-flux cancellation.
- [ ] Independently validate the complete corrected Cartesian NS residual against a force restricted in advance, in physical L-infinity and volume-L2.
- [ ] Mark each completed task `[x]`, record its exact source/receipt/commit and limits, and update all five entry documents.

Full goal **ACTIVE / INCOMPLETE**. Complete pre-repair five-moment values are available; genuine signed function/derivative repair inputs and the positive-order coefficient map are next.
