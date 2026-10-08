# Active kappa derivative refinement and original weighted pressure

Accepted implementation: [0e8af43d](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/0e8af43df16e41943bfcb81886a6c239911ce2be). Predecessor: [CURRENT_FIRSTBRIDGE_PAIRED_BRANCH_C1_O2_REPLAY_2026_10_08.md](CURRENT_FIRSTBRIDGE_PAIRED_BRANCH_C1_O2_REPLAY_2026_10_08.md). Authoritative original source family, pressure datum, candidate N1024 and two strict-sign axial tiles are unchanged. Root implemented mathematics/code/compute; existing read-only reviewer **GPT-5.6 Luna / max** identified the original E_y relation and reviewed source provenance. No new child or Astra child was spawned.

Two executed advances now feed actual cumulative histories. Conditional kappa support reduces the firstbridge derivative overestimate. An independent exact radial-weight integral removes the inner_reference C0 pressure width artifact. Neither advance supplies terminal five-moment closure, a selected numerical field value, globally admitted N, stress admissibility or genuine coefficient recursion.

## Active kappa source derivative

```text
kappa=a+b²/a, Delta=kappa-2, a>0
Delta_Z=a_Z*(1-b²/a²)+2*b*b_Z/a
kappa<=K => a<=K, b²/a<=K, |b/a|<=sqrt(K/a)
|Delta_Z|<=|a_Z|*(1+K/a)+2*|b_Z|*sqrt(K/a)
```

Use K=2 on the negative branch and K=5/2 on the original transition branch (2+eta<=5/2). Only conditional same-source a_C0 and Delta_Z covers are refined in a copied root map. Original a_Z/b_Z, every other y/mixed root and the source functions remain. Flat support stays lazy without an active-branch constraint. Conditional q/q² jets and exact paired Poisson derivatives then drive the unchanged native density owner. The complete13-cell source route was re-evaluated on both tiles; all C0 contributions/backgrounds/P0 remained unchanged at this intermediate stage. The four accepted genuine original O2 integrals are replayed without recomputation.

For display only let L=10^(408906090034569676). Entries below are **log magnitude uppers divided by L**, not physical errors; never exponentiate these enormous bounds.

| Same original derivative cover | Accepted paired stage | Active kappa stage |
|---|---:|---:|
| Active-branch Delta_Z |2.91840497|1.44489658|
| Active-branch q_Z |6.10862216|4.63511377|
| O2 correction m_Z/e_Z |7.41045967|5.93695128|
| O2 correction h_Z |5.89403356|4.42052517|
| O2 correction k_Z |7.32462423|5.85111584|
| O2 correction p_Z |6.09431625|4.62080786|

All20 downstream derivative comparisons strictly improve. These ranges still cannot establish functional terminal identities or stress signs.

## Original pressure weight, not width times maximum

The original inner_reference theta amplitude is E=Utheta/Pstar. Its source log is affine in gap=log(R/Rsh), with slope1/10. Theta has no radial half shift or extra Pstar normalization. AST guards and inherited mixed-row receipts bind this exact source relation:

```text
E_y=E/10, y=log R
deltaE=E*expm1(A/N)
p_density=E*deltaE+deltaE²/2=E²*expm1(2A/N)/2
integral_left^right E² dy=5*E(right)²*(1-exp(-true_width/5))
|integral p_density dy| <= (5/2)*Emax²*(exp(2*Amax/N)-1).
```

The last inequality holds for variable original A(y,Z,phase) using its whole-period magnitude bound; it assumes no mean cancellation. The exact true width has already been integrated, so it is not multiplied a second time. The omitted factor1-exp(-width/5)<=1 only makes the cover conservative. The old independent maximum-times-width estimate lost the exponential E correlation across an enormous radial cell.

**Critical basis detail:** the final signed inner_reference source adapter overwrites the initial generic plain basis. Restore archived E/A with `(0,2*logPstar,0,2*logu,logR)` where logu is `source_provenance.amplitude_adapter.frozen_positive_source_log`. The original E has fourth source power1/2 and unit coefficient. The new checker requires that restored full E log cover equals its archived cover. Reusing firstbridge bases, or using a zero fourth basis from the earlier generic packet, is incorrect. No enormous radius, width or exponential is materialized. Only bounded Amax=60.6 and 2*Amax/N<1 are ordinary scalars.

The new pressure adapter changes only the **inner_reference pressure integral enclosure**, then sums every original rate-zero pressure contribution through all12 active upstream charts. Original background is added once and P0/P0_Z stay independent. Other four C0 corrections and all five Z corrections remain exactly the accepted active-kappa results. The original archives and their old pressure outputs are preserved; the compact weighted-pressure manifest is the authoritative refined pressure path and explicitly references those archives.

| Actual O2-output normalized pressure correction | Magnitude upper |
|---|---:|
| Z[.36,.38],256 cells |0.408698237|
| Z[.36,.38],2048 cells |0.408628958|
| Z[-.38,-.36],256 cells |0.408698237|
| Z[-.38,-.36],2048 cells |0.408648582|

Before this refinement all four pressure correction **log** uppers were about941541067348081227.34. Inner_reference alone now has log upper-3.00170285 (magnitude about0.04970); the remaining pressure budget is primarily Rh_reference, with actual_patch/restoration memory also retained. The table reports cumulative correction, not absolute pressure, physical NS residual, numerical moment defect, or evidence of blowup. Absolute pressure includes the unchanged original pressure background and separate P0. Pressure Z remains enormous.

## Reproducible evidence

Files in `experiments/root_st073/`:

- `lei_ren_part1_paper_compliant_current_native_upstream_active_kappa_C1.py/.json`: original active-branch derivative identity, actual source route and index of two complete source archives; same stem `_positive.json.gz`/`_negative.json.gz` retain all root/primitive/density/geometry/history evidence.
- Same stem `_check.py/_check.json`: all conditional root/q/q²/paired jets, unchanged original C0 source and all actual serial C1 rows;20 stricter O2 derivative covers.
- `lei_ren_part1_paper_compliant_current_original_O2_active_kappa_upstream_replay.json`: four accepted same-N O2 outputs after the derivative step, before the pressure weight step.
- `lei_ren_part1_paper_compliant_current_inner_reference_weighted_pressure.py/.json`: AST-bound radial weight theorem, correct signed bases, two integral proofs, two complete12-chart pressure paths and four final O2 outputs with separate/absolute pressure.
- Same stem `_check.py/_check.json`: independent symbolic antiderivative and exponential magnitude majorant, directed cap reconstruction, all24 pressure memory steps, independent accepted C1 affine operator, other four C0/all Z/P0 exact preservation and four pressure improvements.

Run the two new focused checkers against their saved manifests; no original upstream numerical producer or O2 quadrature rerun is needed. Interval context240 digits, ordinary arithmetic300 digits. Checks PASS. Git-index dependency audit: 1233 source/evidence hashes. Active-kappa producer/checker timings are recorded in their manifests; weighted-pressure producer 1.688s, focused checker 1.828s. Source gates are local and every full-construction flag remains false.

## Exact next six cells

The accepted O2_slope y[0,1] is the full slope chart: `current_native_O2_C1_histories.ROUTE` subdivides it into five adjacent cutoff intervals that cover the same full source domain. Its endpoint is the original O2_axial inlet. The remaining source route is:

| Ordered cell | Original chart coordinate | True logarithmic-radius map |
|---|---|---|
| axial |O2_axial [0,1]|common+exp(Md*x)|
| buffer_0_9 |O2_buffer [0,9]|common+exp(Md)+x|
| buffer_9_11 |O2_buffer [9,11]|common+exp(Md)+x|
| slope_mu |O3_slope_mu [0,1]|common+logP+x|
| power_0_1 |O3_power phase [0,1/Tw]|common+logP+1+Tw*x|
| power_1_2 |O3_power phase [1/Tw,2/Tw]|common+logP+1+Tw*x|

The power endpoint is Rc=Rw*exp(2). Reuse `NativeO2C1Histories`, `serial.serial_cell`, `NativeRcC1Histories.original_O3_cell`, `C1DuhamelOperator`/`append_true_cell`, and `NativeRcParameterTargets.geometry`/`exact_radius_maps`. Preserve source radius/Jacobian, separate P0 and all original seam bindings. The two O3_power cells have source-proved exact zero **own density** increments but must retain decay of inherited histories. Final record names include `actual_Rc_correction_C0/Z` and `actual_Rc_own_history_C0/Z`.

The existing function graph is a source-function definition/transport graph, not a numerical oracle: its `numerical_original_source_oracle_installed` and `actual_original_numerical_integrals_evaluated` gates are false. The narrow missing value interface must evaluate accepted chart roots, phase inverse/A/B densities and true-coordinate integrals at a single common integer N. It cannot return selected points of saved caps as original function values.

## Detailed next tasks

Record completion with source/receipt/results/commit, not a plan or inventory. Prefer extending actual source/transport/target/control interfaces over repeatedly refining conservative caps.

- [x] **ACTIVE-BRANCH-DELTA-Z-SUPPORT:** implement the original conditional derivative majorant from unchanged a_Z/b_Z, preserve all original mixed rows/flatness and feed actual serial histories. Acceptance met: both nonflat branches/tile and20 stricter downstream derivative covers.
- [x] **BRANCH-A-C0-POSITIVITY-CORRELATION:** intersect same-source a_C0 with a<=K and the original positive lower; keep original derivatives. The optional stronger joint b²/a<=K-a bound remains unimplemented.
- [x] **INNER-REFERENCE-WEIGHTED-PRESSURE-C0:** bind original E_y and exact pressure identity, use correct final signed amplitude bases, integrate true radial weight once, retain all inherited pressure/P0. Acceptance met: both genuine tiles and all4 O2 pressure outputs <=0.40870.
- [x] **NEXT-REFINEMENT-TO-O2:** attach both real upstream advances to unchanged original O2 source integrals at N1024 and both256/2048 levels; retain original family/datum/N/units. Genuine O2 quadrature was not rerun.
- [ ] **ACTUAL-O2-TO-RC-C1-TRANSPORT (priority1):** take all5 correction values/Z from `genuine_original_O2_weighted_pressure_replays` as actual incoming to the exact six cells listed above. Do not repeat/reset the accepted slope chart or its inlet. Reuse the exact original cell graph and native provider geometry. Acceptance: an executed same-N refined source-cover route to Rc, every interval retained, actual correction/own outputs in all5 units, original background once, separate analytic P0 and complete rate-zero memory. Distinguish conservative route covers from quantitative original target values; preserve power-cell inherited decay even though their own density is zero.
- [ ] **SIX-CELL-ORIGINAL-NUMERICAL-ORACLE (priority1):** implement the injected value/ordinary-Z provider for the exact six tail cells from accepted original chart roots and phase inverse/A/B density functions. Integrate with true-coordinate Jacobians and independent positive masses at one common N; retain complete native source/family/datum bindings. Acceptance: numerical source-derived signed contributions and reproducible integration errors/refinement, not chosen enclosure endpoints, cross-N cached modulation or a reused slope integral on a different domain.
- [ ] **CORRELATED-RC-TARGET-VALUES (priority1, after actual route):** evaluate the actual source-defined five terminal defects as functions of Z from the full route and original preheat/heat targets. Keep correlated source expressions, exact pressure datum and common N; report directed signed value/Z covers. Acceptance: real finite numerical target enclosures on executed domains, not substituted enclosure endpoints or symbolic graph roots alone.
- [ ] **ACTUAL-FIVE-INDEPENDENT-REPAIR-CONTROLS:** apply the accepted generic five-control inverse to those actual correlated targets, preserving source function dependency and positive denominator/implicit derivative proofs. Verify exact original recovery equations and all functional terminal identities. Acceptance: installed nontrivial same-family control functions and repaired histories, source-dependent derivatives, independent target/repair validation.
- [ ] **FIRSTBRIDGE-JOINT-ROOT-Z:** inspect common original formulas for a_Z/b_Z/p2_Z and Delta_Z before independent magnitude majorants. Retain signs/correlations through cutoff branches where justified; do not differentiate conditional selections. Acceptance: meaningful actual derivative improvement or a quantified limit with source identities; all original mixed rows/seams survive. This remains needed for useful derivative/implicit-control bounds.
- [ ] **NONLINEAR-BRANCH-DENSITY-HULL:** if primitive hulling dominates, pass each nonempty original branch through the exact nonlinear five-density C0/Z graph before a complete same-source union. Preserve density owner type and basis/context/ledger. Acceptance: source-correct branchwise densities and inherited outputs, actual narrower bounds, no derivative erasure or phase samples.
- [ ] **ORIGINAL-INNER-PRESSURE-Z-WEIGHT:** use E_Z=-2Z/(1+Z²)*E from the exact original axial amplitude together with genuine A_Z in `E*E_Z*expm1(2A/N)+E²*exp(2A/N)*A_Z/N`; integrate original E² weight with correct fixed-Z geometry. Acceptance: source/AST identity, full derivative-function cover and real downstream p_Z improvement, unchanged P0_Z. A C0 cap is never differentiated.
- [ ] **RH-REFERENCE-PRESSURE-BUDGET:** after the inner-reference width artifact is removed, assess original Rh_reference/actual_patch/restoration density correlation and true-phase mean integral. Acceptance: original source integrals with inherited memory and quantified pressure value/Z error. Do not infer period cancellation from local antisymmetry.
- [ ] **MICROSCOPIC-CROSSING-Z-ATLAS:** implement regular signed source derivatives across Z=0 where current strict-sign tiles do not cover. Respect original positivity/flatness, exact source basis and radius derivatives. Acceptance: continuous coverage and required mixed jets over the entire axial domain, not finite-point extrapolation.
- [ ] **GLOBAL-CANDIDATE-N:** admit one finite frequency using the actual full route/controls and all source inequalities/higher jets; local N1024 diagnostics are not global admission. Reject reused N7 modulation or phase evidence at another N.
- [ ] **MATCHED-PRESSURE-HEAT-BACKGROUND:** complete five functional moments, axis regularity, annular source joins, compatible analytic preheat pressure and exact/controlled heat exterior with finite-energy tail. Verify high-order joins and stage-appropriate errors.
- [ ] **ADMISSIBLE-STRESS-AND-FLAT-REMAINDER:** construct actual divergence-form stress plus separate remainder, cone signs/margins on every physical region and flat decay with scale dependence. A conservative cumulative bound or symbolic cone statement is insufficient.
- [ ] **GENUINE-N-DEPENDENT-COEFFICIENT-RECURSION:** only after leading/matching/stress acceptance, implement original n=1/n>=2 equations, common core interval, separate moment repair, divergence-preserving cutoff, finite-order remainder and smooth summation evidence. Coordinate rescaling or radial Taylor jets alone do not count.
- [ ] **TWO-FAMILY-OSCILLATORY-CANCELLATION:** implement actual mean and both pulse families, realizable averaged quadratic stress and cancellation, preserving time/space support and forcing regularity.
- [ ] **CORRECTED-CARTESIAN-UVW-AND-DYNAMICS:** recover computable [u,v,w](x,y,z,t) from completed dependent layers; independently measure Cartesian residual, L-infinity/volume-L2<1e-3, finite energy, core contraction/relative elongation, swirl/vorticity scaling and material winding. Distinguish measured dynamics from imposed scales.

The physical goal is still incomplete. The current pressure result clears one major numerical enclosure obstruction; it does not establish full functional matching or genuine scale recursion.
