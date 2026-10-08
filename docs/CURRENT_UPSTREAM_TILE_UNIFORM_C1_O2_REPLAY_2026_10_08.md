# Genuine pre-O2 axial tile refinement and accepted O2 replay

Checked source: [e255eb97](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e255eb97d0eb98176ed8063b7fb77721d9816291). Predecessor: [same-N actual upstream attachment](CURRENT_ORIGINAL_O2_SAME_N_ACTUAL_INCOMING_2026_10_08.md). Existing read-only reviewer **GPT-5.6 Luna / max**; no new child was spawned. Remote main f8e2b02e remains documentation-only relative to this numerical work.

The original sc/2-to-O2-inlet C0 and ordinary-Z correction histories have now been genuinely re-evaluated on **Z[.36,.38] and[-.38,-.36], N1024**. All13 original upstream cells remain present. These are new source covers on those tiles, rather than relabelled whole-Z[-1,1] bounds. Each of the12 active charts evaluates the inherited and accepted parameter-uniform whole-period primitive algorithms on the same original roots, eta, dstar and a lower. A complete certified A/A_Z/B/B_Z enclosure is retained from whichever bound is tighter. Each tile has **34/48 strict primitive magnitude-upper reductions**.

The four accepted genuine N1024 O2 y[0,1] value/ordinary-Z integral records (two signs x256/2048 refinements) are reused unchanged. The new actual upstream correction is propagated with each original own rate, then the original background is added once. Separate P0/P0_Z and pressure's rate-zero upstream memory survive. This improves the real incoming uncertainty without recomputing4608 O2 cells or changing any original source, phase, density, normalization or mass.

## Quantitative advance

Downstream propagated correction magnitude-upper ratios against the preceding whole-Z incoming attachment, positive tile /2048 source cells:

| Original channel | New / old magnitude upper | Approximate upper reduction |
|---|---:|---:|
| m |0.0147960|98.52%|
| h |0.3073050|69.27%|
| k |0.00989192|99.01%|
| e |0.2568862|74.31%|
| p |0.2565790|74.34%|

The negative/2048 ratios for h/e are0.3073016/0.2569826; m/k/p match. All40 downstream correction C0/Z comparisons across the four replay records strictly improve their logarithmic magnitude upper. These are **conservative enclosure diagnostics in the unchanged five units**, not reductions in the Cartesian momentum residual, observed physical amplitudes, or project completion percentages.

The ordinary-Z bounds remain astronomically wide. Let L=10^(408906090034569676) solely as display notation for the logarithmic bound scale. The following are **log magnitude uppers**, and must not be exponentiated:

| Channel derivative | Before: log upper / L | After: log upper / L |
|---|---:|---:|
| m_Z/e_Z |23.6047461|8.1472139|
| h_Z |22.0883199|6.6307878|
| k_Z |23.5189106|8.0613784|
| p_Z |22.2886026|6.8310704|

Pressure C0 correction still has log upper about941541067348081227.34. Its smaller upper remains far from a useful global closure estimate. No terminal identities, repair controls, globally admitted N, cone or recursion are certified.

## Actual chart attribution and next mathematical obstruction

The receipt retains every per-chart true-width signed C0/Z contribution and its formal magnitude bound. The **active_first_bridge** dominates the surviving upstream ordinary-Z uncertainty: its p_Z contribution has log upper6.8310704*L, unchanged at the O2 output because the pressure own rate is0. The m/h/k/e contributions are carried through their genuine downstream decay; none is reset. The separate P0_Z background has not been merged into this correction attribution.

The **inner_reference** contribution dominates pressure correction C0, with log upper941541067348081227.34. The C0 and derivative obstructions therefore require separate source refinements. Increasing O2 source-cell counts cannot resolve these upstream uncertainties.

Existing source-faithful branch-local machinery is available but is not yet wired into this12-active-cell owner:

- `current_native_cutoff_q_cover.conditional_cutoff_branches` enumerates each nonempty Delta<=0, 0<=Delta<=eta and Delta>=eta branch; `conditional_q_jet` retains original a/Delta ordinary derivative rows and requires strict positive a.
- `current_native_collected_q2_transport.conditional_q2_jet` evaluates original q² jets directly, cancelling transition eta/eta before interval evaluation while retaining linear q/q_Z and remaining inverse-eta factors.
- `current_native_paired_C1_transport.paired_C1_support` collects the exact paired identities `(qP)_Z=q_Z*(P+u*P_u)+q²*c_Z*P_u` and `(q²H)_Z=(q²_Z/2)*(2H+u*H_u)+q³*c_Z*H_u`, with c=p2/dstar and genuine nonzero p2_Z. The linear q_Z term remains.
- `current_native_signed_u_phase_cover.NativeSignedUDensityCover.spatial_query` and `current_native_cutoff_density_oracle` show how each conditional branch is carried through nonlinear density and then hulled by `same_source_union`. Conditional branches are not summed.

A mixed source box crossing Delta=0/eta cannot receive active-only a/b restrictions before the branch split. Aggregate q/q_Z rows cannot substitute for the correlated branch rows in the paired formulas. The current history chain enforces the exact `NativeDensityC1LocalIntegrals` owner type; a different paired-density subclass is not a drop-in replacement. A narrow original `whole_period_C1(roots,eta_log,log_a_lower,dstar_log)` adapter or an explicitly reviewed density hook must preserve that owner contract.

## Executable evidence and consumption

All paths below are relative to `experiments/root_st073/`:

- Producer `lei_ren_part1_paper_compliant_current_native_upstream_tile_uniform_C1.py`; manifest of the same stem `.json`.
- Complete positive/negative `_positive.json.gz` and `_negative.json.gz` routes: two lossless archives totalling 2,822,737 bytes. Every source/ordinary-jet/cutoff/primitive/geometry/contribution/inherited/background record is retained.
- Adapter `lei_ren_part1_paper_compliant_current_original_O2_tile_upstream_replay.py`; report of the same stem `.json`, with four complete factored incoming/correction/own/separate-pressure records and quantitative comparisons.
- Focused checker `lei_ren_part1_paper_compliant_current_native_upstream_tile_uniform_C1_check.py`; accepted receipt of the same stem `.json`.

Checker PASS:2 genuine continuous source tiles,13 cells/tile,12 active primitive selections/tile,120 full factored serial C0/Z rows/tile, independent accepted affine O2 replay on4 records, original backgrounds added once, true positive width/mass and rate-zero memory. It reads the saved complete routes at the actual native240-digit context and300-digit ordinary arithmetic. Accepted scalar implicit/parameter derivative proofs are reused. It does not rerun numerical ancestors or the O2 producers. Index closure audit: 1179 hashes. Producer runtime122.359s; focused checker/replay runtime1.672s.

```python
import json
import lei_ren_part1_paper_compliant_current_original_O2_tile_upstream_replay as source

saved = json.loads((source.HERE / source.NAME).read_bytes())
record = saved['genuine_refined_upstream_original_O2_replays'][1]
# Complete source ranges, not point values or a physical velocity evaluator:
correction = record['propagated_actual_correction_C0']
ordinary_Z = record['propagated_actual_correction_Z']
own = record['actual_original_plus_propagated_correction_C0']
```

Use these saved hashes and records for the next task. Recompute a native producer only for a new source enclosure algorithm, source/domain/candidate change, or an unresolved failure. Do not rerun unchanged original O2 contributions to check a changed incoming bound.

## Detailed forward task queue

- [x] **GENUINE-UPSTREAM-STRICT-SIGN-TILES:** re-query every original pre-O2 cell on both genuine axial tiles at N1024; retain the true source-defined flat collar and all12 active cells.
- [x] **UNIFORM-PERIODIC-C1-UPSTREAM:** evaluate accepted parameter-uniform periodic jet bounds and inherited bounds on identical roots; select a complete enclosure without redefining a source, clipping derivatives or differentiating the interval selector. Save all48 decisions/tile.
- [x] **REFINED-ACTUAL-O2-REPLAY:** attach new incoming correction C0/Z covers to all four accepted same-N O2 integral records; retain common factor coordinates, background once, P0 separately and rate-zero pressure memory.
- [x] **QUANTITATIVE-CHART-LOCALIZATION:** retain per-chart true-width contribution attribution and name the firstbridge ordinary-Z and inner-reference C0 pressure bottlenecks. Record40 strict downstream bound reductions with their actual units.
- [ ] **FIRSTBRIDGE-CUTOFF-CORRELATED-BRANCHES (priority1):** in a new source adapter, enumerate every nonempty Delta body/transition/flat branch of the genuine whole firstbridge source box. Keep each restricted Delta C0 with its unchanged original derivative rows and strict a lower. Obtain original conditional q and direct q² ordinary jets with shared eta/dstar/context/family/ledger. Apply active a/b restrictions only after split with their original support proof. Acceptance: complete coverage across Delta0/eta with retained smooth seams and no arbitrary quiet/point branch; explicit branch-source evidence.
- [ ] **FIRSTBRIDGE-PAIRED-PERIODIC-Z (priority1):** use accepted paired_C1_support on each such branch, including genuine linear q_Z and nonzero p2_Z. Hull complete A_Z/B_Z function covers across branches; compare against this stage on the identical source tile. Do not apply paired formulas to aggregate mixed-box q rows. Acceptance: strictly tighter actual firstbridge derivative contributions at fixed N/source, or a recorded non-improvement identifying the next retained term. Preserve the exact density owner contract and original C0 functions.
- [ ] **NONLINEAR-BRANCH-DENSITY-CORRELATION:** if a primitive hull loses the needed q/A/B relation, carry each branch through the original nonlinear five C0/Z density graph before `same_source_union`. Derive an explicit reviewed hook rather than silently replacing the enforced density-owner class. Acceptance: complete same-source coverage and genuine true-width contributions without summing branches or duplicate widths/Jacobians.
- [ ] **FIRSTBRIDGE-PHASE-CONDITIONING:** use the original fractional true-radius phase unions within the dominant bridge; retain tiny positive widths and original periodic inverse derivatives. Refine the dominant source/phase dimension rather than every chart. Acceptance: a narrower genuine derivative contribution and a measured total upstream/O2 improvement.
- [ ] **INNER-REFERENCE-C0-PRESSURE:** isolate the original inner_reference pressure density terms and keep shared E/axial/loop source factors jointly before cancellation and unit export. Treat separately from P0/P0_Z and the firstbridge derivative task. Acceptance: a strictly reduced real pressure-correction bound without changing the analytic datum or deleting the active reference interval.
- [ ] **INDEPENDENT-BACKGROUND-P0-Z:** query/refine the original packet pressure datum and background derivatives only if their full own-pressure budget dominates. q² loop correlation cannot be applied to the independent datum. Acceptance: source-faithful background/absolute-pressure comparisons with original units and independent provenance.
- [ ] **FIXED-N-REFINED-UPSTREAM-O2-ATTACHMENT:** after a new upstream refinement, replay these same genuine O2 contributions and report C0/Z/correction/own/absolute pressure effects. No N7 splice, incoming zero fallback or redundant O2 computation. Acceptance: independent affine composition and lossless source evidence.
- [ ] **MICROSCOPIC-CROSSING-Z-ATLAS:** factor odd p2=Z*Q with even Q and genuine nonzero p2_Z. Cover overlapping negative/regular-|u|<=1/4/positive source coordinates for u=p2*q/dstar, retaining every nonempty source region near Z0. Acceptance: explicit coverage, original ordinary jets and no numerical threshold deleting microscopic sources.
- [ ] **CROSSING-Z-CORRELATED-DENSITY-INTEGRALS:** derive joint exact-midplane/near-midplane derivative expressions before finite export. Acceptance: real C0/ordinary-Z source integrals on a crossing tile, retaining all implicit terms and formal signed bounds wherever finite export is unproved.
- [ ] **WHOLE-AXIAL-SOURCE-ATLAS:** extend genuine strict-sign and crossing tiles through the original required axial domain. Hull overlaps of the same functions without summing duplicate contributions. Acceptance: explicit coverage/unresolved map and compatible overlap bounds.
- [ ] **ORIGINAL-O2-TO-RC-TRANSPORT:** continue from O2_slope y1 through O2_axial/O2_buffer, O3 and actual pulse/end/flatten/heat cells with true maps, masses and inherited history. Acceptance: every original interval/seam covered at one candidate N, P0 kept separate.
- [ ] **CORRELATED-TERMINAL-FIVE-TARGETS:** evaluate the exact Rc graph on the complete source route. Keep C=D_k-A*D_m and C_Z jointly before division by mu*A². Acceptance: quantitative properly normalized functional targets on the axial atlas; graph identifiers or point fits alone are insufficient.
- [ ] **FIVE-INDEPENDENT-REPAIR-CONTROLS:** solve the original support-specific bump/Gram response system from genuine targets, preserving source regularity/nonlinear errors/pressure compatibility. Acceptance: all five terminal identities as axial functions with independent conditioning and error bounds.
- [ ] **ONE-GLOBAL-ADMITTED-N:** combine complete-route amplitude, repair, cone and required derivative inequalities; recompute phase-dependent source integrals if the candidate changes. Acceptance: one finite integer satisfying every original condition. N1024 here remains a candidate.
- [ ] **HIGHER-MIXED-JETS-AND-MATCHED-HEAT:** recover actual required radial/higher ordinary-Z derivatives and inner/annular/pulse/flatten/heat joins, axis regularity, analytic preheat pressure, exact heat exterior and finite energy. Acceptance: functional matching at the required derivative order.
- [ ] **ADMISSIBLE-STRESS-AND-FLAT-REMAINDER:** form original divergence stress T_B and separate E_B from the matched background; measure true global cone margins, high-order flat decay and norm/scale budgets.
- [ ] **TRUE-N-DEPENDENT-COEFFICIENT-RECURSION:** implement distinct n1/n>=2 recovery, common inner domain, independent repair, finite-order remainder and smooth sum. Acceptance: actual order-dependent recursive coefficients/error control, not scaled copies of a field.
- [ ] **TWO-FAMILY-OSCILLATORY-CANCELLATION:** implement original mean corrections and both pulse families; verify averaged quadratic stress realization/cancellation and remainder budgets.
- [ ] **CORRECTED-CARTESIAN-UVW-AND-DYNAMICS:** export [u,v,w](x,y,z,t), independently check forced NS Linfinity/volume-L2 residual<1e-3, finite energy, radial shrinkage/axial aspect, swirl/vorticity exponents and material winding.

For each completed task, check its box and link the actual code/result/receipt/commit with measured changes. A failed improvement is useful evidence but does not complete its stated numerical acceptance condition. Full-route functional closure, stress, actual scale/coefficient recursion and corrected Cartesian field remain open.
