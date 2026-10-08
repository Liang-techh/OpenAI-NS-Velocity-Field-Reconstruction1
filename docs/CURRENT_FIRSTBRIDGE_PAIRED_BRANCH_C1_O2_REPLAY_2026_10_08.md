# Genuine firstbridge paired cutoff-branch derivatives and O2 replay

Checked source: [da661d09](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/da661d093a8d0d43ccff6fb549bb11b8720a4cfd). Predecessor: [genuine upstream tile refinement](CURRENT_UPSTREAM_TILE_UNIFORM_C1_O2_REPLAY_2026_10_08.md). Read-only reviewer **GPT-5.6 Luna / max**, existing agent only. Remote main f8e2b02e remains documentation-only relative to this numerical work.

The actual firstbridge source on Z[.36,.38]/[-.38,-.36], N1024 now uses **complete original cutoff branches, direct q² ordinary jets and paired Poisson derivative identities** before the five density derivatives. The entire13-cell pre-O2 route is re-evaluated; the other11 active cells retain their accepted uniform primitive covers. Original C0 A/B objects, every C0 chart contribution, original backgrounds/P0 and the exact native density-owner type remain unchanged.

The firstbridge A_Z and B_Z/Pstar enclosures strictly tighten on both tiles. Those real source bounds feed the native densities, true-width masses and all subsequent inherited histories, then all four accepted original O2 integral records are replayed unchanged. Every one of the20 downstream ordinary-Z correction comparisons improves against the preceding uniform-tile stage. C0 corrections are exactly unchanged. Full Rc identities/global N/stress/actual scale recursion remain open.

## What changed mathematically

Every nonempty original Delta<=0, 0<=Delta<=eta and Delta>=eta branch is retained, including smooth seams. Only each conditional Delta C0 range is restricted in this stage; all original a/Delta derivative rows stay. `conditional_q_jet` supplies original linear q/q_Z, and `conditional_q2_jet` computes q²/q²_Z directly, cancelling transition eta/eta before enclosure arithmetic. This does not set a source derivative to zero; only the original flat branch has exact zero q and jets.

With P=h^-1*W1, H=h^-2*W2, u=(p2/dstar)*q, c=p2/dstar:

```text
(qP)_Z = q_Z*(P+u*P_u) + q²*c_Z*P_u
(q²H)_Z = (q²_Z/2)*(2H+u*H_u) + q³*c_Z*H_u
```

The original paired bounds `|P+uP_u|<=1`, `|2H+uH_u|<=2pi+3` replace independent q*u_Z products. Linear q_Z and genuine p2_Z remain. Conditional active a/b support bounds are applied only inside the active branch. Complete A_Z/B_Z covers are hulled across branches, then a full inherited or paired derivative enclosure is selected by its certified magnitude upper. No derivative of an interval selector/C0 cap is used. This stage still hulls primitive branches before nonlinear density; it does not claim branchwise nonlinear correlation or full-period cancellation.

## Quantitative result

Let L=10^(408906090034569676) be display notation for the **log magnitude-upper scale**. The following are logarithmic enclosure bounds divided by L. They must not be exponentiated or described as physical momentum-error reductions.

| Actual derivative cover | Accepted uniform stage | Paired branch stage |
|---|---:|---:|
| Firstbridge A_Z |8.2902729|7.5535187|
| Firstbridge B_Z/Pstar |9.7208636|8.9841094|
| O2-output correction m_Z/e_Z |8.1472139|7.4104597|
| O2-output correction h_Z |6.6307878|5.8940336|
| O2-output correction k_Z |8.0613784|7.3246242|
| O2-output correction p_Z |6.8310704|6.0943163|

Both signs and both256/2048 O2 levels have these same formal derivative upper scales. They remain enormous. Pressure C0 correction remains unchanged, log upper about941541067348081227.34, dominated by inner_reference. The firstbridge continues to dominate ordinary-Z correction; rate-zero pressure retains its complete memory.

## Sharper remaining frontier

Saved firstbridge evidence exposes why the derivative bound is still wide. The unconditioned original b and p2 C0/Z logarithmic uppers are about0.0143059*L, t0_Z about2.9040991*L, Delta_Z about2.9184050*L. After conditional q reconstruction, q_Z has log upper6.1086222*L on both active branches. Although the cutoff branch has narrowed Delta C0 and the paired support has narrowed a/b C0, the q-jet consumer still receives the unconditioned Delta_Z enclosure.

The next source-equivalent refinement follows the exact original identity

```text
kappa = a + b²/a,  Delta = kappa-2,  a>0
Delta_Z = a_Z*(1-b²/a²) + 2*b*b_Z/a.
```

On the negative branch kappa<=K=2; on the transition branch kappa<=K=2+eta<=5/2. Therefore a<=K, b²/a<=K and |b/a|<=sqrt(K/a), so

```text
|Delta_Z| <= |a_Z|*(1+K/a) + 2*|b_Z|*sqrt(K/a).
```

This bounds the **same original derivative function** using the unchanged source a_Z/b_Z and accepted positive a lower. It is not a derivative of the conditional selector. It can replace a wider Delta_Z cover in a copied branch-local root map before conditional q/q² jets; the accepted root itself must not be mutated. No active constraint applies to the flat interior Delta>eta. The original eta is globally constant, eta_Z=0; cutoff flatness preserves the Delta0/eta seam derivatives. `current_generic_loop_point_Z.py` evaluates the exact point identity, but no accepted conditional interval helper yet implements this majorant.

For pressure C0, an independent exact source correlation is available:

```text
x=A/N, deltaE=E*expm1(x)
p_density=E*deltaE+deltaE²/2 = E²*expm1(2*x)/2
p_density_Z=E*E_Z*expm1(2*x)+E²*exp(2*x)*A_Z/N.
```

The current density kernel forms its two pressure terms separately. Existing bounded-primitive/q² adapters do not implement this single pressure factorization. A pressure-only cover adapter can preserve the entire E factor and original rate-zero mass; P0/P0_Z are independent and remain separate. Phase antisymmetry alone cannot eliminate the N^-1 pressure term from a local hull: cancellation requires an actual full-period true-phase integral, which this stage has not performed.

## Executable source evidence

Paths are relative to `experiments/root_st073/`:

- `lei_ren_part1_paper_compliant_current_native_upstream_paired_branch_C1.py`: serial runtime adapter for the first primitive call in the exact original active chart order; context bindings are restored in finally. It evaluates original negative/transition/flat q/q² rows, paired supports and complete hulls, then re-runs the genuine native route.
- Same stem `.json`: source family/datum/N, exact13-cell incoming roots, paired-vs-uniform primitive comparisons and two source archive indexes.
- Same stem `_positive.json.gz` and `_negative.json.gz`: complete source/primitive/cutoff/geometry/density/history/background records, losslessly saved, total3,073,394 bytes. Firstbranch records include all native source roots and bases for exact saved-jet reconstruction.
- Same stem `_check.py` and `_check.json`: focused saved-source reconstruction and affine transport receipt. Original scalar paired-parameter/q² checks are reused.
- `lei_ren_part1_paper_compliant_current_original_O2_paired_branch_upstream_replay.json`: four complete factored same-N O2 outputs, corrected bindings, separate/absolute pressure and comparisons against the accepted uniform-tile replay.

Checks PASS: all3 cutoff branches/tile,42 reconstructed q/q²/paired rows/tile, entire original13-cell route, unchanged original source/geometries/backgrounds/C0 contributions on all12 active cells, full inherited C0/Z affine composition and original background once, all20 downstream ordinary-Z bound improvements. Native interval context240 digits, ordinary arithmetic300 digits. Existing accepted scalar reference proofs are not rerun. Index dependency audit: 1222 hashes. Producer127.047s; focused checker/replay2.297s. O2 integrals are unchanged and not recomputed.

## Detailed next actions

- [x] **FIRSTBRIDGE-CONDITIONAL-CUTOFF-Q-Q2:** retain every original cutoff branch/seam and original ordinary derivative, obtain direct q² rows from the same roots, eta/dstar/positive a theorem and context.
- [x] **FIRSTBRIDGE-PAIRED-PERIODIC-Z:** use exact paired identities branchwise, preserve linear q_Z/nonzero p2_Z, hull complete derivatives and keep the original C0/native-owner contract.
- [x] **ACTUAL-SERIAL-C1-AND-ORIGINAL-O2-REPLAY:** propagate genuine improved firstbridge derivatives through all13 upstream cells and all4 accepted same-N O2 records; preserve unchanged C0/background/P0/rate-zero memory.
- [ ] **ACTIVE-BRANCH-DELTA-Z-SUPPORT (priority1):** implement the displayed kappa derivative majorant on a copied nonflat branch root map, from unchanged original a_Z/b_Z and the actual a positive lower. Preserve basis/context/ledger; select a complete valid tighter Delta_Z cover. Keep original y/mixed derivative rows and flat branch laziness. Acceptance: exact source identity and full branch coverage, no selector derivative or source clipping, actual reduced conditional q_Z/primitive/density/history-Z bound at fixed source/N. Record any retained term that still dominates.
- [ ] **BRANCH-A-C0-POSITIVITY-CORRELATION:** use a<=K and the original positive a lower to narrow its C0 cover only within the branch. If using b²/a<=K-a, retain the positive correlated expression rather than subtract independent upper endpoints. Acceptance: valid conditional a/b ratio bounds, unchanged a_Z/b_Z and smooth seam coverage.
- [ ] **SIGNED-DELTA-Z-RANGE:** where the original joint a_Z/b_Z expressions permit, preserve signs before taking the magnitude majorant; use a full certified signed cover. Acceptance: genuine derivative enclosure improvement, not chosen endpoint/root.
- [ ] **NONLINEAR-BRANCH-DENSITY-HULL:** if primitive hulls erase useful q/A/B/E correlations, carry each branch through the original nonlinear five C0/Z density graph before same_source_union. Preserve the enforced density owner type via a reviewed hook. Acceptance: all branches covered, complete signed densities, width/Jacobian applied once and a reduced actual inherited bound.
- [ ] **FIRSTBRIDGE-TRUE-PHASE-CONDITIONING:** refine the original true-radius fractional phase union/source domains only where the remaining firstbridge derivative budget dominates. Acceptance: genuine original inverse derivatives and continuous source coverage, no phase samples replacing integrals.
- [ ] **INNER-REFERENCE-PRESSURE-IDENTITY (priority1):** evaluate E²*expm1(2A/N)/2 in one signed factored source expression and compare with the accepted original pressure cover before its true rate-zero mass. Keep E powers, native units and P0 separate. Acceptance: exact source identity, genuine narrower pressure C0 contribution/output or a quantified non-improvement. Do not infer whole-period cancellation from local hulls.
- [ ] **PRESSURE-ORDINARY-Z-IDENTITY:** after C0 pressure factorization, evaluate its displayed original derivative using E_Z/A_Z with all same-source factors and the genuine exponential. Acceptance: narrower actual p_Z enclosure and pressure memory preserved; never differentiate a C0 cap.
- [ ] **FIRSTBRIDGE-PRESSURE-PERIOD-INTEGRAL:** implement actual original complete-period phase/mean identities with true native radius and cutoff/source variation if cancellation is needed to reach useful pressure bounds. Acceptance: justified integral cancellation with retained errors and all source pieces, not local antisymmetry.
- [ ] **NEXT-REFINEMENT-TO-O2:** replay unchanged accepted O2 records after each real upstream improvement; quantify correction/own/absolute pressure C0/Z and preserve original family/datum/N/units. Rerun O2 only if its actual phase/source/candidate changes.
- [ ] **FULL-AXIAL-AND-ORIGINAL-RC-CLOSURE:** complete microscopic crossing-Z/whole axial atlas, original remaining O2/O3/pulse/end/flatten/heat cells, correlated Rc targets, independent five repair controls, global frequency admission and required higher mixed jets.
- [ ] **MATCHED-BACKGROUND-TO-CORRECTED-UVW:** finish five-moment/analytic-pressure/exact-heat matching, global admissible stress/flat remainder, actual n-dependent coefficient recursion, both oscillatory pulse families and independent Cartesian velocity/residual/dynamics. Full corrected Linfinity/volume-L2 residual<1e-3 remains the final target.

The broader [forward queue](CURRENT_UPSTREAM_TILE_UNIFORM_C1_O2_REPLAY_2026_10_08.md) remains active. Mark a task complete only with actual source/results/receipt/commit and its measured acceptance condition. Genuine scale recursion is still incomplete; narrower local source derivative enclosures do not establish it.
