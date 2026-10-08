# Original reference whole-cell C0/Z five-density integrals

Checked source [b76d6ed6](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b76d6ed614b40a23b4ae7e97d8a1fadc587aeef0). Predecessor: [CURRENT_TWO_CHART_COEFFICIENT_RUNTIME_2026_10_08.md](CURRENT_TWO_CHART_COEFFICIENT_RUNTIME_2026_10_08.md). The existing GPT-5.6 Luna/max worker reviewed source formulas read-only; root implemented, computed, checked and published. No additional agent was spawned.

## What now executes

**OriginalSourceFactorAtlas** collects the original parameter identities before logarithmic evaluation:

```
logdelta = -4logPstar - 30
logR = log110 + 10logCstar + 10logPstar + coordinate
(p,d,ell,0,r) -> (p-4d+10r,10r,ell,0,0)
offset -> offset-30d+r*(log110+coordinate)
```

One fixed-Z atlas owns the directed context, basis and ledger for all eligible cells. It rebases issued Rh_reference/O2_slope values only when family, graph, selected parameter definitions, context and ledger agree. A nonzero extra source-unit slot is rejected. The atlas does not yet handle every chart's frozen amplitude units. Native R is never exponentiated, and astronomical radius logarithms are never numerically subtracted to recover a small coordinate.

Addition compares the **collected relative scale** first. Common astronomical logC factors otherwise hide the ordering of Pstar/delta terms and can destroy a useful signed enclosure. The original positive delta remains a source factor; its finite positive tail budget is also retained in the enclosure of L=1-delta*Z^2.

**OriginalReferenceWholeCells** evaluates the closed original E,V,a,b,p1,p2 and ordinary Z rows on entire exact radial cells. It reuses the actual original pressure datum and its directed alpha enclosure. Pressure-affine source terms keep the sharper all-late error as an explicit **Pstar^-1** factor before physical amplification. A point error floor is not substituted for that factor.

The actual candidate-N radius phase, including its positive original origin budget, is covered across each full cell. The original strict inverse maps a complete phase period onto a complete angle period; nonlinear density coefficients are evaluated on their source piece before union. The unchanged coefficient-pair program produces orders -1 and -2 for m,h,k,e,p and their ordinary Z derivatives.

Each coefficient enclosure is multiplied by the exact positive integral of its own decay kernel. Rates are m/e=1, h/k=3/2 and p=0, with native log-radius Jacobian1. The pressure rate-zero mass is the cell length; the total reference-window mass is5. This computes a **window contribution** only: it neither sets incoming global histories to zero nor resets P0.

Five implementation/evidence files under experiments/root_st073 are committed: the atlas, the whole-cell producer, its lossless compressed report, checker and receipt. The stem for the latter four is **lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals**.

## Implicit derivative improvement

Naively multiplying independent ranges for t and psi_Z can create an astronomical B_Z cover at the original narrow peak, even when the functional product is bounded. The new signed E-chart path retains the original fixed-angle T1,T1_Z,T2_Z and uses

```
psi_Z = -T2_Z/(1+t^2)
t*psi_Z = -T2_Z*t/(1+t^2)
-1/2 <= t/(1+t^2) <= 1/2
B_Z = -a*(E_Z*T1+E*(T1_Z+t*psi_Z))/(4*pi)
```

The bounds follow from (t-1)^2>=0 and (t+1)^2>=0. This encloses an actual rational function; no endpoint, midpoint or saved field cap is selected as a value. The older accepted slow-Z module is unchanged. This local improvement requires strict signed geometry and the original a_Z=b_Z=t0=0 identities.

## Computed contribution ranges

The complete **Rh_reference[-5,0]** source window was computed at exact Z=37/100 and common candidate N=160 with4 and16 cells. The following rounded upper bounds contain the recorded absolute contribution envelopes:

| Row | 4 cells | 16 cells |
| --- | ---: | ---: |
| m ordinary Z |0.116682|0.003682|
| h C0 |0.000999|0.000971|
| h ordinary Z |0.014441|0.003206|
| k ordinary Z |0.074288|0.002068|
| e C0 |0.001235|0.001154|
| e ordinary Z |0.018174|0.004531|
| p C0 |0.004814|0.004394|
| p ordinary Z |0.059375|0.016632|

m C0 and k C0 retain their microscopic original source factors in the lossless report. The m/k Z rows retain L^-1, which is included in the finite ranges above. These are normalized source-window contribution bounds, **not** five terminal defects or a Cartesian momentum residual. Both signs are generally allowed. They do not establish an accurate signed integral or globally sufficient N.

Spatial refinement tightens the Z envelopes substantially, but whole-period range integration still loses oscillatory cancellation. The h C0 bound changes relatively little. More of the same radial subdivision alone is therefore not the main next step.

Producer **23.375s**, focused checker **50.844s**. Evidence:

- 5 exact affine scale identities,5 independent finite-unit diagnostics and an actual huge shared-factor/delta ordering regression;
- 2 exact rational-factor bounds and16 comparisons with independent defining T1/T2 integrals and their Z derivatives, including both signed geometries;
- 16 issued original two-chart rebases and20 actual point coefficient enclosures contained in a whole reference cell;
- 7 frame/Z/ledger/parameter/midplane/partition rejections;
- 10 independent positive own-rate mass comparisons,2 pressure C0/Z identities and19 exact pressure-affine template rows;
- **1200 Git-index dependency hashes PASS**.

Finite-unit fixtures are diagnostic evidence only. They never become original construction parameters, source values or frequency choices. No ancestor producer was rebuilt.

## Executable queue

Read this page before older handoffs. DONE means committed implementation plus a scoped receipt. Update the task checkbox and link its commit when completing it. Preserve the current source family, pressure/error data and unchanged dependency receipts. Give editing workers explicit ownership; keep our repository-scanning Luna worker read-only.

- [x] **ORIGINAL-MACRO-FACTOR-ATLAS/FIXED-Z:** identical source parameters, exact affine radius/delta collection, shared basis/ledger and issued-frame checks for eligible Rh/O2 values. Extra frozen source units remain outside this DONE scope.
- [x] **REFERENCE-WHOLE-RADIAL-CELL-C0/Z-FIXED-NONZERO-Z:** entire closed radial cells, original pressure-affine late error factors and actual phase covers.
- [x] **SIGNED-IMPLICIT-Z-PRODUCT:** bound t*psi_Z before multiplication; preserve the original primitive derivative formulas and source factors.
- [x] **REFERENCE-OWN-RATE-WINDOW-CONTRIBUTIONS:** all five C0/Z densities, orders -1/-2,4/16 exact cells and zero-rate pressure memory. Incoming histories and P0 stay separate.
- [x] **SCOPED-INTEGRAL-EVIDENCE/PUBLICATION:** independent defining-integral, point/whole-cell, memory, factor and index evidence.
- [ ] **REFERENCE-WHOLE-Z-C1-CELLS:** implement source-owned Z intervals across[-1,1], preserving pressure parity and the original L dependence. Split strictly signed p2 cells. Do not extend fixed-Z receipts by sampling or interpolate between them. Acceptance: true C0 and ordinary Z function covers over the complete window, with explicit subdivision/domain proof.
- [ ] **MIDPLANE/SIGNED-TRANSITION-INTEGRATOR:** implement the actual p2=0/nonzero-p2_Z branch and small-u transition around Z=0. Preserve nonzero source derivatives and positive rho/s/hinv. Acceptance: compatible signed/small-u transition enclosures and an actual midplane C1 integral; Z=0 is currently rejected.
- [ ] **PHASE-AWARE-REFERENCE-INTEGRATION:** replace complete-period range multiplication with original phase means plus a controlled remainder for radial coefficient drift and the decay kernel. Retain actual N in both inverse phase and exp(A/N). Do not set order -1 functions or their means to zero without an exact source identity. Acceptance: directed signed integral evidence and a bound separating oscillatory remainder from source/inverse/rounding error, without enumerating an astronomical period count.
- [ ] **ACTUAL-O2-C1-WINDOW-INTEGRATION:** connect accepted O2 defining quadratures, ordered source cells and positive log-q kernels to the same atlas and actual all-N coefficient program. Include ordinary Z, whole-Z sign transitions and separate P0/incoming history. Acceptance: one source family, shared candidate N and complete O2-window contribution exports; prior C0/cell receipts are inputs, not full C1 completion.
- [ ] **O2-AXIAL/11-UNIT-BUFFER-SOURCES:** expose genuine axial and buffer function roles and inherited moment histories. Give each native chart its exact radius offset and Jacobian. Acceptance: issued live source bindings and C0/Z whole-window integrals with seam/incoming-history identities.
- [ ] **ATLAS-REMAINING-SOURCE-UNITS:** extend the typed atlas for original positive-mu and frozen amplitude units needed by O3/reshape/restore. Preserve exact definitions and compatible frames. Acceptance: symbolic factor conversions and rejection of incompatible units; never discard a nonzero fourth source slot or export an opaque giant log as though it preserved correlation.
- [ ] **O3-POSITIVE-MU/POWER/RESTORE-INTEGRALS:** implement remaining source-owned roles from genuine source programs, followed by inner/bridge/patch/restore routes. Acceptance: actual original coefficients and positive parameters, C0/Z integrals and inherited histories. Saved field covers, interval midpoints and shear-only replacements are prohibited as source values.
- [ ] **TYPED-FACTORED-INTEGRAL-PROTOCOL:** add a separately typed partial integral evaluator with owned chart/window, family/graph/hash, exact coordinate/Jacobian, Z domain, candidate N and source/inverse/rounding error records. Acceptance: unsupported roles and charts fail explicitly. Keep the existing scalar FunctionEvaluator guards unchanged.
- [ ] **ALL-17-CHART/24-CELL-ORACLE:** bind complete original C0/Z integrals, pressure and incoming histories to every required control row. Acceptance: actual source/integral callbacks across the complete route; no cap values, zeroed incoming corrections or hidden missing charts.
- [ ] **EXECUTABLE-FREQUENCY/ACTUAL-FIVE-CONTROLS:** distinguish diagnostic N160 from a globally admitted frequency. Use source-uniform phase estimates at immense frequencies. Bind the actual integral oracle to the accepted centered whole-Z/all-N/unit-C1-ball construction and cached Picard evaluation. Acceptance: actual five controls and all terminal identities as Z functions, with numerical oracle errors separate from the Picard tail.
- [ ] **BACKGROUND-JOINS/EXACT-HEAT/STRESS/FLAT-REMAINDER:** retain the previous full queue and verify these conditions on the same completed background and pressure datum. Acceptance: functional joins and actual admissible-cone/flat-remainder margins, not only geometry or local residual.
- [ ] **N-DEPENDENT-COEFFICIENT-RECURSION:** implement distinct n=1 and n>=2 recovery equations, order-specific moment repair and smooth summation. Acceptance: actual recursive coefficients and finite-order remainder diagnostics. Rescaled copies of one field do not satisfy this task.
- [ ] **TWO-PULSE/CORRECTED-CARTESIAN-UVW:** implement averaged quadratic stress cancellation and the final corrected velocity/forcing, then independent Cartesian residual and contraction/elongation/winding diagnostics. This remains the long-term goal.

All17/24 oracle, actual five controls, selected global N and recursive corrected-field flags remain **false**. The current fixed nonzero-Z contribution is one implemented integration layer toward that goal.
