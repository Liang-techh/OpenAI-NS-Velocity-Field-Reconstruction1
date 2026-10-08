# Centered all-N control bridge and unit C1 ball

Checked source [57aec530](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/57aec530daf25bcbdd03c8e9df34c8af3eed33c2). Predecessor: [CURRENT_UPSTREAM_CENTERED_C1_2026_10_08.md](CURRENT_UPSTREAM_CENTERED_C1_2026_10_08.md). Reused **GPT-5.6 Luna / max** for read-only interface and Banach-estimate review; root implemented, computed, checked and published.

## Implemented behavior

The new typed **SavedCenteredWholeZ -> CenteredRcAllNControlsBridge** connects the checked centered whole-Z/all-N target bound to **NativeRcC1IntegralRealization** and the original convergent five-control family. It binds the same original owner, arithmetic context, source family, source graph bytes, 24-cell order and exact geometry. Existing strict NativePairedC1Transport interfaces remain unchanged.

This bridge uses the original centered certificate on **Z[-1,1], all integer N>=160**. The two N1024 strict-sign tile results from the predecessor are retained separately and are not promoted to a whole-domain theorem.

The exact function graph is identical to the accepted control family: **5220 nodes, 168 original coefficient/source references, 336 genuine C0/Z definite integrals**, one shared N in the actual spatial phase frac(N*log(R/r_minus)), source recipes, histories, targets and controls. Original kernels and native radius Jacobians remain applied once. The 35-branch original nonlinear density certificate is retained. P0/P0_Z remain separate from zero-inlet correction histories; quiet power correction is zero and pressure retains rate-zero memory.

Centered averaging bounds the **complete** target r by a coefficient/N^2. It does not set any exact order -1 coefficient function to zero. The original exact F_N=E*A*exprel(A/N), its Z derivative, V/V_Z and B/Pstar remain in the defining graph. No midpoint, range endpoint or logarithmic cap becomes an original function value.

Four new files use stem **lei_ren_part1_paper_compliant_current_centered_all_N_controls_bridge** under experiments/root_st073: producer, manifest, checker and receipt. Producer 0.468s; focused checker 0.531s. No source ancestors or inverse solvers are rebuilt.

## New control estimate

Write the exact original normalized equation as

    B(mu)*h + d(N,Z) + Q(mu,h)/N = 0,  d=N*r.

Read C from all five centered target C0 and Z coefficient covers **before** the N=160 substitution. Then ||d||C1<=C/N. In the same original normalized C1 control coordinates, CA bounds B^-1 and CQ bounds the full quadratic map.

For rho=1, sufficient conditions are:

    N >= 2*CA*C,
    N >= 4*CA*CQ,
    N >= 2*g_max*2^(2/3).

Together they give first Picard norm <=1/2, image radius <=3/4 and Lipschitz constant <=1/2. Therefore the exact source-local C1 limit is unique in the unit ball, with tail <=2^-n. Its stronger frequency-dependent bound is

    ||h_star||C1 <= min(1,2*CA*C/N).

The original h/N bump normalization further yields g_max*min(1,2*CA*C/N)/N for the normalized F/G C1 norm **in Z at fixed repair x**. Spatial x/y derivatives are not admitted by this estimate.

The contract includes all 17 original generic source frequency bounds and conservatively retains the accepted q-flat/active/quiet/repair-band C0 frequency budget. Its previous normalized control ball is checked to contain the unit ball. That local budget is not reoptimized here.

Let H0=10^408906090034569676. Rounded logarithmic quantities are:

| Quantity | Recorded logarithm / H0 |
| --- | ---: |
| New source + centered unit-ball repair sufficient log N | 9.624298732 |
| Retained source/local C0 sufficient log N | 48.16798785 |
| log C1 control upper bound at the retained local floor | -38.54368912 |
| log normalized F/G fixed-x Z-C1 upper bound at that floor | -86.71167698 |

These are formal source-bound inequalities, not evaluated defects, physical residuals, executable integer values or completion percentages. The old all-N repair threshold was 7.410459670*H0 for a much larger ball; this change **does not claim a lower sufficient frequency**. It proves a stronger control-size conclusion at the already retained local frequency. The enormous finite integer remains represented by its exact definition/log threshold.

## Verification and scope

The focused check independently proves the Banach/positivity/N-factor identities, compares the entire exact graph and nonlinear roots with the accepted family, checks every genuine kernel/coefficient/Jacobian integral and its common-N phase, and checks all five coefficient/floor bounds. It preserves the quiet cells and zero-rate pressure memories. **1176 Git-index dependency hashes PASS.** Accepted primitive, cutoff, IBP, source regularity and local cone proofs are reused.

The adapter supplies CachedC1ControlEvaluator with the new contract when an actual oracle is explicitly provided. An original source point/integral oracle is still missing. No numerical controls, finite terminal function closure, globally admitted common N, axis/whole-field higher regularity, exact heat match, global stress cone, coefficient recursion, oscillatory correction or corrected uvw are declared complete. All corresponding flags remain false.

## Executable next queue

This is the current queue. A task is DONE only after its concrete implementation and focused evidence are committed. Preserve the same source family and document conditional versus actual results.

- [x] **CENTERED-TO-ALL-N-TYPED-ADAPTER:** checked whole-Z averaged full-target sidecar bound to the same original C1 integral owner, graph, family, geometry and hashes.
- [x] **CENTERED-SOURCE/LOCAL-UNIT-BALL-CONTRACT:** original nonlinear equation, all 17 source bounds, original bump positivity, retained q-flat/band C0 budget, smaller C1 control and fixed-x profile bounds, exact factored tails.
- [x] **EXACT-GRAPH/INTEGRAL/N/P0-REGRESSION:** full graph unchanged, 336 genuine integrals, 168 common-N source refs, nonzero exact order -1 functions retained, quiet pressure memory unchanged.
- [ ] **ORIGINAL-SOURCE-POINT-ORACLE (next implementation):** install real same-source parameter/E/V/A/B and ordinary slow-Z callbacks for the 17 charts. Use the exact original cutoff and monotone implicit phase inverse; bind chosen N, source geometry and original P0/P0_Z. Cover flat/support/phase seams and tiny nonzero arguments. Range or log-cap values are forbidden as point data. Reuse the typed bridge's evaluator entrypoint.
- [ ] **ORIGINAL-24-CELL-INTEGRAL-ORACLE:** evaluate the exact current coefficient/target functions with the correct physical log-radius maps, kernels, one native Jacobian and same phase N. Keep original incoming histories and P0 distinct. Return error bounds for source evaluation, quadrature and roundoff independently of the Picard tail. Avoid direct enumeration of astronomically many periods; any averaging/asymptotic quadrature must include the actual mean and controlled remainder.
- [ ] **EXECUTABLE-FREQUENCY-STRATEGY:** distinguish the enormous proved sufficient threshold from a necessary frequency. Either derive source-correlated tighter sufficient conditions for a materializable N or evaluate a rigorously defined logarithmic/oscillatory integral representation. Never call an arbitrary convenient finite N globally admissible. Record which control, cone and source conditions that N satisfies.
- [ ] **LOCAL-C0-BUDGET-REOPTIMIZATION:** inspect how the old large control ball enters active/q-flat/quiet/band state-error polynomials. Substitute only the proved smaller same-coordinate control bounds and retain source terms. Recompute signed margins and frequency inequalities if this removes a meaningful obstruction; do not rerun unchanged whole-source producers.
- [ ] **GLOBAL-N-CONTRACT:** join source/repair/local C0 inequalities with the missing global cone, spatial join, higher total-y/phase-held derivative, exterior/energy and later recursion/pulse hierarchy conditions. Mark each condition proved, conditional or missing, and bind the identical N through every layer. Current local threshold is insufficient for global admission.
- [ ] **ACTUAL-FIVE-CONTROLS/TERMINAL-FUNCTIONS:** feed the actual oracle to arbitrary-depth cached Picard evaluation; map h_star/N to the five original compact bumps. Independently replay all five cumulative moment identities as Z functions, including derivative and pressure compatibility. A conditional C1 limit or interval linear response is not an installed numerical field.
- [ ] **SOURCE-LOCAL-SENSITIVITY:** tighten only the eta/a/source correlation or mixed-jet terms that block a useful frequency or oracle. Preserve genuine q_Z, p2_Z and cutoff sensitivity; precision alone cannot eliminate them.
- [ ] **AXIS/SPATIAL-JOINS/HIGHER-JETS:** establish whole-Z/axis/edge source atlases and actual spatial velocity, pressure and stress seam derivatives. The C1(Z) integral realization does not prove spatial seam equality.
- [ ] **ANALYTIC-P0/EXACT-HEAT/ENERGY:** connect the original datum to the preheat field and exact finite-energy heat exterior with legal collar and joins.
- [ ] **STRESS/FLAT-REMAINDER/TRUE-RECURSION:** certify residual=-div(T_B)+E_B, signed cone margins and flat remainder; implement distinct n=1/n>=2 equations and moment repair before controlled smooth summation.
- [ ] **TWO-PULSE-FAMILIES/CORRECTED-UVW:** build both oscillatory families and mean corrections against admitted stress; verify averaged quadratic cancellation, Cartesian NS/divergence/energy, contraction/elongation and actual material winding. Animation remains secondary.
