# Original all-N coefficient functions and uniform C1 repair conditions

Checked source: [b4ed26be](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/b4ed26be1ced3e86f30c6a797e6b40f9f9107309). This follows the [exact five-control map and finite-iterate stage](CURRENT_NATIVE_FUNCTIONAL_CONTROLS_2026_10_07.md). The latest source graph contains 4,858 nodes. Focused checker: PASS, 18.594s. Working/index dependency audit: 1,135 matching hashes. Read-only reviewer metadata: **gpt-5.6-luna / max**. No new worker or source ancestor constructor was needed.

The five-moment repair now has exact **N-dependent coefficient functions**, continuous all-N C0/Z histories and fresh uniform C1 target bounds over the complete original domain. All five transformed C1 bound upper limits strictly improve against the preceding original all-N target baseline. This is a reduction in conservative bounds, not a measured improvement in the actual corrected field or a percentage of project completion.

The new bounds feed the repair ball/contraction/positivity conditions. The recorded lower threshold also includes all 17 original generic source-frequency requirements. A conditional C1 Picard tail follows for the same source family after the stated conditions hold. **A globally compatible N, instantiated fixed-point controls/tail, corrected terminal closure and coefficient recursion remain open.**

## Exact coefficient functions and original phase

Every source reference retains one shared positive integer N and the actual phase `fractional_part(N * original_log_radius_minus_inlet)`, the original native coordinate and Z variable. Write `mathcal_A, mathcal_B` for the original periodic primitive pair, and `E,V` for the original axial/swirl background inputs. These are function handles from the accepted original graph, including their first Z derivatives; they are not selected values from a cover.

```text
F_N = E * mathcal_A * exprel(mathcal_A/N)
F_N_Z = E_Z * mathcal_A * exprel(mathcal_A/N)
        + E * mathcal_A_Z * exp(mathcal_A/N)

                N^-1 coefficient       N^-2 coefficient
m               mathcal_B              0
h               F_N                    0
k               V*F_N + E*mathcal_B     F_N*mathcal_B
e               2*V*mathcal_B - E*F_N   mathcal_B^2 - F_N^2/2
p               E*F_N                  F_N^2/2
```

First Z derivatives use the complete product rule, with nonzero original V/V_Z retained. The formula is continuous at `mathcal_A=0`, including nonzero `mathcal_A_Z`. Although the powers N^-1 and N^-2 are extracted, their coefficients still depend on N through exprel and the spatial phase. They are not an N-independent polynomial expansion.

Each coefficient is integrated through its own original cumulative kernel (rates 1, 3/2 or 0) and the true `dy/dcoordinate` once. The two coefficient histories start from the genuine original zero inlet and pass through all 24 continuous radial cells. Quiet pressure memories remain; histories do not reset at a cell join. Exact original unsplit density histories are retained separately for comparison.

Normalization uses the same original positive amplitude A/A_Z and strictly positive mu. In particular, the divided row is formed from the joint numerator `k-A*m` and its complete Z derivative before division by `mu*A^2`. The exact targets reconstruct as

```text
N*r(N,Z) = coefficient_target[-1](N,Z)
           + coefficient_target[-2](N,Z)/N.
```

The five-control Picard map, residual, implicit Jacobian and band recipes now take this exact reconstructed function family. A support cap, midpoint or previously saved N2048 range never defines a coefficient function.

## Fresh uniform ranges and frequency conditions

The range route queries the original q/q^2 and paired primitive derivative supports over the whole fractional period, full `Z in[-1,1]` and original continuous native coordinates. It retains 35 actual conditional cutoff branches. Nonlinear coefficient ranges are formed inside each original branch before overlapping branches are united. Original factor bases, ledgers, distinct true masses and incoming memories are carried through all 24 cells.

The local factor cover is valid for integer N>=160: `|mathcal_A|<=5/4` gives `|mathcal_A/N|<=1/128`. A positive exp outer cover encloses both exp and exprel factors. This local safety condition is not a global frequency selection. Raw total-y fast-phase derivative terms of order N^0 remain in the original higher-jet sidecars; a uniform N^-1 bound for all y derivatives is not claimed.

Fresh weights and positive inverse bounds at the actual source mu turn the new five target C1 caps into repair ball, half-contraction and positivity lower conditions. The combined directed log-N threshold is the maximum of log(160), this repair/positivity threshold and the 17 original generic source thresholds. It records the finite-integer recipe `N>=ceil(exp(the directed threshold))`; the huge exponential/integer is not materialized.

These admitted requirements are lower bounds. Remaining global cone, higher-jet and function-join N constraints have not been assembled or checked for compatibility. The record therefore keeps `current_whole_N_selected=false`.

For the exact same-source C1 family, after the ball conditions, `L<=1/2` and `||h1-h0||_C1<=rho/2`, the standard tail

```text
||h-h_n||_C1 <= L^n/(1-L) * ||h1-h0||_C1 <= rho*2^-n
```

holds conditionally. This establishes the tail algebra and a family contract. It is not yet an instantiated numerical tail for a selected globally compatible N or an evaluated original fixed point.

## Evidence, APIs and scope

- 210 exact comparisons with the original unsplit full signed density DAG, C0 and Z, over 21 nonquiet cells.
- 168 genuine N-dependent C1 coefficient integral functions; 40 exact history/target C1 and N reconstruction checks.
- 90 independent manufactured original increment C0/Z scalar comparisons at N=160,257,2048, including nonzero V/V_Z and zero primitive with nonzero derivative; the same 90 references lie inside the fresh uniform ranges. These references are not evaluations of the current original field.
- 168 actual-phase coefficient source references; 336 nonzero signed C0/Z integral rows and 144 exact quiet/zero rows, over all 24 original cells.
- 35 original conditional branches, eight quiet pressure coefficient C0/Z memory rows and all 17 original generic source requirements.
- 24 exact Picard/residual/band C1 identities and 25 independent implicit Jacobian entries; conditional uniform tail algebra.

Files use prefix `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_native_Rc_all_N_function_controls.py`: exact phase-bound coefficient functions, original all-N range transport, combined source/repair conditions and artifact recording.
- `.json`: defining graph, all-N cell/branch histories and targets, improved C1 caps, directed log-N conditions and explicit incomplete gates.
- `_check.py` / `_check.json`: independent exact original density comparisons, scalar/range references, source provenance, C1 identities and conditional tail evidence.

Reuse `NativeRcAllNFunctionControls(controls_owner, paired_owner)` with the accepted live owners. `.build(iterations=3)` produces the exact graph; `.route(Z=(-1,1))` computes fresh all-N ranges. `record(owner,built,live,warm_original_route_reused=True)` can publish an already-computed route. The current artifact used this last path: `execution_seconds=null`, so no full producer timing is claimed. The route is complete and checked; do not rerun it merely to obtain a benchmark.

Original P0/P0_Z preservation is inherited from the accepted source; this turn did not replay pressure reconstruction. Interval exprel uses an outer enclosure, not a newly certified point residual oracle. Actual original point/full corrected physical integrals, solved controls, an instantiated fixed-point tail, terminal Z-function closure, global stress cone, genuine coefficient recursion, oscillatory correction and full corrected NS gates remain false.

## Detailed executable continuation

- [x] **SOURCE3d-functions:** bind both coefficient orders to original E/V/periodic primitive functions, their C1 derivatives and actual shared-N spatial phase; preserve exact original unsplit density roots.
- [x] **SOURCE3d-transport:** transport coefficient histories through the full original 24-cell domain with true kernels/Jacobians, zero inlet and quiet pressure memory; reconstruct N*r and the exact five-control graph.
- [x] **SOURCE3d-uniform:** produce original full-period/whole-Z uniform C1 ranges before selecting N; record all five strict bound reductions and all branch/mass provenance.
- [x] **CONTROL1b-lower-contract:** combine the improved repair/positivity conditions with all 17 admitted original source lower thresholds; state the conditional rho*2^-n tail without setting global-N or solved-control gates.
- [ ] **CONTROL1b-global-inventory:** read original global cone/higher-jet/heat/join conditions and classify each N dependence. Keep raw N^0 fast-y terms. Record lower, upper, equality and domain constraints, the exact source family and whether they have an accepted proof. Do not replace the inventory by a maximum of known lower bounds.
- [ ] **CONTROL1b-common-N:** assemble a typed same-family parameter definition meeting every required condition on the complete original domain. Check actual compatibility, including any upper/mixed constraints, before admitting one finite integer N. Store directed inequalities and retain the original positive mu/A/P0 provenance.
- [ ] **SOURCE3c1-first-bridge:** reduce the dominant active_first_bridge Z loss with genuine linear q_Z and conditional a/Delta correlations or source-coordinate subdivisions over the complete original domain. Keep actual microscopic widths, original sc and derivative coverage at branch crossings. Compare the active bound frontier after this production change.
- [ ] **SOURCE3c2-joint:** preserve complete signed m/k source correlations through the distinct rate-1/rate-3/2 kernels, incoming histories and A_Z, before the divided target is enclosed. Prove any cancellation from common functions; independent cap subtraction or a shared artificial kernel cannot establish it.
- [ ] **CONTROL1b-point-or-averaged-oracle:** implement a directed original point or rigorously averaged integral oracle for the same N-dependent graph. Resolve actual spatial phase without aliasing; give explicit quadrature/remainder errors and positive-denominator guards. Missing source/parameter/integral oracles must fail closed. Never use saved bounds as point coefficients.
- [ ] **CONTROL1b-instantiated-tail:** at the compatible N, obtain actual same-source iterate increments and a useful certified C1 tail. Retain the exact function roots alongside enclosures, propagate the tail through h_Z and band recipes, and choose depth from the requested residual tolerance instead of an arbitrary iterate count.
- [ ] **CONTROL2-band-field:** install convergent axial/swirl coefficient functions on the original power band, including h_Z. Recover radial velocity from the original divergence constraint and analytic pressure with separate original P0/P0_Z; carry inlet, axis and support conditions.
- [ ] **CONTROL2-five-moments:** independently integrate corrected original physical cumulative moments and Z derivatives on a fresh complete Z cover, including the certified tail and all five original weights. Require terminal identities and quantitative residual errors. A range containing zero or manufactured coefficients is insufficient.
- [ ] **CONTROL2-joins:** establish actual corrected inner/annulus/support/terminal pressure and heat function joins at the required derivative orders; verify incoming histories are carried rather than restarted.
- [ ] **HIGH:** build the required higher y/Z/parameter derivatives of the same corrected functions, retaining phase N powers and exact cancellations. Prove axis regularity and all function seams needed by the global stress and exterior steps.
- [ ] **OUTER:** assemble the remaining whole-domain admissible stress cone, background residual decomposition, analytic pressure/heat exterior and finite-energy tail for the corrected field, with a single compatible parameter family.
- [ ] **REC-n1:** recover the genuine first coefficient using the original shear loop and cumulative moment equations; perform its own independent moment repair on a common inner interval.
- [ ] **REC-n2plus:** implement the actual n>=2 recurrence, per-order moments/repairs, divergence-preserving cutoffs and finite-order residual estimates. Reusing or rescaling the present background field does not close these equations.
- [ ] **REC-flat:** choose summation/cutoff scales from the proved coefficient bounds and construct the smooth flat remainder with required derivative estimates.
- [ ] **WAVE:** install the two original oscillatory pulse families and verify averaged quadratic stress cancellation and remainder bounds against the completed stress cone.
- [ ] **PHYS:** produce the corrected Cartesian [u,v,w](x,y,z,t), independent NS/divergence/forcing residuals and measured vortex-core contraction, relative axial thinning and material winding as t approaches T. Animation is secondary to the field and diagnostics.
- [ ] **HANDOFF:** mark a task complete only with its exact source commit, same-family receipt and actual accepted gate. Preserve previous stages as history, update these four leading handoffs after substantive changes and push the checked work.

Gate: `current_original_N_dependent_coefficient_functions_and_uniform_C1_repair_contract_bound`. The full long-term goal remains active and incomplete.
