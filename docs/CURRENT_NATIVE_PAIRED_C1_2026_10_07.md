# Exact paired Poisson Z derivatives and the current weighted target frontier

Checked source: [38500227](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/3850022715a6469a620ec50fc794284ba0c188db). Same original source family, candidate N2048, full Z[-1,1] and24 continuous radial cells. This extends [whole-period Z support](CURRENT_NATIVE_PERIODIC_C1_SUPPORT_2026_10_07.md), using the original q² mixed jet and C0 primitive supports without changing the defining field. All five terminal Z absolute upper bounds strictly improve again; no strict C0 reduction is claimed in this step.

## Exact original paired factors

Let h=sqrt(1+u²), r=u/h, s=h^-2=1-r² and D=1-2r*cos(psi)+r². The original period integrals give P=h^-1 W1 and H=h^-2 W2. For the original Mobius angle E,

```text
P+u*P_u = s^(3/2)*sin(psi)/D
2H+u*H_u = E + r*s*sin(psi)/D
            + s²*sin(psi)*(cos(psi)-r)/D²
```

The exact original identity D²-s² sin²(psi)=((1+r²)cos(psi)-2r)² implies |P+uP_u|<=sqrt(s)<=1. Also D=(cos(psi)-r)²+sin²(psi), so2|sin(psi)(cos(psi)-r)|<=D. With D>=(1-|r|)² and s²/D<=(1+|r|)²<=4, the second paired factor has absolute bound2pi+3. Both statements cover signed r and extend continuously through r=0. The original signed/small Fourier evaluator remains the defining function; these expressions are uniform outer-bound proofs, not a standalone singular r=0 evaluator.

For c=p2/dstar and u=c*q, with dstar constant in slow coordinates,

```text
(qP)_Z  = q_Z*(P+u*P_u) + q²*c_Z*P_u
(q²H)_Z = (q²_Z/2)*(2H+u*H_u) + q³*c_Z*H_u
c_Z = original p2_Z/dstar
```

The accepted whole-period estimates |P_u|<=20pi and |H_u|<=1000pi then give

```text
paired_qP_Z_upper  = QZ + 20pi*R*CZ
paired_q²H_Z_upper = (2pi+3)*RZ/2 + 1000pi*R*Q*CZ
```

Here Q,QZ,R,RZ,CZ bound the absolute values of the original q,q_Z,q²,q²_Z,c_Z. They use the same original context/factor basis/ledger and branch-local q² jets. No division by q or p2 is introduced, and no derivative/cutoff-crossing source is zeroed. At u near0, P+uP_u is generally O(1); the linear q_Z sensitivity is genuine and remains.

The original fixed-angle second-integral derivative becomes

```text
firstZ = 4pi*TZ*Q + T*paired_qP_Z_upper
T2Z_upper = 4pi*T*TZ + 4*firstZ + 4*paired_q²H_Z_upper
```

The B_Z fixed-angle product rule uses AU*paired_qP_Z_upper/(2pi), replacing the former separate AU*(2QZ+10Q*UZ). Its original moving-angle term AU*(T+1)*L/(4pi), a_Z contribution, A_Z derivative, E_Z product rule and source cross terms remain. The original nu_Z and implicit phase rules are unchanged. C0 support caps are not differentiated.

## Implementation and independent evidence

Files use `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_native_paired_C1_transport.py`: exact paired theorem, guarded three-site adapter of the original whole-period bounds, original branch-local density/oracle, complete route and applied-support weighted frontier.
- `.json`: actual24-cell transport, all ten terminal comparisons and current/preceding weighted contribution frontiers.
- `_check.py` / `_check.json`: original GenericShearLoop integral references, paired identities/bounds, original scalar phase/Z support comparisons, actual source traces, arithmetic/memory guards and contribution ranks.

The shared AST replacement helper rejects any count other than one at every specified support/query site; the firstZ insertion is also checked. Only equivalent derivative estimates change. Original C0 A/B, radius phase, E/V/cross terms, source units, true microscopic Duhamel widths, incoming history and separate P0/P0_Z remain. Nonlinear source evaluation precedes every overlapping branch hull.

Independent focused checks:3 exact symbolic identities;210 original P/H paired identity/bound/partial-derivative comparisons across signed, zero, small and large u;54 original scalar C0/y/Z/phi comparisons, including12 direct A_Z/B_Z support comparisons before selecting a tighter first-jet range. Scalar references use the original GenericShearLoop integrals/differentiation,120 dps and1e-85 allowance for paired tests. Synthetic fixtures are not current field data. Actual route checks cover69 paired support traces,138 derivative decisions,240 live integral rows and quiet pressure memory. Read-only review: **GPT-5.6 Luna / max**, source-faithful paired formulas and frontier ordering checked. Working/index dependency audit:1125 hashes.

## Actual improvement and remaining scale

Of138 source derivative decisions, 82 use tighter paired support and56 preserve a tighter original range. The preceding unpaired stage used80/58. All five terminal Z bounds strictly improve; no terminal C0/Z bound worsens. Producer 101.828s / checker 7.125s.

The following numbers are approximate logarithms of absolute upper bounds; ratios compare logarithms, not target values or actual errors:

| Rc target | Previous log upper | New log upper | Previous/new logarithm |
|---|---:|---:|---:|
| M_Z | 8.147214e+408906090034569676 | 7.41046e+408906090034569676 | 1.099421 |
| I_Z | 6.630788e+408906090034569676 | 5.894034e+408906090034569676 | 1.125 |
| S_Z | 8.147214e+408906090034569676 | 7.41046e+408906090034569676 | 1.099421 |
| Cp_Z | 6.83107e+408906090034569676 | 6.094316e+408906090034569676 | 1.120892 |
| D=(J-M)/mu_Z | 8.147214e+408906090034569676 | 7.41046e+408906090034569676 | 1.099421 |

Upper logarithms are about9–11% lower than the previous whole-period stage. Remaining bounds are still enormous and insufficient for useful repair contraction. These improvements neither choose actual controls nor establish compatible N, terminal closure or coefficient recursion.

## Refreshed actual source frontier

Each current cell contribution is transported to Rc with the exact product of downstream decay factors before applying the original target/first-Z normalization. The current mass is already included once in each contribution; no second mass/Jacobian is added. Source primitive diagnostics consume the C0/Z ranges actually applied before density, with exact-zero records handled separately. Pre-support raw A/B estimates no longer determine the reported primitive frontier.

| Target contribution | Current dominant original cell |
|---|---|
| M | transition |
| I | axial |
| S | axial |
| Cp | inner_reference |
| D=(J-M)/mu | axial |
| M_Z | active_first_bridge |
| I_Z | active_first_bridge |
| S_Z | active_first_bridge |
| Cp_Z | active_first_bridge |
| D=(J-M)/mu_Z | active_first_bridge |

All Z contributions remain dominated by active_first_bridge. C0 ranks have moved to transition/axial/inner_reference; repeating the old all-rows first-bridge C0 claim would be wrong. Per-cell hull ranks do not prove joint cancellation or total target values.

## Executable next tasks

- [x] **SOURCE3c1-paired-q/u:** exact original paired P/H derivative bounds before density, retaining genuine q_Z and c_Z; all24 continuous transport and all five strict terminal Z reductions.
- [x] **SOURCE3c1-range-loss-refresh:** current/previous24-cell weighted contribution ranks with actual applied primitive supports. Reuse the current warm route and recorded frontiers; avoid another coverage inventory.
- [ ] **CONTROL1a-defining-functions:** consume `current_generic_moment_repair_operator.py` and `current_native_Rc_parameter_targets.py` contracts. Provide the exact source-bound five bump matrix, Gram/quadratic function action, inverse action and complete first-Z Picard map. Explicitly distinguish defining functions from existing sidecar/log-norm bounds. Preserve `B_exact*h+N*r(N,Z)+Q_exact(h)/N=0`; do not replace source functions by cap endpoints.
- [ ] **SOURCE3d-all-N:** keep original N^-1/N^-2 coefficient functions before evaluating N2048, with the exact same source parameter, primitive phase and incoming history dependence. Build compatible log-N inequalities from those functions; fixed-N ranges cannot be rescaled into them.
- [ ] **SOURCE3c1-first-bridge:** separate genuine linear q_Z from remaining interval loss at active_first_bridge. Use original a/Delta/source derivative identities and exact conditional source correlations, or split its actual source coordinates while retaining the full original domain, selected sc and true widths. Never infer zero derivatives from a C0 cap or choose branch endpoints as field values.
- [ ] **SOURCE3c2-joint-source:** form `Dk*Mk*f_k-Ac*Dm*Mm*f_m` and its full Z derivative before independent key hulls. Keep distinct m rate1/k rate3/2 masses, Ac_Z, C_Z and division by mu*Ac². One common mass for the difference is invalid. This follows the actual axial-dominated C0 joint contribution.
- [ ] **CONTROL1b/2:** actual convergent five C1 controls at a proved compatible N, repaired original field and whole-Z terminal identities with support/pressure/heat joins. Source bounds and finite log-N certificates are not solved control values.
- [ ] **HIGH/OUTER/REC/WAVE/PHYS:** higher jets and joins, remaining global cones/pressure/heat/energy, genuine n=1/n>=2 recovery and independent moment repair, smooth summation/flat remainder, two pulse families and averaged quadratic stress cancellation, corrected Cartesian NS and measured material winding/scale recursion.

Gate: `current_original_paired_Poisson_Z_support_and_fixed_N_transport_executed`. Actual controls/global compatible N/terminal closure/coefficient recursion/oscillatory/corrected NS gates remain false. Keep the full long-term goal active.

## Concrete CONTROL1a handoff from the read-only source scan

`current_generic_moment_repair_operator.py` lines36–110 currently computes directed bump/Gram/matrix/inverse and quadratic **bounds**; lines113–172 derive sufficient C1/log-N conditions. These sidecars do not define coefficient functions. `current_native_Rc_parameter_targets.py` lines92–123 and165–232 evaluate scaled target covers and all-N bounds. Its N1024 inlet query is not a selected repair N.

The defining source graph already exists in `current_native_Rc_function_transport.py` lines212–319: `build()` returns the24-cell original history roots, original A, joint numerator, normalized targets and `N_scaled_targets`. Its `evaluate()` lines334–388 requires an injected original source/parameter/integral oracle and fails closed when unsupported. These exact roots are the control-map input; range endpoints cannot replace them.

- [ ] Own a new `current_native_Rc_functional_controls.py` plus focused checker/receipts. Consume `NativeRcFunctionTransport.build()` and an explicit original C1 source oracle. Expose `controls(Z_box,N)` with defining graph/action roots and C1 enclosures; distinguish definitions, finite Picard iterates and certified fixed-point controls.
- [ ] Define exact original bump integral roots, with `g_i(x)=beta((log(x)-c_i)/ell)/(ell*J0*x)`, original centers/width and disjoint supports. Use exact integrands, not `fresh_weights()` interval endpoints: `D_i=int((x^-mu-1)/mu*g_i dx)`, `HI_i=int(x^(1/2)*g_i dx)`, `HS_i=int(x^(-1/2-mu)*g_i dx)`, `HCp_i=int(x^(-3/2-mu)*g_i dx)`, `C_i=int(x^-1/2*g_i² dx)`, `E_i=int(g_i² dx)`, `P_i=int(x^-1*g_i² dx)`. Retain the source-dependent mu and stable divided difference.
- [ ] Build the exact matrix/action with target row order `(M,(J-M)/mu,I,S,Cp)` and control order `(axial0,axial2,swirl0,swirl1,swirl2)`. The five nonlinear rows are `0`, `(C0*a0*e0+C2*a2*e2)/mu`, `0`, `sum(E_i*(a_i²-e_i²/2))`, `sum(P_i*e_i²/2)`. Disjoint bump supports remove cross-bump integrals; do not erase same-bump products or the divided joint row.
- [ ] Implement `h_next=-B(mu)^-1*(d_scaled(N,Z)+Q(mu,h)/N)`, `h0=0`, where `d_scaled` is the exact `N_scaled_targets` root. Preserve the defining equation `B*h+N*r+Q/N=0`. Feeding r as if it were N*r is a normalization error. Bind an exact inverse action and its separate enclosure proof, not a point chosen from an inverse interval.
- [ ] Implement the derivative of the same map, or exact implicit action `(B+D_hQ(h)/N)*h_Z=-d_scaled_Z`. A_Z/A already belongs to the exact target graph. Maintain separate P0/P0_Z inputs for field recovery.
- [ ] Prove compatible finite N and a directed Picard tail before labelling controls solved. Existing fixed-N2048 bounds are too large to support contraction. Store the actual iterate/tail provenance, source family, namespace, domains and hashes.
- [ ] Independently replay the five **corrected physical moment integrals** and their Z derivatives on a fresh full source Z interval. Require both `B*h+N*r+Q/N` and its derivative to enclose zero, with a certified tail, then verify functional terminal identities and seams. A synthetic graph callback or point fixture cannot establish actual closure.
- [ ] Do not use the legacy `SharedFiveMomentRepair.coefficients()` from `lei_ren_part1_paper_shared_five_moment_repair.py`; it belongs to a different source family/defect map. Keep `actual_five_controls_installed` and `actual_terminal_Z_function_closure_installed` false until the current original functions and corrected field satisfy their scope.

Read-only source inventory model metadata: **GPT-5.6 Luna / max**. This section is an implementation handoff; it does not claim CONTROL1a complete.
