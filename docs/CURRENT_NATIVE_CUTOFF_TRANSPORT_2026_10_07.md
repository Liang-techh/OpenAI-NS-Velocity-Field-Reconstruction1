# Original cutoff-local source and fixed-N 24-cell C1 range transport

Checked cutoff/source implementation: [102c8fc0](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/102c8fc071de9e3babb5c891eee6fc36523d151a). Checked final transport implementation: [5aa2f83c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/5aa2f83c4ef2412175deb867f5680cacc35ece21). This supersedes the three-provider obstruction in [CURRENT_NATIVE_SIGNED_U_COVER_2026_10_07.md](CURRENT_NATIVE_SIGNED_U_COVER_2026_10_07.md).

All seventeen declared original full-Z interior providers now resolve at the same N=2048. More substantially, all twenty-four original continuous radial cells, with full Z[-1,1], produce fixed-N signed density and integral ranges. The five C0/Z correction histories are transported from the exact inlet through the full route to Rc. Five normalized Rc target C0/Z ranges are installed. These are conservative enclosures of the original source-function integrals; their very large widths do not establish small defects, convergent controls or terminal closure. Genuine coefficient recursion remains unfinished.

## APIs and source files

Files use `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_native_cutoff_q_cover.py` / `.json` / `_check.py` / `_check.json`.
- `current_native_cutoff_density_oracle.py` / `.json` / `_check.py` / `_check.json`.
- `current_native_cutoff_range_transport.py` / `.json` / `_check.py` / `_check.json`.

`NativeCutoffQCover(existing_NativeQSlowJets)` restricts the original C0 cutoff geometry into Delta<=0, 0<=Delta<=eta and Delta>=eta, where Delta=kappa-2 and kappa=a+b²/a. The negative and transition branches retain the original ordinary a/Delta derivative rows through y2/Z1. The transition uses the original smooth sigma and factorial-correct ordinary derivatives. Only the proven flat branch sets q and all q derivatives exactly zero. Original positive gamma/root lower proofs remain. Every possible branch is retained; shared endpoints are harmless overlaps.

The aggregate q query retains unconditioned packet/root metadata and a hull of q jets. It explicitly requires nonlinear phase/density consumers to use the conditional branch query objects. Each branch source q is the same object as its ordinary ZERO jet. The source-family, graph, context and original ledger are preserved.

The decisive primitive identity is

`A = a/2*(Phi-psi_fraction)`.

At the original inverse, Phi=phi and both phase fractions lie in[0,1]. When q!=0, Delta<eta and a<=kappa=2+Delta<2+eta<=5/2. When q=0, A=0. Therefore |A|<=5/4 on the original source and |A/N|<=1/128 for every integer N>=160. Both original small and Mobius formulas are checked symbolically; the latter uses s=1-r².

`supported_expm1(A/N,N)` retains the original signed formal A/N factor and bounds only its positive exponential average. `supported_exponential(A/N,N)` bounds only the exponential factor. The original A_Z, E/E_Z, V/V_Z and every signed density cross term remain. A sign-definite input range that contradicts the support theorem fails explicitly. No capped scalar defines A or a velocity field.

`NativeCutoffDensityCover(existing_NativeDensityC1LocalIntegrals)` runs the accepted signed-u phase/density body separately for every original cutoff branch. Its explicit AST adapters replace exactly two source calls and exactly two exponential-factor calls, each once. The nonlinear density is evaluated before either cutoff or signed-u hull. Phase-cell unions follow; overlapping branch integrals are never added.

`NativeCutoffFactoredOracle(role_owner,built=None)` retains the existing original graph-role/namespace/root/hash dispatch and separate P0/P0_Z and Rc E/E_Z. Its density frame is `{'C0': {m,h,k,e,p}, 'Z': {m,h,k,e,p}}` at the supplied candidate N. The existing scalar `transport.evaluate()` still rejects range frames.

`NativeCutoffRangeTransport(role_owner)` consumes the checked density receipt. `cell(label,chart,left,right,Z=...,N=...)` uses `NativeRcParameterTargets.geometry()` and `true_width_kernel()` for original coordinates, special endpoints and microscopic widths. `route(Z=(-1,1),N=2048)` propagates all24 cells. `run(role_owner,return_live=True)` returns the report, owner and live route output for focused acceptance; the normal return is the report.

The source ranges are in d/d(log R) units. Each positive true-width Duhamel mass is applied once after branch/phase hulls. There is no additional native-coordinate Jacobian. Microscopic widths are collected before absolute endpoint subtraction and stay nonzero formal factors. All source ranges rebase into the accepted common arithmetic coordinates; this does not restore lost correlations.

For each rate lambda_j,

`H_out = exp(-lambda_j*w)*H_in + mass_j*source_j`,

with the same equation for the ordinary first Z derivative. Rates are m=1, h=3/2, k=3/2, e=1 and p=0. Endpoints, weights and radius phase are Z independent. The exact initial correction is zero by the original inlet theorem; background histories and P0/P0_Z remain. Original inlet theorem flags are recorded explicitly. O3_power requires the original exact flat cutoff branch and identical q/ZERO source, but still attenuates every incoming nonzero-rate correction. Pressure correction memory is unchanged across quiet cells. Any unresolved source cell blocks all downstream history/target admission.

`fixed_N_target_rows()` accepts already evaluated fixed-N histories. It forms C=k-Ac*m and C_Z=k_Z-Ac_Z*m-Ac*m_Z before division by mu*Ac², and uses the full quotient derivative. Ac and Ac_Z come from the same original Rc amplitude; mu is checked against the original repair mu. This API is deliberately separate from all-N `target_rows()` coefficient maps indexed by powers -1/-2. No N2048 output is promoted to an all-N coefficient.

## Actual results and limits

- Original declared interior queries:17 enclosed /0 unresolved, N2048 and full Z[-1,1];160 density-role dispatches.
- Former blockers bridge_first and O2_axial now use conditional q slow jets. actual_patch tail exponential factors consume the original support bound at the same N2048. The generic selected sc/2 collar range also resolves, without claiming it is an exact-flat zero query.
- Full exact route:24 enclosed continuous cells across17 original charts;23 exact function-radius joins. This extends interior provider probes to source-cell integral enclosures.
- Five original C0/Z histories reach Rc; ten normalized target ranges are present. Positive kernel masses:120; live common-basis contribution ranges:240; contribution/incoming/outgoing C0/Z rows:720. Six quiet pressure C0/Z memory rows remain identical.
- Focused references:42 independent q ordinary derivative comparisons,30 independent original density C0/Z comparisons,30 independent normalized Rc target/derivative comparisons. Additional actual broad O2 source/integral and original Rc amplitude dispatch, and a continuous O2 cell at N160. All212 inherited namespace/root/role bindings and nine inherited rejection guards remain.
- Working/index dependency hashes:1103 for cutoff density acceptance and1107 for transport acceptance. Read-only scoped mathematical review: **GPT-5.6 Luna / max**.

The normalized target upper ranges are extraordinarily broad in logarithmic scale. This is a limitation of independent range propagation and source dependency loss, not a measured enormous physical NS residual. In particular, forming k-Ac*m from independently transported ranges cannot prove its quantitative cancellation. The current ranges establish computability and source coverage. They do not establish useful repair contraction or a globally compatible N. N2048 is a candidate for these range queries.

## Detailed next-agent tasks

Read this document and the leading section of CURRENT_CHECKPOINT.md first. Preserve the checked original source artifacts; add a separate layer. Keep completion checkboxes scoped and attach source commit plus producer/checker evidence when closing a task.

- [x] **SOURCE2e1:** original cutoff-conditioned ordinary q jets, active/transition/flat branch cover, same source/ledger and original derivative rows.
- [x] **SOURCE2e3:** exact original primitive/support identity and supported exponential factors; actual_patch resolves at the same N2048 without changing A/A_Z.
- [x] **SOURCE2d-all17-declared:** all seventeen original declared full-Z provider queries enclosed with explicit original role dispatch. This checkbox refers to the declared query set.
- [x] **SOURCE2e4-fixedN:** original continuous radial/full-Z source ranges on all24 exact route cells; every required cutoff/u/phase branch retained, original seams bound.
- [x] **SOURCE1b2-range:** full24-cell fixed-N C1 integral range transport, original exact inlet, microscopic true widths, incoming/quiet pressure memory and separate original P0.
- [x] **SOURCE3b-range:** five fixed-N Rc target C0/Z ranges and joint numerator algebra, original Rc amplitude/mu. Quantitative joint cancellation is a separate open task.

- [ ] **SOURCE3c1 dominant range-loss inventory:** retain live shared source recipes for every24-cell contribution and its subsequent decay. Report the dominant contribution to each normalized target before and after common-basis rebase. Separate native root dependency loss, cutoff/phase hull loss, mass/decay width loss, large Ac/Ac_Z division, and joint k-Ac*m subtraction. Do not label a broad range as an actual residual. Deliver a small machine-readable frontier and identify the smallest production change that improves a target bound.
- [ ] **SOURCE3c2 correlated joint target:** build the original shared-expression numerator at source/integral level, with correct unequal m/k Duhamel weights and the same Ac/Ac_Z. Collect common formal powers and near-equal exponential factors before independent range conversion. Preserve C_Z=k_Z-Ac_Z*m-Ac*m_Z and all original E/V cross terms. Establish a quantitative range for the divided row /mu with directed remainder evidence. Do not impose cancellation between unrelated interval caps.
- [ ] **SOURCE3c3 adaptive range transport:** partition only the dominant original radial/Z cells and required phase/cutoff crossings. Use exact native endpoints and true subcell widths, preserve the original common global phase and N, and compose all incoming memories in order. Compare the refined ranges with the current full-cell enclosure and report enclosure quality/cost. No point samples may define a source function or integral; do not refine all24 cells by default.
- [ ] **SOURCE3c4 useful target tolerance:** define a target-range width sufficient for the original repair contraction inequalities. Keep absolute target range, enclosure uncertainty and all-N sidecar cap distinct. Report whether a selected candidate N and source refinement meet the required bound. Do not declare useful controls from finite but enormous ranges.
- [ ] **SOURCE1b1 general range graph evaluator:** add a companion API to scalar evaluate for typed original range functions and exact parameter/radius/Jacobian nodes. Dispatch by graph namespace/root/role, family/hash and live context. Use directed formal exponential tails and collected radius differences. Preserve both C1 paths and reject unknown/unresolved sources. The specialized24-cell adapter already works; do not rebuild it merely to check the same stage.
- [ ] **SOURCE2e2 exact collar provider (optional separate layer):** bind the existing checked exact-flat initial-collar theorem to its admitted selected-sc domain and source hash. The generic sc/2 range now resolves but is not declared exact flat. Keep original E/V/P0 and incoming histories; never blanket-zero bridge_first or a missing gap.
- [ ] **SOURCE3d all-N correlated coefficients:** use original uniform_density_orders and separate N^-1/N^-2 density coefficient functions through the route. Retain their true N dependence, source phase and derivative rules. Apply true-width mass to each order; do not rescale fixed-N histories into coefficients. Preserve the shared joint numerator and show bounds sufficient to select a finite compatible N without materializing astronomical integers.
- [ ] **CONTROL1a exact repair graph (parallel):** extend the typed source-function graph with original five bump functions, exact B(mu) entries, Gram/quadratic Q(mu,h), inverse action and full C1 Picard map h_next=-B_exact^-1*(N*r+Q_exact/N). Keep sidecar inverse/weight bounds separate from defining matrix entries. Preserve the divided /mu row and all nonlinear cross terms. Deliver a function-map API and focused independent moderate-source reference, not midpoint matrices.
- [ ] **CONTROL1b controls:** consume actual typed original r/r_Z and the same original Ac,Ac_Z,mu,N. Solve the original nonlinear five-control system with contraction/error bounds. Compute the full Z derivative of its implicit equation, including parameter/source derivatives. h=-r and B=I are invalid substitutes.
- [ ] **CONTROL2 terminal closure:** evaluate all five original moment terminal functions over whole Z using corrected fields and independently bounded control error. Carry the repair through Rc->2Rc; keep post-repair2Rc->Rb admission separate. Preserve exact divergence structure, pressure datum and earlier incoming histories.
- [ ] **HIGH/OUTER:** extend required ordinary higher mixed jets, axis/seam smoothness, pressure/exact heat matching and full upstream stress cones on the selected corrected source. Existing scope-limited cone receipts do not prove the whole corrected tensor.
- [ ] **REC1:** implement the genuine n=1 recovery equations with original domain, inner datum and independent five-moment repair. Produce a coefficient field and finite-order remainder evidence; do not relabel coordinate scaling or fixed-N range transport as recursion.
- [ ] **REC2/SUM:** implement n>=2 equations with their actual order dependence, per-order repairs, compatible cutoff placement before curl, finite-order estimates and smooth summation/flat remainder.
- [ ] **WAVE/PHYS:** actual two pulse families, positive amplitudes and signed mean correction, finite-error averaged quadratic stress cancellation, corrected u/v/w and Cartesian NS residual, energy/tails and measured physical core contraction/aspect/swirl/material winding.

New gates: `current_original_conditional_cutoff_q_slow_jet_cover_executed`, `current_original_cutoff_local_supported_density_all17_declared_queries_executed`, and `current_original_cutoff_local_fixed_N_24_cell_C1_range_transport_executed`. The transport report separately records full24 range success. Actual controls, terminal closure, globally compatible N, coefficient recursion, oscillatory correction and full NS completion remain false. The long-term goal remains active.
