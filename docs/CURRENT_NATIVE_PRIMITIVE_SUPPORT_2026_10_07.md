# Original primitive C0 support and improved fixed-N Rc targets

Successor: [CURRENT_NATIVE_COLLECTED_Q2_2026_10_07.md](CURRENT_NATIVE_COLLECTED_Q2_2026_10_07.md) adds direct original branch-local q squared mixed jets and a five-site original phase consumer, with complete24-cell continuous transport. Paired q/u derivatives, actual controls/compatible N/closure and recursion remain open.

Checked source implementation: [289be044](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/289be044808854a6a70042c284588d069f48a0db). This extends [CURRENT_NATIVE_CUTOFF_TRANSPORT_2026_10_07.md](CURRENT_NATIVE_CUTOFF_TRANSPORT_2026_10_07.md) on the same original24-cell route, full Z[-1,1] and candidate N2048.

Five original Rc target C0 absolute upper ranges now strictly improve; S_Z, Cp_Z and the divided joint target Z range also improve. The original M_Z and I_Z upper ranges remain unchanged. All24 continuous source cells and their C1 Duhamel transport still resolve. This is quantitative source-range improvement, not actual controls, terminal closure, a globally compatible N or coefficient recursion.

## Source theorem and implementation

The original modulation primitives satisfy

`A=a/2*(phi-psi/(2*pi))`,

`B/E=-(a*T1/(2*pi)+b*phi)/2`,

where T1 is the integral of the original direction t from0 to psi. For the original signed/small Poisson shape,

`mean(w)=0; mean(w²)=1/(2*(1-r²))=h²/2; mean(t²)=t0²+2q²=nu-1`.

The Fourier coefficients r^(n-1) give the same normalization for either sign of r. Cauchy-Schwarz implies |a*T1/(2*pi)|<=sqrt(a*(v-a)); b²<=a*(v-a), and sqrt(a*(v-a))<=v/2. On active support,

`v=kappa+sigma²*(2*eta-Delta)<=2+2*eta<=3`.

The flat branch has exact A=B=0. Therefore the original source has |A|<=5/4 and B=E*beta with |beta|<=3/2, uniformly over the original p2 and phase. These are source-function outer-range theorems, not values assigned to A, B or the velocity field.

New files use `experiments/root_st073/lei_ren_part1_paper_compliant_`:

- `current_native_Rc_range_loss.py` / `.json`: live original24-cell weighted target contribution inventory, source-to-common-basis diagnostics and original primitive maxima. The old target and primitive range bottleneck was active_first_bridge in every row.
- `current_native_bounded_primitive_transport.py` / `.json` / `_check.py` / `_check.json`: original C0 support ranges before nonlinear density, complete new continuous transport, comparisons and focused acceptance.

`signed_C0_support_range(original,bound)` keeps the original object when its absolute upper enclosure is tighter. Otherwise it uses the proved outer support range, retaining the original known sign and the factored B=E*beta source. A sign-definite range contradicting the theorem raises an error. Exact zero and genuinely tighter tiny formal factors are preserved. No bound endpoint or midpoint defines a source value.

`NativeBoundedPrimitiveDensity` retains the accepted cutoff/signed-u phase body. Before its nonlinear density calculation, it restricts only the original A and B C0 ranges. A_Z and B_Z are the same original objects. E/E_Z, V/V_Z, all cross terms, original phase/radius and q derivative rules stay intact. The report explicitly records when an A C0 arithmetic coordinate changes; it does not falsely claim its old formal coordinate remains untouched.

`NativeBoundedPrimitiveOracle` preserves original role/namespace/root/hash dispatch and separate P0/P0_Z. `NativeBoundedPrimitiveTransport` consumes the checked complete fixed-N baseline, uses unchanged exact geometry/true-width mass and carries every incoming correction through the same24 cells. The baseline source artifacts are unchanged.

## Measured range improvements and remaining limit

| Original Rc target | New C0 log absolute upper, approximately | Strict improvement |
|---|---:|---|
| M | 1.1769263e17 | yes |
| I | 1.1769263e17 | yes |
| S | 2.3538527e17 | yes |
| Cp | 1.1769263e18 | yes |
| D=(J-M)/mu | 1.1769263e18 | yes |

These are logarithms of conservative upper bounds, not actual target values. The previous corresponding upper logarithms themselves were of order10^408906090034569676. The improvement removes much of the avoidable C0 primitive-range explosion. The remaining bounds are still far too large to establish useful repair contraction at N2048.

The selected eta logarithm is approximately -3.490641234e408906090034569676, while the dstar logarithm is approximately -5.298317367. The remaining M_Z/I_Z bottleneck involves the extremely small original cutoff eta and separately bounded q/q_Z/u_Z factors. Some derivative sensitivity is genuine. Shared eta and q/u factors must also be collected before interval evaluation, so unrelated log-coordinate copies do not create artificial growth.

Focused acceptance:120 independent original A/B support comparisons, six independent signed/small Poisson mean-square comparisons and thirty independent density C0/Z formula comparisons. Synthetic fixtures are labelled separately from original field data. All240 live integral C0/Z ranges retain the original common context/ledger. Of138 branchwise C0 range decisions,64 use a tighter theorem range and74 retain the original formal range. Five C0 and three Z target upper ranges strictly improve; quiet pressure memory remains. Producer96.094s/checker1.875s. Working/index audit:1113 dependency hashes. Read-only reviewer model metadata: **GPT-5.6 Luna / max**, scoped PASS.

## Detailed next-agent work

- [x] **SOURCE3c1-baseline:** actual24-cell weighted target range-loss inventory and dominant original primitive/source cell identified. Native/common-basis ranges and downstream attenuation remain explicit.
- [x] **SOURCE3c1-C0-support:** original |A| and B/E support proofs applied to C0 ranges before nonlinear density. Original tiny factors and derivatives retained; actual five C0/three Z target upper ranges improve.
- [ ] **SOURCE3c1-C1-support:** add a separate branchwise ordinary Z support layer. Consume existing `serial.whole_period_C1`, `periodic_parameter_theorem` and `whole_cutoff_C1` prerequisites, with their actual original source/hash/family. Intersect only derivative outer ranges, retaining tighter original A_Z/B_Z objects. Do not infer derivative bounds by differentiating a C0 cap or zeroing a cutoff-crossing cell.
- [ ] **SOURCE3c1-q²-Z:** derive q²_Z directly before multiplying separate q and q_Z bounds. For the transition theta=Delta/eta and sigma=sigma(1-theta), use `q²=sigma²*eta*(2-theta)/(2a)` and `q²_Z=-(2*sigma*sigma'*(2-theta)+sigma²)*Delta_Z/(2a)-q²*a_Z/a`. Thus q*q_Z=q²_Z/2 and the defining eta/eta factor cancels exactly. Negative and flat branches need their own original formulas and endpoint consistency. Avoid interval cancellation of two independent copies of the same huge eta logarithm.
- [ ] **SOURCE3c1-paired-q/u:** before hulling, retain `(q*P)_Z=q_Z*(P+u*P_u)+q²*c_Z*P_u` and `(q²*H)_Z=q*q_Z*(2H+u*H_u)+q³*c_Z*H_u`, with u=c*q, c=p2/dstar, P=h^-1*W1 and H=h^-2*W2. Use original signed/small source identities, not caps as coefficient values. Preserve any genuinely linear q_Z sensitivity and report it separately from avoidable bound loss.
- [ ] **SOURCE3c2-joint-source:** build paired m/k target contributions before their independent per-key hulls. With F=E*A*g_N and B=B/Pstar, `f_k-Ac*f_m=[V*F+(E-Ac)*B]/N+F*B/N²`; its derivative includes V_Z,F_Z,E_Z-Ac_Z and B_Z. Keep the correct distinct k rate3/2 and m rate1. At source/integral level use `Dk*Mk*f_k-Ac*Dm*Mm*f_m`, and the corresponding Ac_Z term. One common mass for this difference is invalid. Preserve the final C_Z and division by mu*Ac².
- [ ] **SOURCE3c3-dominant-refinement:** refine active_first_bridge or a newly evidenced dominant cell only. Keep exact selected-sc/native endpoints, original global N*y phase, true microscopic widths and ordered incoming memory. Compare source/primitive and target bounds before/after; no point samples may define a continuous source or integral.
- [ ] **SOURCE3d-compatible-N:** retain separate original N^-1/N^-2 coefficient functions through the all-N route. Apply new support bounds with their same source dependence; fixed-N2048 ranges cannot be rescaled into coefficients. Report the original repair contraction inequalities and whether their required finite log N can be proved. A finite but huge upper range alone is insufficient.
- [ ] **CONTROL1a (parallel):** exact typed bump matrix, Gram/quadratic functions, inverse action and complete first-Z Picard map. Preserve `B_exact*h+N*r+Q_exact/N=0`; sidecar bounds remain separate from defining functions. Then solve actual five C1 controls and establish whole-Z terminal identities through the repaired field.
- [ ] **REC/HIGH/OUTER/WAVE/PHYS:** genuine n=1 and n>=2 recovery with independent repairs, higher jets/seams, remaining corrected cones/pressure/heat/energy, smooth summation/flat remainder, oscillatory stress cancellation and corrected Cartesian NS/physical winding diagnostics. Source-range improvement does not complete these tasks.

New gate: `current_original_supported_C0_primitives_before_fixed_N_C1_transport_executed`. All actual control/terminal closure/global compatible N/coefficient recursion/oscillatory/full corrected NS gates remain false. Keep the long-term goal active; avoid repeating the accepted24-cell coverage stage.
