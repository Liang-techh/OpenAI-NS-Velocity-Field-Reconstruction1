# Original conditioned implicit phase and A/B first derivatives

> Successor: [CURRENT_NATIVE_DENSITY_C1_LOCAL_INTEGRALS_2026_10_07.md](CURRENT_NATIVE_DENSITY_C1_LOCAL_INTEGRALS_2026_10_07.md) ([d9ab0b82](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d9ab0b822a8766fdac95482003f626263fcc4770)) now executes same-source V_Z, all five signed density Z functions and actual local C1 integral contributions on the declared domains. Common cell basis/transfer, original incoming/P0_Z and full inlet-to-Rc histories remain open.

Checked implementation: [9be86964](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/9be869643f5379284db002c73a4330878f503bd8). **A_y, A_Z, B_y/Pstar, B_Z/Pstar and both fast phase derivatives now execute from the same original q/source jets. Actual spatial derivatives also execute at candidate N=1024.** There are14 original native query cells,112 primitive value/first-derivative enclosures and24 total spatial derivative enclosures. These include the whole local O2-slope radial/Z integral cell. Density C1 integrals and inlet-to-Rc cumulative histories remain open.

## Formula and numerical coordinates

The source comes directly from NativeQSlowJets, including the refined Delta=(a-2)+b^2/a and original eta. Values and first y/Z derivatives share one formal factor basis and directed arithmetic ledger. Original y denotes log radius; native width and Pstar conversions are not applied a second time.

For the conditioned loop, u=p2*q/dstar, h=sqrt(1+u^2), r=u/h, s=h^-2 and nu=1+t0^2+2q^2. Derivatives use (h^-1)_lambda=-u*h^-3*u_lambda and r_lambda=h^-3*u_lambda. The retained h^-1 and rho=1/[h(h+abs(u))] remain positive logarithmic factors. No rounded1-abs(r) defines a denominator.

The original direction uses the cosine Poisson profile, t=t0+2q*h^-1*(cos(psi)-r)/(1-2r*cos(psi)+r^2). The implicit phase derivative uses an analytic positive expression instead of the potentially cancelling dual derivative. In the E fraction chart,

Phi_E = [s+(t0*h^-1+2q*(cos(2pi*E)+r))^2]/[nu*(1+2r*cos(2pi*E)+r^2)] > 0.

The conservative positive lower bound s/(4nu) suffices. In the psi fraction chart Phi_psi=(1+t^2)/nu >=1/nu. The angle-map derivatives retain rho^2+4abs(r)*sin^2 or cos^2 denominators. At fixed fractional phase, x_lambda=-Phi_lambda/Phi_x; every primitive derivative is P_lambda|x+P_x*x_lambda and P_phi=P_x/Phi_x.

The small-r branch differentiates the original48-term W1/W2 Fourier sums. Explicit geometric bounds enclose the value and independent r,s,psi first partial tails, followed by the chain rule. It handles p2=0 across a changing source without replacing genuine derivatives by zero. A q-flat source is lazy/exact zero only with verified zero q slow rows. At phases0,1/2,1, primitive values and slow derivatives are exactly zero by symmetry, while fast phase derivatives are still calculated.

## Actual spatial chain and executed scope

The existing radius binder supplies phi=frac(N*log(R/r_minus)), with the unchanged original radius operator and original selected/analytic constants. Its phase is independent of Z. The candidate spatial chain is

```text
A_y_total = A_y|phi + N*A_phi
B_y_total/Pstar = B_y|phi/Pstar + N*B_phi/Pstar
A_Z_total = A_Z|phi
B_Z_total/Pstar = B_Z|phi/Pstar
```

Native queries:

- Five fixed-phase boxes at phi=.137 and Z=.5: inner_reference/O2_slope coordinate.1337, O2_buffer offset5.337, O3_slope_mu offset.537 and O3_power original offset.537 divided by the same Tw.
- The same five boxes with the actual derived spatial phase at N=1024.
- Three exact symmetry phases0,1/2,1 on O2_slope coordinate.1337,Z=.5.
- One whole actual O2-slope radial/Z cell: y=[.13369999,.13370001], Z=[.49,.51], using the actual spatial phase union on the full cell.

The two unresolved full-Z q boxes remain open. These successful queries do not establish whole-chart inverse/derivative coverage or a globally certified common N. Only first derivatives are installed; q's available y2Z1 rows do not automatically supply higher implicit phase jets.

```python
from lei_ren_part1_paper_compliant_current_native_phase_first_jets import NativePhaseFirstJets

# q_owner is the accepted same live NativeQSlowJets owner.
backend = NativePhaseFirstJets(q_owner)
result = backend.spatial_query('O2_slope', ('.49', '.51'),
    backend.ctx.mpf(('.13369999', '.13370001')), 1024)
for cell in result['cells']:
    A_Z = cell['values']['A_Z']
    B_Z_over_Pstar = cell['values']['B_Z_over_Pstar']
    A_y_total = cell['chain']['A_y_total']
```

Evidence: fresh same-source native replay of14 query cells/112 primitive rows/24 spatial chain rows;54 independent original scalar C0/first-derivative comparisons covering both signed r, small nonzero r, p2=0 crossing, transition, flat and symmetry fast derivatives. Scalar comparisons use100 dps original phase/primitive functions and finite differences with the disclosed numerical allowance; they are independent comparisons, not a theorem that the native global field is complete. Explicit first derivative tail bounds and positive phase formulas supply the enclosing arithmetic. 1024 working/index source hashes pass. Producer 30.781s; focused checker 33.500s on the already-live seed. Read-only review: **GPT-5.6 Luna / max**.

## Detailed next production tasks

- [x] **LEFT4c2-implicit-phase-first-jets on declared boxes:** same-source q rows; both angle coordinates; positive Phi_x; implicit y/Z derivatives; small-r derivative tails; genuine p2=0 crossing; lazy flat branch.
- [x] **LEFT4c2-A/B-first-jets on declared boxes:** original A/B value, slow y/Z and fast phi derivatives; exact symmetry values/slow derivatives with fast derivatives preserved.
- [x] **LEFT4c2-actual-fast-phase-first-chain on declared boxes:** actual radius phase and N=1024 chain, including the whole local integral radial/Z cell; no selector-width double conversion.
- [x] **LEFT4c3-original-normalized-V_Z on declared domains (see successor):** extract the ordinary axial Z derivative from the same packet/RadiusPolynomial and rebase it onto exactly the same basis/ledger as E,V,A,B. Preserve the original nonzero V; never infer V_Z from sampled values or reuse a different source query.
- [x] **LEFT4c3-exponential-increment-Z on declared domains (see successor):** deltaE=E*expm1(A/N); deltaE_Z=E_Z*expm1(A/N)+E*exp(A/N)*A_Z/N. Retain the existing factored expm1 identity for tiny nonzero A; bound exp(A/N) as a directed finite factor. deltaV=B/(N*Pstar), deltaV_Z=B_Z/(N*Pstar). Pstar is the same constant source normalization, with no Z derivative.
- [x] **LEFT4c3-five-signed-kernel-Z on declared domains (see successor):** differentiate m=deltaV; h=deltaE; k=V*deltaE+E*deltaV+deltaE*deltaV; e=2V*deltaV+deltaV^2-E*deltaE-deltaE^2/2; p=E*deltaE+deltaE^2/2. Keep every original cross term and sign. Return five C0 and five Z rows on a whole actual phase/radial/Z box, with explicit unresolved status if any source derivative is missing.
- [x] **LEFT4c3-local-C1-positive-weight-integrals on declared domains (see successor):** reuse the actual local cell and unchanged Duhamel rates1,3/2,3/2,1,0. Since original radius bounds/phase/positive weight are Z independent, integrate each signed density Z cover against the same positive mass. Store function-domain provenance and distinguish the local contribution from an entire incoming-to-terminal history.
- [ ] **LEFT4c3-common-adjacent-cell-basis:** define exact rebase identities for source factors and inherited Duhamel weights before adding neighboring cell contributions. Enclosures from different basis/ledger objects cannot simply be added. Do not replace a basis mismatch with a midpoint or unrelated radius cap.
- [ ] **LEFT4c3-original-incoming-histories/P0:** read the same five original incoming functions at r_minus and their Z jets. Bind the original absolute analytic pressure datum P0 once. Carry quiet-region memory unchanged, including pressure; no resets at chart seams.
- [ ] **LEFT4c3-full-chart-subdivision:** cover the current unresolved bridge_first/O2_axial full-Z cutoff boxes using retained Delta/eta and original coordinates; the known Z=[.49,.51] flat slices do not close full-Z coverage. Subdivide remaining active angle/phase boxes as needed, keeping smooth cutoff boundaries and true microscopic widths.
- [ ] **LEFT4c3-inlet-to-Rc-C1-functions:** integrate all five signed original kernels from actual inlet through original chart lengths to Rc; propagate incoming values and Z derivatives. Bound oscillation/phase/integration error without allocating astronomically many cycles. Report actual endpoint functions and both C0/C1 enclosures.
- [ ] **LEFT4d-actual-target/control-functions:** build actual Rc target A,A_Z and divided(J-M)/mu from these same cumulative functions. Apply the reserved independent repair inverse/control geometry. Check unique controls and all five terminal identities as Z functions, not selected samples.
- [ ] **LEFT4c2-higher-phase/A/B-orders:** extend ordinary implicit recurrences to the mixed orders required by C4 continuation, including cutoff and angle tail derivatives. Available q_y2Z1 alone does not establish these higher primitive jets.
- [ ] **HIGH/CONT/OUTER/ENERGY/LEFT4e:** establish higher source coverage, same-function seams, analytic pressure/heat/finite-energy compatibility, common finite N with derivative constants, and the global modified admissible cone.
- [ ] **REC/WAVE/PHYS:** actual n-dependent coefficient recovery with independent moment repairs and smooth sum; two-family oscillatory averaged quadratic stress cancellation; corrected Cartesian NS residual and multi-time scale/particle diagnostics.

For every task, mark completed only with the exact code/result/check commit, original domains and remaining limitations. Keep all global gates false until their own function-level requirements pass. Next production is signed-density C1 and local C1 integration; do not rerun unchanged q/C0 phase stages. Prior dependencies remain in the [q slow-jet handoff](CURRENT_NATIVE_Q_SLOW_JETS_2026_10_07.md) and [local integral handoff](CURRENT_NATIVE_LOCAL_SIGNED_INTEGRALS_2026_10_07.md).

Scoped gate: `current_original_native_conditioned_phase_and_A_B_first_slow_jets_executed`. Full coefficient recursion, global repair/cone and corrected Navier-Stokes field are not complete.
