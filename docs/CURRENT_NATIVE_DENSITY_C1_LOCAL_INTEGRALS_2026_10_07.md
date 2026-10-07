# Original signed density Z functions and actual local C1 integrals

Checked implementation: [d9ab0b82](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d9ab0b822a8766fdac95482003f626263fcc4770). **All five original signed density kernels now have genuine Z derivative functions, and the actual local signed Duhamel contributions now have C1 Z enclosures.** Seven native source queries give35 signed-density values,35 signed-density Z rows,10 local integral values and10 local integral Z rows. The full inlet-to-Rc histories and terminal repair remain open.

## Same-source velocity and derivative construction

The backend consumes the accepted NativePhaseFirstJets owner. Its spatial query returns one original q/source packet shared by every derived phase cell. E/E_Z come from that same packet's original roots. V/V_Z are extracted from packet.velocity['axial'][0] with ordinary_axial_coefficient(k=0,1), wrapped in RadiusPolynomial and evaluated on exactly the same formal source basis/ledger as E,A,B.

The packet's axial velocity already has the original common Pstar normalization. Taylor coefficient_k is multiplied by k! to obtain the ordinary Z derivative. No native selector derivative, radial half shift or second Pstar normalization is introduced. Source factor and radius covers enclose analytic functions; no interval endpoint/midpoint becomes a field value.

For candidate N, let dE=E*expm1(A/N) and dV=(B/Pstar)/N. The new derivative functions are

```text
dE_Z = E_Z*expm1(A/N) + E*exp(A/N)*A_Z/N
dV_Z = (B_Z/Pstar)/N
m_Z = dV_Z
h_Z = dE_Z
k_Z = V_Z*dE + V*dE_Z + E_Z*dV + E*dV_Z + dE_Z*dV + dE*dV_Z
e_Z = 2*V_Z*dV + 2*V*dV_Z + 2*dV*dV_Z
      - E_Z*dE - E*dE_Z - dE*dE_Z
p_Z = E_Z*dE + E*dE_Z + dE*dE_Z
```

Every original nonzero-V and quadratic cross term is retained. The factored expm1 identity preserves unmaterializable tiny nonzero A. The directed exp(A/N) coefficient is evaluated only after the existing abs(A/N)<=1 finite-source guard succeeds. This does not install a general huge-positive-log exponential backend. Exactly flat primitive/slow jets yield exactly zero density changes and Z rows.

The actual phase is frac(N*log(R/r_minus)). The original radius phase is independent of Z, so the phase-held A_Z/B_Z rows are already the total spatial Z derivatives. No extra N*phi_Z term is added.

## Actual local C1 integral functions

The local region remains the original O2-slope coordinate[y_left,y_right]=[.13369999,.13370001], with exact width1/50000000, dy=dcoordinate, and Z=.5 or the whole Z=[.49,.51]. Candidate N is1024. The actual phase union and density/source bounds cover the full radial/Z cell.

```text
I_j(Z)   = integral_left^right exp(-lambda_j*(right-s))*f_j(s,Z) ds
I_j_Z(Z) = integral_left^right exp(-lambda_j*(right-s))*f_j_Z(s,Z) ds
```

The source and cutoff are smooth on the compact covered cell. The original integration endpoints, radius phase and positive weights are Z independent. Thus differentiation passes under the finite integral with no boundary terms. Rates remain1,3/2,3/2,1,0 for m,h,k,e,p.

A separate derivative-kernel hull is constructed from the same phase cells and basis, then multiplied by the same positive kernel mass. The C0 hull itself is never differentiated. Mass times a signed source hull is an outer integral enclosure; it is not an exact constant-density integral identity. C0 and Z covers are separate and do not provide a joint sign/correlation theorem.

All ten original C0 local integrals remain nonzero, with m/h/k/p negative and e positive. At Z=.5 the derivative enclosures have m/h/p positive and k/e negative. On the whole Z interval, h/p remain positive and k/e negative; the m derivative cover spans zero and is recorded as sign-undetermined. That is an enclosure limitation, not an invented zero derivative.

The p row covers only the pressure-history increment. The original absolute P0 is stored unchanged and not reset. Recovering a total absolute-pressure Z derivative later must add the independent original P0_Z path.

## Executed scope and API

- Five actual spatial boxes at Z=.5: inner_reference/O2_slope coordinate.1337, O2_buffer offset5.337, O3_slope_mu offset.537 and O3_power original offset.537 divided by the same Tw.
- Two whole local radial cells at Z=.5 and Z=[.49,.51].
-35 signed C0 plus35 signed Z density enclosures;10 signed C0 plus10 signed Z integral enclosures.
- Exact lazy zero density changes in the two declared O3 flat queries; tiny nonzero O2 buffer changes retained.

```python
from lei_ren_part1_paper_compliant_current_native_density_C1_local_integrals import NativeDensityC1LocalIntegrals

# first_owner is the accepted same live NativePhaseFirstJets.
backend = NativeDensityC1LocalIntegrals(first_owner)
result = backend.contribution(Z=('.49', '.51'),
    left='.13369999', right='.13370001', N=1024)
I_k = result['contributions']['k']
I_k_Z = result['Z_derivatives']['k']
```

Evidence: fresh same-source replay of seven queries;35 native density Z rows and10 actual integral Z rows;100 independent scalar C0/Z comparisons formed from full modified velocity products/squares instead of the expanded density formulas;36 independent analytic Z-dependent integral C0/Z comparisons; tiny factored increments/derivatives and flat zeros; exact original width/positive mass/source-basis provenance. 1030 working/index hashes pass. Producer 18.875s; focused checker 19.578s on the already-live source. Read-only source/math review: **GPT-5.6 Luna / max**.

## Detailed next production tasks

- [x] **LEFT4c3-original-normalized-V_Z on declared domains:** same packet and ordinary Taylor conversion; original nonzero V retained; no double unit/width conversion.
- [x] **LEFT4c3-exponential-increment-Z on declared domains:** factored tiny expm1, bounded exponential chain and original normalized B_Z/N.
- [x] **LEFT4c3-five-signed-kernel-Z on declared domains:** all five original C0/Z kernels with every signed/cross term.
- [x] **LEFT4c3-local-C1-positive-weight-integrals on the actual local cell:** C0 and separate genuine Z covers, actual phase union, exact width and original positive weights.
- [ ] **LEFT4c3-common-adjacent-cell-basis:** build a checked common source-family coordinate basis for neighboring cells. Bind identities for the four fixed source logs and original varying radius factors; preserve directed logarithmic covers. Rebase whole C0/Z enclosures before addition, without selecting midpoint fields or treating a radius cap as the true radius.
- [ ] **LEFT4c3-local-serial-C1-Duhamel-transfer:** on a bounded adjacent-cell partition, accumulate I_next=exp(-lambda*width)*I_incoming+I_cell for each original rate and the same transfer for Z rows. Test both signed contributions and quiet cells, keep rate0 pressure memory. Record exact endpoint/partition/domain provenance and explicit remaining integration error.
- [ ] **LEFT4c3-original-incoming-histories:** at the actual left inlet r_minus=Ra*exp(hb*s_c/2), extract the original five incoming functions and ordinary Z rows directly from the unchanged native packet. Do not initialize all histories to zero unless the exact original function theorem gives zero. Keep incoming correction support separate from original background memory.
- [ ] **LEFT4c3-P0-and-P0_Z:** bind original absolute analytic preheat datum and its ordinary Z derivative from the same packet. Add it only in pressure recovery, with correct common normalization; local p increments alone are not absolute pressure or its derivative.
- [ ] **LEFT4c3-original-chart-length-and-transfer:** install actual positive log-radius lengths/Jacobians for bridge/switch microscopic charts, macro transition, reshape, restoration/patch and O2/O3 pieces. Use existing original radius-offset identities. Quiet segments still transport nonzero inherited histories.
- [ ] **LEFT4c3-whole-source-cutoff-subdivision:** resolve remaining bridge_first/O2_axial full-Z active/flat cutoff coverage with retained Delta/eta and original coordinates. Existing Z=[.49,.51] flat slices do not close full-Z coverage. Return explicit unresolved function domains without fabricated derivatives.
- [ ] **LEFT4c3-broad-integral-oscillation-control:** replace microscopic-only cell coverage with efficient whole-chart original oscillatory integration or paper-consistent cancellation estimates. Preserve N*y correlations and true lengths; do not allocate astronomical cycles or average away signed terms without a proof.
- [ ] **LEFT4c3-actual-inlet-to-Rc-C1-histories:** propagate all five original incoming functions plus the modified signed contributions to Rc with both C0/Z enclosures. Account for source/phase/integration error and chart seams; report function-domain results, not sample-only terminal values.
- [ ] **LEFT4d-actual-target-functions:** construct Rc target A,A_Z and divided(J-M)/mu from these same histories and original axis datum. Keep inverse-mu scaling factored and separate reference/target roles.
- [ ] **LEFT4d-reserved-independent-repair-controls:** use the accepted reserved annulus geometry and repair inverse to produce unique controls as functions of Z. Verify all five actual terminal identities and propagate changed radial velocity/pressure/stress downstream.
- [ ] **LEFT4c2-higher-phase/A/B-and-C4-source:** extend ordinary mixed phase/primitive recurrence and series tails to orders required for continuation. First Z density functions do not establish C4 regularity.
- [ ] **HIGH/CONT/OUTER/ENERGY/LEFT4e:** same-function higher seams, analytic pressure and exact heat tail, finite energy, derivative constants/common finite N and global modified admissible cone.
- [ ] **REC/WAVE/PHYS:** n-dependent coefficient equations with independent repairs/smooth sum; two-family oscillatory averaged quadratic stress cancellation; corrected Cartesian NS/scale/particle diagnostics.

Mark a task complete only with its exact code/result/check commit and declared function domain. Next production is common adjacent-cell C1 transfer and original incoming histories, then full inlet-to-Rc accumulation. Do not rerun unchanged q/phase/local C0 stages. Global completion gates remain false. Dependencies and higher derivative requirements remain in the [first-jet handoff](CURRENT_NATIVE_PHASE_FIRST_JETS_2026_10_07.md).

Scoped gate: `current_original_native_signed_density_Z_and_local_C1_Duhamel_functions_executed`. Local C1 contributions do not close five terminal moment identities or implement actual coefficient recursion.
