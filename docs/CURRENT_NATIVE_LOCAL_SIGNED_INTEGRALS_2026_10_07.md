# Actual native local signed Duhamel integrals

Checked implementation: [9863712f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/9863712f8c0df7f14eb0fc353dc036745b38a5fc). **Five actual nonzero signed local integral contributions are now executable over a whole active O2-slope source cell, including the whole axial interval Z in[.49,.51].** True spatial phase and original source functions feed the integrals. This advances beyond point density queries. It is one local C0 integral contribution, not the full cumulative history from r_minus to Rc.

## Exact domain, source and formula

Original O2_slope coordinate s in[.13369999,.13370001], exact log-radius width w=1/50000000=2e-8. The original chart has dy=ds, so no radius/Jacobian factor is missing. Two queries cover Z=.5 and the entire Z interval[.49,.51]. Candidate N=1024 is used in both geometry and density evaluation.

At each source query, the unchanged original radius offset gives phase frac(N*log(R/r_minus)). All returned phase cells are consumed. Native E,V,q, conditioned inverse and A/B enclose the source over the entire radial/Z cell, and the five signed density kernels enclose all those functions and phases. No sample or midpoint defines the integral.

For each original rate lambda_j in(1,3/2,3/2,1,0), the actual function is

```text
I_j(Z) = integral_left^right exp(-lambda_j*(right-s))*f_j(s,Z) ds.
```

The whole-cell signed range f_j is multiplied by the positive kernel mass w*integral_0^1 exp(-lambda_j*w*t)dt. Existing original exp_average supplies the directed positive mass enclosure without subtracting1-exp. Pressure has lambda_p=0, so its mass is w exactly. Candidate density kernels already contain their N/Pstar normalization; no extra1/N, radius, or unit multiplier is added.

When spatial phase wraps, all phase-cell kernels are hulled in the same original source/radius basis and arithmetic ledger. A shared upper log used to enclose a union is only an arithmetic coordinate; it does not select a field value. Different source-query bases must be explicitly rebound before future cumulative sums.

## Executed results and interface

Both actual radial/Z queries yield five nonzero contributions,10 integral enclosures in total. Signs on the whole tested sets: m,h,k,p are negative; e is positive. Thus the local source contribution is nonzero across Z in[.49,.51], not merely at Z=.5. Its magnitude remains in the existing logarithmic source factors. These local signs do not determine the sign/magnitude of a global oscillatory sum, where cancellation is possible.

```python
from lei_ren_part1_paper_compliant_current_native_local_signed_integrals import (
    NativeLocalSignedIntegrals)

# density_owner and binder share the already-live original seed.
integrator = NativeLocalSignedIntegrals(density_owner, binder)
result = integrator.contribution(
    Z=('.49', '.51'), left='.13369999', right='.13370001', N=1024)
I = result['contributions']  # five signed log-factored C0 integral enclosures
```

The same original common-unit P0 is recorded unchanged. The local integral has no assumed incoming defect history and does not reset any real incoming state. To obtain a history, transport the actual incoming value and add this local contribution. C1 Z derivatives and modified radial velocity/pressure/stress recovery are not supplied by this C0 API.

Evidence: fresh same-source two-query native replay with10 nonzero signed contributions;6 independent analytic positive-mass comparisons;12 independent signed affine-source integral comparisons; signed phase-union hull, microscopic positive/zero width, exact pressure mass and invalid input checks. 1027 working/index hashes pass. Producer 10.125s; focused checker 9.641s on the already-live native source. Read-only review: **GPT-5.6 Luna / max**; no material math/source-scope error found.

## Detailed next production tasks

- [x] **LEFT4c2-original-spatial-phase:** original17 affine radius offsets/Jacobians,16 radius/periodic-phase seam identities, exact inlet phase and candidate N*y binding on declared requests; see [spatial handoff](CURRENT_NATIVE_SPATIAL_PHASE_2026_10_07.md).
- [x] **LEFT4c3-local-active-C0:** actual nonzero signed local Duhamel contribution functions on the stated O2 radial/Z cell. Whole cumulative histories remain open.
- [ ] **LEFT4c2-q-slow-jets:** execute q_y,q_Z and required mixed/higher derivatives from the original correlated a,b,Delta roots. Use exact flat cutoff derivatives, positive active gamma and original cutoff derivatives. Retain tiny Delta=0 q; do not replace its function by0.
- [ ] **LEFT4c2-implicit-phase-jets:** use original Phi_psi positivity to recover inverse derivatives at fixed actual phi. Separate fast N*y derivatives from slow y/Z dependence. Recover A_y,A_Z,B_y,B_Z and needed higher orders in the conditioned coordinates with original hb/Pstar conversions exactly once.
- [ ] **LEFT4c3-density-C1:** differentiate every signed kernel and original expm1 factor in the same native basis. Supply whole-cell C1 Z covers; do not insert zero Z derivatives into GenericMomentRecovery.
- [ ] **LEFT4c3-common-cell-basis:** derive an explicit common factorization for adjacent cells so matching huge original amplitude/radius terms cancel symbolically before numeric rebasing. Sum transported cell contributions only after one shared basis/ledger is proved. Preserve signed intervals and tiny nonzero terms.
- [ ] **LEFT4c3-adaptive-active-integrals:** extend through a nontrivial O2 interval using true source/phase/conditioning subdivisions; transport each local contribution to its actual endpoint and enclose the total. Treat wraps as unions, unresolved q/u boxes as refinement requests and oscillatory cancellation with a controlled signed error. Record exact physical chart lengths and refinement/error evidence.
- [ ] **LEFT4c3-actual-incoming-history:** attach NativeBridgeSourcePackets.left_inlet() original five histories and unchanged P0 at the same r_minus, with the exact zero defect inlet. Quiet/flat collars transport these histories, including unchanged pressure memory. A zero local contribution must never reset incoming data.
- [ ] **LEFT4c2-whole-source-coverage-and-seams:** cover all17 original Z/radial charts for q, phase/A/B, velocities and kernels; bind same-function C1/higher traces at every source seam. Point queries and radius-only seam identities cannot close velocity/stress traces.
- [ ] **LEFT4c3-cumulative-Rc-functions:** compute all five actual signed histories and C1 Z jets from the true inlet through Rc. Preserve rates1,3/2,3/2,1,0 and absolute pressure datum. Commit executable functions and certified integration error, not only a bound table or phase fixture.
- [ ] **LEFT4d-Rc-targets:** actual five terminal target functions, A/A_Z normalization and divided(J-M)/mu in the reserved band. The division by tiny mu is part of the original recovery and cannot be skipped.
- [ ] **LEFT4d-unique-controls-and-field:** run the original reserved repair inverse, enclose unique control functions, prove the five terminal identities as Z functions and recover changed radial velocity/pressure/stress/continuation beyond2Rc.
- [ ] **HIGH / CONT / OUTER / ENERGY:** needed higher derivatives, modified same-function seams, heat-exterior pressure compatibility and finite energy across the changed family.
- [ ] **LEFT4e:** one genuine common finite N satisfying every source, derivative, integration and repair condition, then the global modified admissible stress cone. N=1024 here is not that admission.
- [ ] **REC:** actual n-dependent coefficient recovery with independent moment repairs and smooth sum. Coordinate scaling or a local integral is not coefficient recursion.
- [ ] **WAVE / PHYS:** mean/two-family oscillatory stress cancellation, flat forcing and corrected Cartesian uvw/p/f; independent residual and vortex/trajectory/scale diagnostics.

New scoped gate: `current_actual_native_C0_local_signed_Duhamel_contributions_executed`. Every global completion gate remains false. The important remaining blocker is now source-wide C1 cumulative histories and terminal repair, especially q/implicit-phase slow derivatives and stable integration across the original long charts. Continue bounded function production, mark only the achieved scope complete, and preserve unrelated edits.
