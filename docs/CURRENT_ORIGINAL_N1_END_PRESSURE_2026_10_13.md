# Original end-pulse n1 pressure values — 2026-10-13

The actual end pulse, s=log(R/Rv) in [-4,0], now has source-bound m13/m15 value enclosures and 256 interval affine models covering every real Z in [-1,1]. The maximum absolute model errors in fixed positive source units are 2.4674448847e-5 for m13 and 2.3006530659e-5 for m15, both below .001, including through axial zeros. Seven independent queried positions are recorded: -1, -.639, 0, .211, .371, .639, 1.

These are errors of one radial segment in source-normalized pressure units. Complete physical pressure defects, n1 repair and full Navier–Stokes residual accuracy remain open. Full pressure derivative approximation accuracy is not installed.

## Original source and numerical meaning

The provider consumes the checked refined coefficient owner from [the selected coefficient checkpoint](CURRENT_ORIGINAL_RP_REFINED_COEFFICIENTS_2026_10_13.md). The same original incoming functions, complete future energy, beta normalization and nonzero mu/delta are used. Genuine C5 coefficient ranges on all four original axial cells supply values and derivative ranges at every model location; the model does not choose new profile coefficients or interpolate a replacement field.

Write Aphys=Pstar*U0*exp(logE-13/(2mu)-13), ubar=1/(1+Z^2), and Phi=10.99*sqrt(2pi)*(exp(5)-exp(3)). Reported m13 values use the positive unit Phi*Aphys*Rv; m15 values use Phi*Aphys*Rv^2. The original absolute logRc/logRp identities and original Rv=Rp*exp(13/mu) are bound to the accepted native source. Actual q=V/R is recovered from the physical U and own mean m=Mz/R; it includes the derivative of ubar and carries no extra sqrt(R/2).

Shared inverse-mu powers are combined symbolically before boxing. After expanding the same original logE, the linear source units have inverse-mu coefficients 11/2 and 37/2. The quadratic units have -4/2 and 22/2. In particular, before expanding logE, Rv^2*exp(-13/(2mu)-13) contributes 39/(2mu)-13. Physical Rv, N and the enormous end exponential are never materialized.

The source first row is exactly W_j(lambda), lambda=.5-mu. Its weighted coefficient sum is recovered from `-mu*(incoming1+ap*F1)`. Independent interval sums of C1 and C2 would destroy that small correlated result. The own-mean Fubini weight uses the finite-mu identity X_j=4A_j-D_j+rho_j; the positive nonzero rho_j integral and its proved bound are retained. W_j(lambda+1) uses the original beta J/J_r moments and a nonzero finite-mu Taylor remainder.

At both finite end points U is exactly zero. At the lower point q is nonzero and q_R=-q/R, so q+R*q_R cancels exactly in the unweighted viscous boundary. The weighted boundary remains in m15 as -qbase/Rv. The full quadratic reduced-S sources, including q^2 boundaries, U*q_Z, U^2 and U*U_Z integrals, are kept as nonzero directed error bounds based on genuine whole-band field ranges. The actual logA and logRv prove the finite positive error caps; neither contribution is replaced by zero.

Each affine approximation uses the exact midpoint of its intercept and slope endpoint intervals. Its error includes coefficient ambiguity, the genuine linear f''/2 range times half-width squared, and every nonzero full-source perturbation. The linear Taylor derivative ranges are auxiliary bounds; they do not establish uniform full-pressure derivative accuracy.

Reproduction: `python scripts/run_original_n1_end_pressure_source.py`. This loads checked ancestors without rerunning their math suites, computes this source quartet, runs its independent checker, and performs an accepted public query at global mp.dps=15. The provider is `CurrentOriginalN1EndPressureSource`; `evaluate(Z)` returns a live packet, `report(packet)` returns its source view, and `affine_models()` returns the continuous error models. Axial interval queries must fit inside one of the four original coefficient cells; the 256-model caller splits the whole domain accordingly.

The source quartet has stem `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_end_pressure_source`, with `.py`, `.json.gz`, `_check.py` and `_check.json`. Its receipt installs only `current_original_n1_end_pressure_source_installed`.

## Completed tasks

- [x] Bind the actual checked end coefficients, current angular inlet and original absolute radial scales.
- [x] Preserve exact first-row cancellation, finite-mu X remainder and both physical integration measures.
- [x] Recover the physical q, linear pressure values and finite lower viscous endpoint.
- [x] Retain genuine nonzero quadratic and weighted-viscous contributions as proved directed source bounds.
- [x] Build 256 adjacent absolute-error models over the entire real axial domain, including sign changes.
- [x] Run independent source/physical integration fixtures, model checks and a default-precision public query.

## Next tasks for repository agents

- [ ] Prefer a combined Rp-to-Rv source calculation where terminal own-mean identities avoid cancellation between separately boxed charts. With terminal M(Rv)=0, use `I0(m)=-integral Uz*log(R/Rp)dR` and `I1(m)=-Rp*M(Rp)-I1(Uz)`; independently verify the finite lower terms before accepting them.
- [ ] For the combined unweighted mean, evaluate the genuine saddle index derivative `F1'=integral (2+u)*original_kernel_1 du`. Its shared source-row formula is `I0(m)/(Aphys*Rv)=ubar*(ap*F1'+13*incoming1-G2-sum Cj*rho_j)`, with `G2=incoming2-incoming1+ap*(F2-F1)`. Prove the positive full kernel/tail enclosure; do not assume F1'=2F1 exactly.
- [ ] For the combined R-weighted main contribution, retain the original kernel with radial exponent `3/2-mu`, its saddle and exact common scale. Its contribution relative to the end unit is nonzero even when a source-log bound makes it tiny. Keep nonlinear main/gap terms in their actual correlated units instead of using a pointwise end-unit bound across the entire pulse.
- [ ] Evaluate the remaining entrance, main, exit, long gap and gap-end pulse contributions with the same positive selected branch, own mean and integer source-power ledger.
- [ ] Carry the original moment row cancellations across gap/end interfaces. Preserve linear, quadratic and viscous components in their own exact shared units; do not add independently boxed enormous logarithms or falsely cancel wide value boxes.
- [ ] Evaluate the remaining active pre-repair charts with source-normalized absolute models where signs change. Already accepted unchanged lower, O3, compact repair band, quiet and end pieces must each appear once.
- [ ] Restore complete signed m13/m15 with the accepted m11/m12/m14 functions. Distinguish a partial radial integral, its remainder bound and a closed global defect.
- [ ] Perform the original n1 two-plus-three repair with positive I2/Astar controls and the same normalized bump sources; recover its own M/V/P after the cutoffs and repairs.
- [ ] Establish functional terminal identities, axis regularity, interface smoothness and compatible analytic preheat pressure on the repaired family.
- [ ] Continue to genuine n2 and temporal recursion only after that n1 dependency. Admissible stress, oscillatory cancellation, flat remainder and independent corrected NS residual validation remain unfinished.
- [ ] Mark tasks done only with their actual evaluated source, independent receipt, dependency hashes and commit; update the five checkpoint/task entry files after the next physical layer advances.

Full original goal **ACTIVE / INCOMPLETE**. End-pulse pressure values are now available; complete n1 tails and repair are still required.
