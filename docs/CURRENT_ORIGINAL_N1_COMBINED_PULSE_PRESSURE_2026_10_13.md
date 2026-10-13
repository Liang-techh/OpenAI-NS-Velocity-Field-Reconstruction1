# Original combined Rp-to-Rv n1 pressure — 2026-10-13

The original main pulse, long gap and end pulse are now evaluated together over the full radial interval Rp<=R<=Rv. There are 256 adjacent interval affine models covering every real Z in [-1,1]. Maximum absolute errors in the fixed positive original source units are 1.7214685424e-5 for m13 and 2.3006530659e-5 for m15, both below .001. Seven point queries are recorded: -1, -.639, 0, .211, .371, .639, 1.

This advances the whole outer pulse, superseding the isolated end-segment pressure contribution when composing global defects. It does not close the remaining pre-Rp charts, complete n1 repair, pressure derivative accuracy, temporal recursion or the corrected Navier–Stokes residual. Full original goal **ACTIVE / INCOMPLETE**.

## Source and exact identities

The same checked C5 ap/C1/C2 functions, original incoming rows, finite nonzero mu/delta, original beta and main profile, energy and absolute radii are retained. With Aphys=Pstar*U0*exp(logE-13/(2mu)-13), ubar=1/(1+Z^2) and Phi=10.99*sqrt(2pi)*(exp(5)-exp(3)), m13 uses the unit Phi*Aphys*Rv and m15 uses Phi*Aphys*Rv^2. These are source-normalized pressure integral errors, not global physical NS residual errors.

Write I0(f)=integral_Rp^Rv f dR and I1(f)=integral_Rp^Rv R*f dR. From the actual terminal M(Rv)=0 and own mean m=M/R, finite-interval Fubini gives I0(m)=-integral Uz*log(R/Rp)dR and I1(m)=-Rp*M(Rp)-I1(Uz). The lower term is retained. Shared source identities give:

```text
G2 = rows[1]-rows[0]+ap*(F2-F1)
I0U/(Aphys*Rv) = -ubar*mu*rows[0]
I0m/(Aphys*Rv) = ubar*(ap*F1prime+13*rows[0]-G2-sum(Cj*rho_j))
I1U/(Aphys*Rv^2) = ubar*(sum(Cj*Wplus_j)+ap*Kmain)
I1m/(Aphys*Rv^2) = -I1U/(Aphys*Rv^2)-ubar*mu*rows[0]*(Rp/Rv)
```

The new checked saddle-index derivative evaluates F1prime=2*F1+integral(u*kernel_1) over the full original positive kernel support u in [0,11]. Both the positive finite-u core remainder and the accepted positive saddle tail remain. F1prime is not replaced by 2*F1.

The R-weighted main source uses its original exponent 3/2-mu. Its genuine support bound in the common end unit is below exp(-1000); the directed nonzero source enclosure is retained. The finite lower ratio Rp/Rv=exp(-13/mu) also remains nonzero. The unweighted viscous endpoint cancels exactly because q+R*q_R=0 at the flat pulse inlet. The m15 weighted endpoint remains -qinit/Rv.

Whole-pulse nonlinear errors use physical energy estimates over t=log(R/Rp) in [0,13/mu]. The actual initial mean and its first two Z derivatives are bound to the original incoming source. Unweighted Hardy and weighted finite-interval Cauchy estimates control the own mean, including its initial value. Bounds for q and q_Z give constants 288 and 7960; reduced pressure identities retain every quadratic volume and boundary source. All twelve shared-unit physical energy log bounds are below -1000. Pointwise bounds in the end unit are not extended across the main pulse.

Each affine error includes coefficient ambiguity, its genuine linear second-derivative Taylor remainder and the nonzero full quadratic/weighted-viscous source errors. Full-pressure derivative approximation accuracy remains false. Independent finite-interval moment/Fubini, main-kernel, Hardy/Cauchy and reduced-pressure fixtures, all models, saved/live source binding and a public query at global mp.dps=15 are checked. Previously accepted ancestor math suites are not rerun.

## Reproduction and completed tasks

Run `python scripts/run_original_n1_pulse_pressure_source.py` from the repository root. The checked index-derivative prerequisite can be reproduced with `python scripts/run_original_Rp_pulse_saddle_index_derivative.py`. The combined provider is `CurrentOriginalN1PulsePressureSource`; `evaluate(Z)`, `report(packet)` and `affine_models()` expose the source query and continuous models. Interval queries must stay within one of the four original coefficient cells.

The two source quartets have stems `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_index_derivative` and `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_source`, each with `.py`, `.json.gz`, `_check.py` and `_check.json`. Only the corresponding source-input and combined-pulse gates are installed.

- [x] Enclose the actual original saddle index derivative with its finite-u remainder and full positive tail.
- [x] Combine the original main/gap/end contributions using the exact terminal own-mean identities.
- [x] Keep the R-weighted main source, finite lower ratio, finite-mu end remainder and viscous endpoint.
- [x] Derive whole-pulse nonlinear energy bounds from the actual physical inlet and field derivatives.
- [x] Evaluate seven queries and 256 continuous axial models below .001 absolute error in fixed source units.
- [x] Bind independent checks and default-precision public queries to the exact source/raw dependency hashes.

## Next tasks for repository agents

Proceed in dependency order. Check the current source receipts and hashes once; reuse unchanged accepted parents. Do not replace a source function with a cap or its interval midpoint.

- [ ] Inventory the original pre-Rp radial domain and map every chart to one integration contribution. Identify exact joins and which existing lower/O3/repair-band/quiet pieces are already accepted.
- [ ] Bind each missing chart to the checked current physical U, own m and q, original finite N, original phase and native radial Jacobian. Preserve separate inner-mean and outer-pressure coordinates.
- [ ] Keep shared N and inverse-mu powers correlated before converting to log bounds. In U*q_Z and boundary terms, combine compensating powers symbolically rather than subtracting separately boxed enormous quantities.
- [ ] Determine whether a rigorous entire pre-Rp magnitude enclosure is below the error budget in the combined pulse's fixed positive units. If it is, retain the nonzero directed remainder; do not label it an evaluated signed center.
- [ ] Otherwise evaluate the missing pre-Rp pressure source integrals with their own means, original phase and both dR/R*dR measures. Produce absolute-error models across axial zeros.
- [ ] Prove finite q/q_R boundary terms at every retained endpoint. Do not discard a gap contribution solely because U is zero while its own mean remains nonzero.
- [ ] Add accepted pieces exactly once. The new combined Rp-to-Rv result already includes the end segment; do not add the older isolated end result again.
- [ ] Use exact meridional pressure support beyond Rv only for the quantities for which it was established. Preserve angular/pressure compatibility data independently.
- [ ] Compose complete signed m13/m15 with the accepted m11/m12/m14 functions in the original common family and physical units.
- [ ] Separate each evaluated integral, finite error enclosure and actual global defect. Keep the complete-tail gate false until full-domain composition is checked.
- [ ] Recover any axial derivative accuracy required by the original n1 repair using the genuine source jets and refined phase/own-mean bounds. Value accuracy alone is not derivative closure.
- [ ] Implement the original n1 two-plus-three repair with its positive I2/Astar controls and normalized source bumps, preserving all five terminal function identities.
- [ ] Recompute that repaired family's own M/V/P after cutoffs and repairs; do not inherit the unmodified pressure or means.
- [ ] Check axis regularity, all physical joins and high-order smoothness on the repaired family. Recover the compatible analytic preheat pressure datum and exact heat exterior.
- [ ] Install n1 repair acceptance only with functional terminal identities over the full axial domain and independent source-bound numerical evidence.
- [ ] Implement genuine n2 recovery equations and its own moment repair after n1. Retain n-dependent equations, coefficient definitions and the correct common inner domain.
- [ ] Extend coefficient and temporal recursion with finite-order remainder bounds; do not call coordinate rescaling alone a completed recursion.
- [ ] Establish admissible stress-cone margins and the separate flat remainder across the inner exit, matching/pulse/flatten/heat regions.
- [ ] Implement the two original oscillatory pulse families and their mean corrections; independently check averaged quadratic stress cancellation before claiming corrected NS accuracy.
- [ ] Measure the final corrected physical residual independently in Cartesian coordinates against the original force. Do not define force as the computed residual.
- [ ] When a task is completed, mark its item `[x]`, record the exact source quartet/receipt, error units, relevant commit and remaining open gates, and update all five checkpoint/task entry files.

For recovery after an interrupted agent turn, read this file and `CURRENT_CHECKPOINT.md`, inspect `git status`, and resume the first unchecked physical task. Keep existing computed raw artifacts and checked ancestors; do not restart completed source calculations merely because the conversation stream disconnected.
