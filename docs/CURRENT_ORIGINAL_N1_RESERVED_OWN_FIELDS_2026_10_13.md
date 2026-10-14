# Actual reserved n1 own fields and radial histories, 2026-10-13

Status: the source-bound construction installs C1 mean/M/K/P and five cumulative moment functions, plus genuine value functions for Q/V, on Rc <= R <= 2Rc. The full research goal remains **ACTIVE / INCOMPLETE**. Terminal zeros, higher derivative norms and global n1 closure are not certified by this milestone.

The actual selected alpha1/alpha2 and beta1/beta2/beta3 functions now enter their physical compact velocity patches and signed radial integrals. The previous coefficient numerical ranges remain enclosures of these functions; no interval center defines a field.

Let s=log(R/Rc), theta=Utheta1 and theta0=Astar*(2R)^(-1/2-mu). The leading axial velocity vanishes on each new patch support. The five patch densities per physical dR are:

```text
delta_m11 = Uz1
delta_m12 = 2R*F1 = sqrt(2R)*theta
delta_m13 = theta0*theta/R
delta_m14 = sqrt(2R)*theta0*Uz1
delta_m15 = -theta0*theta
F1 = theta/sqrt(2R)
dR = R ds
```

Every patch history is the actual signed integral from s=0 to the query s. All earlier bumps remain in the history after their local velocities vanish. The pure-power leading swirl is used only beneath the new patch supports, where the checked original old compact corrections and their derivatives vanish. It is not asserted to be the entire leading swirl on old supports.

The incoming velocity moment histories are the actual complete inner m11/m12/m14 functions. The incoming pressure histories are the true inner endpoint plus the lower-collar and active pre-Rc contributions, each once. In particular, the core pressure datum is preserved rather than reset or reintegrated.

The original fixed N1p forcing depends on the leading n=0 field. Its Rc-to-R contribution comes from the accepted compact pressure integrals with only their upper limit replaced by exp(s). The original physical Pstar node, dR/R*dR measures, nonlocal leading mean and finite viscous endpoints remain unchanged. At exp(s)=2, all twenty integral rows, sixteen endpoint rows and four resulting pressure rows reproduce the original accepted complete compact graph handles exactly.

The fixed source uses the original leading recovery weight 1-delta. The repaired n1 field has the distinct recovery:

```text
M1 = incoming_m11 + integral(Uz1 dR)
mean1 = M1/R
K1 = mean1-Uz1
mean1_y + mean1 = Uz1
M1_y = R*Uz1
Q1 = [2Z*Uz1-(1+delta)Z*mean1-(1-Z^2)*mean1_Z]/(1-delta*Z^2)
V1 = R*Q1
P1 = incoming_P1 + original_fixed_pressure_prefix + integral(theta0*theta/R dR)
```

The pressure patch changes the pressure source; it does not change the fixed N1p datum. `P1_y_patch` names only the patch contribution, never the total pressure radial derivative.

The current C1 coefficients provide C1 mean/M/K/P and the mean_Z needed for Q/V values. Q_Z/V_Z need genuine coefficient Z2 functions and are not fabricated. The unchanged checked cutoff continues to supply core/first fields; this milestone does not certify every global chart or join.

Public `at_band_fraction(q)` queries require an exact rational 0 <= q <= 1, with s=q*log(2). They return source-function handles in Z, not physical numerical values. Out-of-band queries fail. Direct symbolic handles are valid only on the documented chart.

Artifacts:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_reserved_own_fields.py`
- Same stem `.json.gz`, `_check.py`, `_check.json`.
- `scripts/run_original_n1_reserved_own_fields.py` checks the source-pinned archive without rebuilding ancestors.

Focused evidence includes the actual ten ordinary value/Z1 signed density identities, fourteen own recovery identities and the exact fixed-pressure prefix replay. An independent finite-mu fixture uses Rc=3.7, mu=.08, nonconstant coefficients/Astar, nonzero incoming mass/pressure and translated partial bump integrals at five radial points. Its 75 comparisons cover a support interior, inactive gap and full band endpoint. The fixture is mathematical validation, not original physical numerical accuracy.

Completed tasks:

- [x] Insert actual selected coefficients into the physical reserved velocity patches.
- [x] Construct ten signed value/Z1 history integral functions with the true dR measure.
- [x] Preserve all incoming mass and pressure histories without resetting between patches.
- [x] Construct the original partial fixed-pressure prefix and reproduce its full-band graph exactly.
- [x] Recover actual C1 mean/M/K/P and value-only Q/V with the n1 weight 1+delta.
- [x] Define five terminal integral functions using band endpoints plus the actual quiet/pulse pressure tails.
- [x] Enforce the exact public band domain and distinguish patch from total pressure derivative.
- [x] Check production identities and independent partial-history/FTC recovery cases.

Next tasks:

1. [ ] Reduce the actual terminal patch integrals by exact translated-bump changes of variables, then identify them with the accepted physical response matrix. Keep this as a proof about the integral fields, not merely the already-zero linear response expressions.
2. [ ] Prove all five transported terminal cumulative functions and their ordinary Z1 rows vanish. Derive the resulting mean/M/K/Q/V and pressure terminal statements only after those integral identities.
3. [ ] Extend genuine complete defect and selected coefficient Z2/Z3 functions, preserving existing C1 handles and all original source/forcing histories.
4. [ ] Upgrade own mean derivatives and recover Q/V Z1 from genuine mean_ZZ. Do not differentiate numerical interval endpoints.
5. [ ] Bind the complete pressure radial derivative to the same original N1p source plus patch contribution; retain the actual pressure gauge and finite endpoints.
6. [ ] Transfer required mixed/radial C3/C5 numerical bounds through the compact bumps and true radius factors.
7. [ ] Build and check global inner/reserved/quiet/pulse continuation and all required joins, then assess complete n1 closure.
8. [ ] Attach the compatible analytic preheat pressure datum, finite-energy tail and exact heat exterior.
9. [ ] Construct the admissible stress and flat remainder with separate physical norms and cone margins.
10. [ ] Implement genuine n>=2 coefficient equations, n-dependent recovery and independent moment repair on the common core interval.
11. [ ] Bound finite-order remainders and smooth summation; implement oscillatory/mean corrections and averaged stress cancellation.
12. [ ] Independently validate corrected Cartesian NS residuals in physical Linf and volume-L2 against 1e-3 with restricted forcing.

The milestone provides actual constructed field functions. It does not claim original numerical own-field accuracy, terminal zeros, complete n1 closure, temporal recursion or corrected NS validation.
