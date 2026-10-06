# Current modified O2/O3 signed paper stress and source remainder

Implementation and actual source receipts: commit [db84cfa3](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/db84cfa3d1c17e5667bf403fac0307a3c691cac8).

## Constructed result

The actual repaired O2 buffer/O3 variable-slope/quiet-power field now has **callable full signed original-paper theta and axial stress columns through ordinary logR order3, and all three source-remainder components through ordinary logR order2**. The construction uses its own five histories, retained incoming data and same-axis absolute pressure. It retains15 stress sectors (theta7/axial8) and11 remainder sectors (radial9/theta1/axial1).

The source is the preceding [own repaired-history milestone](CURRENT_O3_REPAIRED_HISTORIES_2026_10_06.md), including the same exact implicit repair controls. Physical Uz=Pstar*Vhat is retained. No Pstar or radius cap becomes a field value. All original variable amplitude derivatives, linear moment terms, nonlinear products and radial/axial viscosities remain signed.

At q=2, the new signed stress/remainder source equals the original source for all Z[-1,1] through all common retained axial rows and advertised ordinary/mixed derivative orders. Higher available axial rows of the split radial Pstar^0 product are retained; a zero term in the old combined expression lowered its Taylor projection order. The checker compares the common projection without modifying stored source rows.

This constructs the local paper columns/remainder. It does **not** yet install the completed diagonal/Cartesian tensor, physical source packaging/divergence, open q>2 continuation, modified source/physical/heat joins, independent O2 taper or modified closed O3 cone/common N. It also does not prove remainder temporal flatness or the full corrected NS equation. Original counts remain15 strict nonzero whole regions plus separate exact zero exterior, with17 original regions open. Inventory33/32/14 and primitive atlas14/8 are unchanged and separate from modified-source admission.

## API, files and source lineage

`CurrentModifiedPreStress(source=checked_repaired_histories).stress(region,Z,coordinate)` supports the parent's exact local domains: O2_buffer offset[0,11], O3_slope_mu offset[0,1], quiet_O3_power q[0,2]. The focused controller stage is `currentmodifiedprestress`.

Files use the prefix `experiments/root_st073/lei_ren_part1_paper_compliant_current_modified_pre_stress`: `_operator.py`, `.py`, `_check.py`, `.json`, `_check.json`, `_views.json.gz`. The complete compressed views retain the actual source, signed stress/remainder, whole O2/O3/repair/full-Z boxes, three fresh bump interiors and exact exit. A further fresh checked query is evaluated directly.

The checked parent chain remains: original current graph -> finite-N local modulation -> own five modulated histories -> independent source-correlated implicit repair -> own repaired histories -> signed modified stress/remainder. Reuse live session64009 on this host when available; session IDs are host-local. Avoid cold all-stage rebuilding. On another host, reconstruct only the required checked dependency chain once.

## Exact amplitude sectors and original paper columns

Let u=Utheta/Pstar, Vhat=Uz/Pstar, and use the actual ordinary source rows. The modified normalized moments are `m=m0+Pstar*m1`, `k=k0+Pstar*k1`, while h,e,p retain normalized Pstar^0 rows. They represent physical moments:

`Mz=R*m`, `Mtheta=sqrt(2)*R^(3/2)*Pstar*h`,

`Mtheta_z=sqrt(2)*R^(3/2)*Pstar*k`, `Mztheta=R*Pstar^2*e`, `pressure=Pstar^2*p`.

For d=1-Z^2, L=1-delta*Z^2 and transport(m)=(1-delta)*Z*m+d*m_Z, the ordinary stress source has the same full original formulas as `current_pre_pulse_stress_operator.py:raw_pre_stress_rows()`. Only exact radius/constant-Pstar factors are split. Each part's physical source factor is `R^mode[0]*Pstar^mode[1]/sqrt(2)`. Shifted logR rows include the exact radius exponent; they do not replace local source derivatives by a fixed slope.

| Column and term | Radius power | Pstar powers | Count |
| --- | --- | --- | --- |
| theta local transport | 1/2 | 1 | 1 |
| theta retained angular moment | 1/2 | 1 | 1 |
| theta retained mixed k moment | 1/2 | 1,2 | 2 |
| theta u*transport(m) | 1/2 | 1,2 | 2 |
| theta actual variable radial shear | -1/2 | 1 | 1 |
| axial local Uz | 1/2 | 1 | 1 |
| axial Uz*transport(m) | 1/2 | 1,2 | 2 |
| axial retained linear m | 1/2 | 0,1 | 2 |
| axial retained own Mztheta | 1/2 | 2 | 1 |
| axial actual absolute pressure | 1/2 | 2 | 1 |
| axial radial shear | -1/2 | 1 | 1 |

The source m/k Pstar^1 sector is essential: dropping it loses the extra Pstar^2 theta/mixed transport and the axial Pstar^1/Pstar^2 contributions. Feeding Vhat as raw Uz loses a Pstar factor. Applying another Pstar to an already factored mode double counts it.

## Full signed remainder sectors

The actual radial rows are Ur0+Pstar*Ur1; theta and axial rows have one Pstar factor. The new operator applies the original full remainder coefficient programs to those rows before collecting equal powers. Linear radial terms retain powers0,1. Both ordered radial cross terms are retained in the power1 nonlinear coefficient. Local axial transport `Uz*d_z(Ur)` contributes powers1,2.

| Remainder term | Pstar powers | Original physical beta | Half normalization |
| --- | --- | --- | --- |
| radial time | 0,1 | -3 | yes |
| radial radial viscosity | 0,1 | -3 | yes |
| radial quadratic/axial transport | 0,1,2 | -3 | yes |
| radial axial viscosity | 0,1 | -3+2delta | yes |
| theta axial viscosity | 1 | -3+delta | no |
| axial axial viscosity | 1 | -3+delta | no |

Exact source radius powers and beta/half modes match the original remainder. Theta/axial stress/remainder vectors remain signed throughout; absolute-value bounds are not new field values. Actual frequency-dependent local rows remain, including their N derivative growth. There is no cone or flatness claim from small profile amplitude alone.

## Evidence and scoped acceptance

The symbolic proof replays the actual source programs, not handwritten replacement formulas. It binds the checked parent `.history()` to both new operators, consumes the original arbitrary-variable full-paper stress AST theorem (20 identities), proves 81 new signed stress/remainder source identities, and compares the unchanged original full-remainder coefficients. It retains the correct physical Uz, all radial powers and products, and every original radius/beta/half mode.

Evidence also includes complete whole/full-Z producer views, q=2 common source projections, retained extra axial rows, fresh API query, checked constructor/receipt, focused controller and compilation. Working/index hashes matched 846 dependency files. Read-only reviewer: GPT-5.6 Luna / max; no missing term or incorrect Pstar mode found. The initial representation-length mismatch was fixed by comparing all common retained coefficients while preserving the longer source rows.

Only three new gates are admitted: full signed theta/axial paper stress3 source, full signed three-component remainder2 source, and actual-source Pstar modes. Own physical energy, completed tensor/divergence, physical/heat joins, modified cones/common N/global corrected NS remain false. N=10^12 still meets only the repair threshold27,303,666; no common cone-frequency admission is made.

## Correct next construction dependency

The completed diagonal tensor is determined by the physical axial stress: `Ttheta_theta = r*d_z(Trz)`. In `pulse_end_physical_C2.py:lift_physical_packet()`, `diagparts['completed_radius']=logR/2+log(2)/2` is retained exactly. It is not reconstructed from physical kinetic energy, and it is not a radius cap. The diagonal cancels radial tensor divergence by the original source identity.

Physical total kinetic energy is a separate required computation. The profile moment Mztheta is already in the own repaired five histories and is not the total kinetic energy. Postpulse complete-future swirl energy uses the published future/2 normalization in `current_postpulse_energy_history.py`; it does not account for the modified O3 radial/axial kinetic terms. Therefore TENSOR4 physical tensor completion need not wait for TENSOR3 global kinetic-energy integration.

## Detailed next tasks

Mark a task complete only after connecting the source operator, artifact/receipt and callable path. Record the commit and remaining scope. Preserve unrelated files. The full NS/recursion/wave/dynamics objective remains active.

- [x] **REPAIR1–4:** fixed independent compact supports, correlated new inverse, complete five-moment nonlinear map and unique exact implicit constant controls.
- [x] **REPAIR5–8:** actual bump ordinary logR4/axial5 rows, continuous partial primitives, own repaired five histories, same-axis pressure and own-moment radial source.
- [x] **REPAIR9a:** actual-source q=2 functional five-moment/velocity/pressure/radial exit through retained orders.
- [x] **TENSOR1:** full signed original theta/axial paper columns with every actual Pstar sector and variable source row.
- [x] **TENSOR2:** full signed source remainder2, all radial/axial quadratic products, actual beta/radius/half modes and common functional q=2 trace.
- [ ] **TENSOR4a — Source-bound physical packet.** Package these15 column sectors as `full_meridional_stress_log_sectors`, with ordinary source mixed3 rows and exact radius/Pstar logarithms. Use the original geometry owner from the same checked registry for O2/O3; quiet q is log-radius offset, not normalized O3 phase. Convert quiet q to phase=q/Tw only when invoking the original radius dispatcher. Keep the original actual pressure and same implicit repaired-source hash. No radius/amplitude cap substitution.
- [ ] **TENSOR4b — Physical remainder adapter.** Adapt `lift_physical_packet()` only at its `source_errors` input to consume these11 factored source sectors; leave physical operators/factor logs/cylindrical-Cartesian maps unchanged. Each mode[1] multiplies the actual logPstar. Do not collapse Pstar sectors into interval values or feed an unmodified velocity into the remainder program. Prove the adaptation structurally and retain actual beta/half exponents.
- [ ] **TENSOR4c — Completed diagonal/divergence.** Consume the unchanged arbitrary-source full physical operator proof and build `Ttheta_theta=r*d_z(Trz)`, with exact completed-radius factor and all axial sectors. Export physical cylindrical stress3, divergence2, diagonal2, Cartesian completed tensor and all three remainders. Prove radial divergence cancellation on the modified source. Do not claim the tensor already admissible merely because completion is installed.
- [ ] **TENSOR3a — Actual kinetic density.** Compute radial, theta and axial squared velocity from the real sectors. Include both Ur0*Ur1 cross terms and Pstar^2 radial/axial terms, the original radial sqrt(R/2) factor, exact physical map and volume Jacobian. Distinguish this nonnegative total density from signed Mztheta.
- [ ] **TENSOR3b — Continuous integrals and finite energy.** Implement actual partial/future kinetic integrals across modulation/repair, retain all radial/axial/cross terms, and combine them with unchanged downstream swirl energy only through the proven exact source exit. Preserve future/2 normalization when that published profile is consumed. Prove forward/backward definitions are one function and report physical radial-tail/global-time bounds; original interior kinetic energy cannot be copied after velocities change.
- [ ] **TENSOR5a/REPAIR9b — Open q>2 continuation.** Define exact original power continuation to actual Tw, justified by common q=2 source traces and flat differences. Parent quiet API presently ends at2. Complete the modified source/physical dispatch on all affected charts; keep original inventory counts distinct from modified chart admission.
- [ ] **TENSOR5b — Completed exit and downstream heat joins.** Prove completed tensor/remainder/physical spatial4/time1 exit equality, then inherit original downstream analytic pressure/heat source through exact same-function continuation. Do not infer whole heat compatibility from a local numeric endpoint comparison.
- [ ] **TENSOR6 — Affected O2/O3 interfaces.** Prove buffer entrance/taper boundaries, O2-to-O3 transition, slope-to-quiet power, all support edges and repair exit, using one actual modified source and actual lambda/time maps. Split only where charts/flat source definitions require it; sampled matches are insufficient.
- [ ] **TENSOR7 — Phase-aware derivative estimates.** Bound true NlogR modulation derivatives and independent repair jets, pressure/history terms and signed shear contributions. Retain their N growth and small positive mu. Identify exact correlations/cancellations before enclosing; pointwise O(1/N) profile smallness does not control derivatives.
- [ ] **REPAIR10/MOD7 — One finite common N.** Combine repair, derivative/shear, pressure, completion, O2 taper, interface and full two-vector cone inequalities on fixed supports. Publish one sufficient finite integer and its actual-source proof. Repair-only threshold is not enough.
- [ ] **REPAIR11a — Whole O2 taper cone.** Use complete actual signed stress/pressure/moments/remainder and all Z including flat edges; prove independently on the modified field.
- [ ] **REPAIR11b — Closed modified O3 cone.** Address the original zero-shear endpoint obstruction with the new source. Prove both strict original signed two-vector inequalities throughout transition/repair/continuation, with full components and common N. Do not inherit original-source cone flags.
- [ ] **GLOBAL1:** remaining O2/Rh, patch/restore/reshape/switch/bridge/core cones; one smooth global admissible completed tensor/physical lift and compatible heat exterior.
- [ ] **REC1:** actual n=1 coefficient recovery, correct source equation/common core domain, independent per-order five-moment repair and nonflat leading-origin cancellation.
- [ ] **REC2:** actual n>=2 equations and independent repairs; cut off streamfunction/vector potential before curl to retain exact divergence-free fields.
- [ ] **REC3:** finite-order error bounds, temporal flatness and valid smooth coefficient sum; coordinate scaling of the leading field is not coefficient recursion.
- [ ] **WAVE1/WAVE2:** two real oscillatory families with positive amplitudes/flat weights, covariance and finite-frequency bounds, divergence-preserving lifts, mean corrections and averaged quadratic cancellation of the admissible background stress.
- [ ] **PHYS1/PHYS2:** resolved common u/v/w/p evaluators in similarity/physical coordinates; independent corrected Cartesian forced NS Linfinity/L2<1e-3, finite kinetic energy/tail, axis regularity and smooth forcing.
- [ ] **DYNAMICS:** measure core radial contraction, axial/radial aspect ratio, swirl/axial/vorticity scaling, true coefficient-level recursion and material-particle winding; distinguish streamlines from accumulated winding. Animation remains secondary.

Immediate next owner: TENSOR4a–c actual-source physical/completed tensor packaging. TENSOR3 physical energy is a separate parallel branch. The full long-term objective is not complete.
