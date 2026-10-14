# Actual analytic core C5 and lower pressure C3, 2026-10-13

Status: required derivatives of the same admitted analytic nonlinear core now feed original lower-collar U/mean C5, Q/Q_y/Q_yy C4, and signed m13/m15 C3 integral functions. All earlier C3/C2/C1 source handles and original histories remain. Together with the [pulse pressure C3 milestone](CURRENT_ORIGINAL_N1_PULSE_PRESSURE_C3_2026_10_13.md), this removes another real dependency for complete five-defect/coefficient derivatives. Full goal remains **ACTIVE / INCOMPLETE**.

The higher core components use the already admitted unique analytic solution of the actual ten angular and ten axial nonlinear terms with the same Green inverse. Only ordinary slow-Z derivative order changes; quantity, rho, radial order, family, source providers, axis conditions and nonlinear definition are guarded. Phi/Psi/mean_Psi radial orders 0, 1 and 2 expose Z0..Z5, preserving all existing C3 handles. An altered component definition is rejected. No second core solution or numerical enclosing center is selected.

The anchored gauge retains its exact H-root lower endpoint and moving upper endpoint Z. Its higher derivatives follow the true FTC:

```text
G_Z = gradient
G_Z2 = gradient_Z
...
G_Z5 = gradient_Z4
F0 = exp(-logC-Lambda*G)
```

The higher F0 rows are checked against independent exponential differentiation. The exact nonzero Lambda and source root are retained. Fixed-endpoint differentiation is not incorrectly applied to G.

The original collar continues its ratio ODE, Psi transport and own mean from the true rho=4 inlet. The nonzero incoming mean and nested integrals remain; none is reset at the collar. Exact binary sc and original nonzero hb determine the original sc/4..sc/2 domain. Physical radius and measures are:

```text
rho=4*exp(hb*s), R=epsilon*rho
dR=epsilon*hb*rho*ds
R*dR=epsilon^2*hb*rho^2*ds
```

U/U_y/m/m_y/Utheta and radial source rows through rho order 2 gain C5. Q/Q_y/Q_yy and Q_rho rows gain C4 from the full mean_Z5. The leading recovery uses 1-delta. It does not use the repaired n1 1+delta weight.

The direct pressure forcing preserves both viscous radial derivatives:

```text
Q_y=rho*Q_rho
Q_yy=Q_y+rho^2*Q_rhorho
N1p=-S/2+(Q_yy+Q_y)/R=-S/2+2*Q_R+R*Q_RR
m13=integral N1p dR
m15=-integral R*N1p dR
```

The new pressure Z2/Z3 functions differentiate the actual signed source bodies at the same finite, Z-independent endpoints. Physical units are applied once. The original minus in m15 is applied once. The original analytic pressure datum P0 remains; no higher P0 rows or changed gauge are introduced by this adapter.

Independent mathematical validation uses a polynomial fifth-order axial Psi input with a nonzero mean_Z5, nonunit epsilon/hb, true transported mean and direct physical-R differential operator. Higher U/mean/Q rows are compared symbolically. Both complete signed pressure integrals through Z3 are compared exactly at Z=-2/3, 0 and 3/5: 24 integral rows, retaining the full radial domain and measures. This validation supports the function construction, not original-scale numerical containment. In particular, some radial/axial combinations exceed the existing numerical core oracle's total-order-five interface; their new numerical bounds remain open.

Artifacts:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_analytic_core_lower_pressure_C3_functions.py`
- Same stem `.json.gz`, `_check.py`, `_check.json`.
- `scripts/run_original_n1_analytic_core_lower_pressure_C3_functions.py` inspects accepted source archives without rebuilding.

Completed within this scope:

- [x] Bind required core Z4/Z5 to the same unique analytic twenty-term nonlinear solution and Green inverse.
- [x] Preserve original radial component definitions and lower-order source handles.
- [x] Extend anchored gauge/F0 with the exact moving-endpoint FTC.
- [x] Extend actual collar ratio/Psi/mean histories, U and mean through C5.
- [x] Extend Q and its two radial derivatives through C4 with leading weight 1-delta.
- [x] Extend direct signed lower N1p and m13/m15 integral functions through C3.
- [x] Retain exact binary sc, nonzero hb, physical measures, incoming history and m15 sign.
- [x] Save scoped source acceptance and an archive inspection entry point.

An upstream gap was identified for the active/compact branch: p2 is exported as `full_signed_p2_axial4` or `exact_full_p2_axial4`, with exactly five Taylor entries 0..4. The underlying `current_original_Rm_patch_mixed4_cells.py` explicitly labels `raw_radial_Z_order4_only_no_selected_Z5=True`. Its internal `Qpad` zero supports product bookkeeping and is not a genuine fifth derivative. The switch/reference-restore/long-reshape producers keep this order. A C5 source adapter must first construct the missing radial/inertial fifth source row from actual formulas, including any additional mean/profile/P0 derivatives required by the dependency chain.

Next tasks:

1. [ ] Add genuine higher radial Q/full inertial source Taylor construction. Export new versioned axial5 keys, preserve existing axial4 outputs, and require six genuine entries. Never expose internal Qpad zero as row 5.
2. [ ] Determine and provide every extra axial mean/profile/P0 row needed by the radial/inertial fifth derivative. Extend original defining kernels and source equations; do not differentiate numerical caps.
3. [ ] Bind switch, reference-restore and long-reshape p2_Z5 to those actual source functions with the same family, context, basis, ledger, P0 and selected N. Apply factorial(5)=120 once.
4. [ ] Extend actual active loop/inverse/A/B through C5 with original flat predicates and positive inverse Jacobian. Then extend finite-N density and five-rate Duhamel histories with all original memories.
5. [ ] Normalize genuine complete n0 C4/C5 target functions with actual amplitude derivatives and the correlated divided row. Add same-source quantitative D4/D5 target bounds.
6. [ ] Extend the same compact implicit-control branch with the original Jacobian and selected N:

```text
J=B+DQ(h*)/N
J*h4=-d4-[4*D2Q(h1,h3)+3*D2Q(h2,h2)]/N
J*h5=-d5-[5*D2Q(h1,h4)+10*D2Q(h2,h3)]/N
```

The repository's `quadratic_action` is D2Q, not half of D2Q. Derive quantitative control and Picard tails before claiming numerical C4/C5 control bounds.

7. [ ] Recover true active/compact U/mean C5 and Q C4, then complete original pre-Rp pressure C3, retaining all charts, viscous endpoints and quiet histories.
8. [ ] Compose complete five signed defect C3 functions, retaining current C1 handles and each inner/lower/pre-Rp/pulse contribution exactly once.
9. [ ] Extend selected two-axial/three-angular coefficient C3 with actual Astar derivatives, finite mu and correct physical units.
10. [ ] Recover repaired mean_ZZ and interior Q/V Z1 with n1 weight 1+delta, then required higher recovery rows.
11. [ ] Transfer new core/lower/pulse numerical derivative containment and required mixed/radial C3/C5 norms. Mathematical fixtures and analytic definitions do not supply these bounds.
12. [ ] Bind total own pressure radial derivative and every required inner/band/quiet/pulse join; establish global n1 regularity and quantitative closure.
13. [ ] Complete same-field preheat/finite-energy/exact heat matching, admissible stress and flat remainder.
14. [ ] Implement true n>=2 equations, n-dependent recovery and independent moment repair on the common inner interval.
15. [ ] Bound remainders/smooth summation and implement oscillatory/mean corrections with averaged quadratic stress cancellation.
16. [ ] Independently validate corrected Cartesian NS Linf and volume-L2 residuals below 1e-3 with restricted forcing.

The new analytic core/lower and earlier pulse function results do not establish numerical C3 accuracy, complete pre-Rp/defect C3, selected coefficient C3, mixed norms, global n1 closure, true scale recursion, a blow-up theorem or the corrected NS solution.
