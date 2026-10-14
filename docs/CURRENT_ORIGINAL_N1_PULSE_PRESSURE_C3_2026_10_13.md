# Actual native/Gamma C5 and pulse pressure C3, 2026-10-13

Status: genuine native/Gamma C5 functions now feed leading main/gap/end pulse velocity and own mean C5, auxiliary Q C4, and signed n1 pulse pressure C3. Existing lower-order source handles, physical units, finite endpoints and nonzero histories are retained. This completes one higher-derivative dependency after the [actual C1 terminal repair](CURRENT_ORIGINAL_N1_TERMINAL_INTEGRAL_IDENTITIES_2026_10_13.md). The full goal remains **ACTIVE / INCOMPLETE**: complete five-defect/coefficient C3, quantitative derivative containment, mixed/radial norms, global n1 joins and true n>=2 recursion remain open.

The exact checked native source is used as a whole function:

```text
Q0=1+Z^2
u=U/Q0
m_j=C_j*Z*Q0
energy=C0+C_E*Z^2*Q0^2
X=constant
```

U, C_j, C0, C_E and X are Z-independent source scalars. The new native derivatives are ordinary derivatives. Their comparison with the actual six-row native Taylor callback applies factorial(k) once. All existing C3 FunctionRefs are preserved. The cubic m_j have proved zero fourth/fifth derivatives; energy retains its fourth and sixth polynomial powers. Its higher rows are `C_E*(48+360*Z^2)` and `720*C_E*Z`. C_E is recovered from the actual energy_ZZ function divided by `2*(1+12*Z^2+15*Z^4)`, which is at least 2 on the real domain. The quotient's equality to the source scalar is checked before it is treated as constant. The actual X source and its earlier derivatives are checked before adding higher zeros.

The Gamma extension keeps the same full positive improper kernels, actual nonzero S and finite offset t:

```text
a=delta/2, T=exp(-t), xi=2*(1-Z^2)*S*T
M_k = integral[0,infinity] exp(-v)*v^(a+k)/(1+xi*v)^(a+k) dv / Gamma(1+a)
M_k,Z = 4*(a+k)*Z*S*T*M_(k+1)
0 <= M_k <= Gamma(1+a+k)/Gamma(1+a), k=1..5

Dhat_Z4 = -48*(a+1)*S*T^2*M2
           -384*(a+1)*(a+2)*S^2*T^3*Z^2*M3
           -256*(a+1)*(a+2)*(a+3)*S^3*T^4*Z^4*M4
Dhat_Z5 = -960*(a+1)*(a+2)*S^2*T^3*Z*M3
           -2560*(a+1)*(a+2)*(a+3)*S^3*T^4*Z^3*M4
           -1024*(a+1)*(a+2)*(a+3)*(a+4)*S^4*T^5*Z^5*M5
```

The original M1/M2/M3 handles remain. M4/M5 use the same normalization and endpoints. There is no division by 1-Z^2 at Z=+/-1. The kernels have an integrable Gamma bound at both v endpoints and exponential decay in t. Exact kernel differentiation and the Dhat recurrence are checked independently. Finite mathematical fixtures differentiate the original Dhat integral through orders four and five and exercise proper-integral/product/quotient calculus through order five. These fixtures do not certify the original astronomical-scale numerical field.

Selected ap/C1/C2 functions and complete main/gap/end histories are differentiated as actual source bodies. Leading U/U_y/m/m_y/Utheta now have C5 rows, while Q/Q_y have C4. The leading recovery remains:

```text
Q = [2*Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z]/(1-delta*Z^2)
```

Its 1-delta weight is distinct from the repaired n1 1+delta weight. The inactive gap retains its cumulative mean. Absolute leading swirl pressure, with the original P0/Mp, remains separate from the reduced n1 forcing pressure.

The thirty signed pulse pressure integrals retain physical measures `R*ds`, `R^2*ds`, `(R/mu)*dxi` and `(R^2/mu)*dxi`. Their Z2/Z3 rows differentiate the actual integrands at fixed original endpoints. Finite viscous endpoint functions are extended too. Eight pressure C3 functions are exposed: m13/m15 in each of main/gap/end, and their full Rp-to-Rv totals with the original outer endpoint telescope. The original m15 sign is applied once. Pstar is not reapplied.

The checked public owner requires the exact source archive, proof receipt, source hashes, native owner/family and function bindings. All referenced graph nodes belong to an exact pinned prefix. A shared live graph may append unrelated nodes after that prefix; those nodes are excluded from this result, and cannot replace its function handles or alter protected nodes.

Artifacts:

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.py`
- Same stem `.json.gz`, `_check.py`, `_check.json`.
- `scripts/run_original_n1_pulse_pressure_C3_functions.py` inspects accepted archives without rebuilding ancestors.

Completed within this scope:

- [x] Extend genuine native u/m1/m2/energy/X ordinary Z4/Z5 from the same source identities.
- [x] Bind reconstructed U/C_E and the constant X to the original source and Taylor factorial convention.
- [x] Construct full Gamma M4/M5 and exact Dhat Z4/Z5, preserving the earlier source handles.
- [x] Propagate genuine selected pulse ap/C1/C2 through C5.
- [x] Extend complete leading main/gap/end U/U_y/Utheta/m/m_y through C5, retaining nonzero gap histories.
- [x] Extend leading Q/Q_y through C4 with the original 1-delta weight.
- [x] Extend the signed pulse pressure integrals, endpoint functions and eight pressure functions through C3.
- [x] Save source-pinned scoped acceptance and a lightweight archive entry point.

Next tasks, in dependency order:

1. [ ] Expose genuine leading active-chart U/mean Z4/Z5 FunctionRefs from the actual unprojected source rows. Bind context, family, native coordinates, AST/source hashes and factorial conversion; magnitude axial5 rows alone are insufficient.
2. [ ] Extend the compact n0 control system to C4/C5 with actual target derivatives and the same invertible implicit Jacobian/selected N. Derive the fourth/fifth differentiated equations and quantitative tails; do not invent opaque derivative leaves.
3. [ ] Recover active/compact Q/Q_y through C4 from those true U/mean C5 functions.
4. [ ] Extend complete pre-Rp m13/m15 through C3, preserving original physical measures, pressure datum, selected compact controls and active/compact/quiet partition.
5. [ ] Extend the unique analytic inner/lower pressure functions through C3 using their actual core/ratio/Psi/own-mean kernels and fixed domain.
6. [ ] Compose complete five signed defect C3 functions with existing velocity moments and the new pressure segments. Preserve current C1 handles and consume each pressure contribution once.
7. [ ] Extend the selected two-axial/three-angular inverse through C3 with finite mu, correlated axial divided differences and genuine Astar derivatives. Use independent physical coefficient units.
8. [ ] Recover reserved-band own mean_ZZ and interior Q/V Z1 with the repaired n1 1+delta weight; then add the required higher rows from actual selected coefficients.
9. [ ] Transfer numerical derivative containment to the new pulse and other source functions. Retain directed errors and source provenance; mathematical fixtures do not supply this containment.
10. [ ] Bound required mixed/radial C3/C5 repair norms through physical compact bumps, actual radius factors and original nonzero mu.
11. [ ] Bind total own pressure's radial derivative and all required inner/band/quiet/pulse joins. Retain the original gauge and fixed leading N1p.
12. [ ] Establish complete global n1 regularity, required norms and closure. Existing terminal C1 identities remain useful but do not imply this gate.
13. [ ] Complete same-field preheat/finite-energy/exact heat matching, admissible stress and flat remainder.
14. [ ] Implement true n>=2 equations, n-dependent recovery and independent moment repair on the common inner interval.
15. [ ] Bound finite-order remainders and smooth summation; implement oscillatory/mean corrections and averaged quadratic stress cancellation.
16. [ ] Independently validate corrected Cartesian NS Linf and volume-L2 residuals below 1e-3 with restricted forcing.

This milestone supplies actual derivative functions on the pulse branch. It does not certify original numerical C3 accuracy, complete five-defect C3, reserved coefficient C3, mixed norms, global n1 closure, genuine scale recursion, a blow-up theorem or the final corrected NS solution.
