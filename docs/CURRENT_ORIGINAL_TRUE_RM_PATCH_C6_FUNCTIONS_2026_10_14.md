# Original leading Rm patch C6 functions, 2026-10-14

The [finite-N continuation milestone](CURRENT_ORIGINAL_TRUE_RM_PATCH_FINITE_N_C5_2026_10_14.md) now completes patch correction/modified C5, total radial functions and signed pressure C3 to Rh. The leading-only construction described below remains valid; its first six follow-up tasks are completed within function scope.

The true prescribed Rm inlet now supplies the original five leading implicit control functions through Z6, actual Rm-to-Rh background E/V/m/h/k/e/p/P0/pressure and E_y/V_y C6, and full signed Q/inertial/p2 C5. The map uses exact original finite bump integrals and the same uniquely admitted leading branch. Exact source normalization, raw-to-normalized map/Jacobian conjugacy and the inherited nonzero-J theorem are checked. This extends the leading background function construction. Finite-N patch corrections, six outer charts, complete preRp pressure C3 and new original-scale numerical containment remain open. Full goal **ACTIVE / INCOMPLETE**.

Inspect `python scripts/run_original_true_Rm_patch_C6_functions.py` without rebuilding ancestors. On the original host reuse session55770/cache26, `public_modified4`, `Rm_patch6_owner`, `Rm_patch6_receipt`, `public_Rm_patch6`, and the checked core/lower/pulse/terminal owners. The [preceding ten-chart modified Q/pressure milestone](CURRENT_ORIGINAL_MODIFIED_Q4_PRESSURE_C3_FUNCTIONS_2026_10_14.md) remains valid.

The true input is `public_modified4.background.endpoints['Rm']`, with seven genuine ordinary axial rows per prescribed field. The leading patch uses the prescribed background inlet; its finite-N correction history is a separate nonzero incoming object, `public_loop5.Rm_correction`. Neither the older C5 interval adapter nor a midpoint defines a new field value.

The actual parent normalization is

```text
W=sqrt(2R)*F0/Pstar, FF=R*F0^2/Pstar^2=W^2/2
E=W*phi, V=U/Pstar, m=M/Pstar
h=W*H/2, k=W*K/(2Pstar)
e=A/Pstar^2-FF*B, p=FF*C, pressure=P0+p
```

At the original offset -6, E=Utheta/Pstar=am=exp(-3/5)/(1+Z^2). W retains its F0/radius dependence. With Theta=H/(2phi), Theta_mixed=K/(2phi), and Am=Pstar*am, the five input defects are exactly

```text
d0=M-4Z
d1=Theta_mixed-4Z*Theta
d2=Theta-5/8
d3=(A-8Z*M+16Z^2)/Am^2-(B/phi^2-5/6)/2
d4=(C/phi^2-5)/2
```

The pressure datum P0 is retained in the total field and excluded as a directly added datum in these five defects. Existing source moments may still depend on P0 through their prescribed dynamics. The checker proves the defining normalization identities on the actual parent nodes with arbitrary Z/raw moments/positive Pstar,Rm,F0, then checks all seven ordinary rows against independent nonconstant raw inputs. The swirl/energy and pressure halves occur once.

For centers 5/4, 3/2, 7/4 and radius1/40, gamma_i(x)=40*beta(40(x-center_i))/int_-1^1 beta, with beta(u)=exp(-1/(1-u^2)) on |u|<1 and exact0 outside. All full/partial weights are actual finite integrals of x^p gamma_i(x)^m. The checker compares their integrands, bounds, normalization and dx measure to the original directed integral definitions, pins that proof, and checks the exact unit-mass entries against its certified records. Beta derivatives use their original flat extension at support boundaries.

The original leading map is F(h,Z)=Lh+Q(h,alpha)+d, with alpha=Am^-2. Its values are the same unique original common-box solution of F=0; an enclosure or a residual interval containing zero does not define them. The exact branch theorem is bound by identical fixed N0, source family, prescribed source/initial-data chain, source normalization and exact integral identities. The normalization is encoded as function graphs:

```text
S=diag(t,t,a,a,a), a=10^10*t^2/Pstar^2
G(u)=S^-1 F(Su)
J_G=S^-1 J_F S
alpha_normalized=alpha*t^2/a=exp(6/5)*(1+Z^2)^2/10^10
```

L commutes with S because its original2+3 blocks do not mix. The normalized cross/swirl/pressure weights acquire a; axial-square weights remain unchanged. Symbolic checks cover all five map and25 Jacobian entries with arbitrary controls, defects, weights and positive units. The original strict whole-Z common-box Neumann bound applies to this same normalized map, so J_F is invertible and det(J_F)^2>0. The exported adjugate inverse uses that strictly positive squared determinant.

The transfer also compares the predecessor archive's actual source-binding dictionary by exact equality, binds each full exact weight to its individual original directed interval, and verifies that all four old angular-weight covers contain the actual exact angular unit. Its coefficient is checked as an exact rational10^10*t^2; the positive source scale is Pstar^-2. Directed comparisons in log space retain the original microscopic unit without using an enclosure endpoint as a field value.

For each ordinary derivative order1..6, the known lower control rows and every alpha derivative build the full nonlinear known RHS. The next row solves the same J_F. Only the unknown row is temporarily zero inside this coefficient equation. The independent fixture has nonzero sixth derivatives for all five controls and alpha and recovers all35 ordinary rows exactly; it also checks the full inverse algebra.

At x=R/Rm in[1,e], the original partial residuals d+L_partial*h+Q_partial reconstruct all five cumulative primitives. The patch retains E=am*H, V=(4Z+axial_correction)/Pstar, independent P0, radius=Rm*x and physical log-radius Jacobian1/x. At x1 every bump/partial integral is zero, recovering the actual inlet. Beyond the last support71/40, partial weights equal full weights and the same unique-map identity closes the five residual functions; it is not a new inlet reset. The original Rh targets, including amplitude Z products and physical scales, are checked through C6. Generic recovery consumes genuine C6 field/pressure inputs and supplies all signed C5 sectors.

Completed within this function scope:

- [x] Bind the true prescribed Rm C6 input and original independent P0.
- [x] Construct exact original beta normalization, full/partial integrals and flat first/second radial bump derivatives.
- [x] Prove actual defining normalization identities with arbitrary Z/raw primitives and the original energy/pressure halves.
- [x] Check all seven ordinary defect/amplitude/alpha rows with independent nonconstant raw functions.
- [x] Encode and prove the original raw/normalized five-map and Jacobian conjugacy.
- [x] Bind the original unique whole-Z common-box branch and nonzero-J theorem to this exact map.
- [x] Construct genuine ordinary leading control rows0..6, including all amplitude products.
- [x] Recover actual leading background fields/radial fields C6 and complete generic C5 from Rm to Rh.
- [x] Retain original inlet/terminal matching identities, physical factors and pressure datum.
- [x] Save source-pinned archive/receipt, checked public restoration and inspection command.

Next tasks, in dependency order:

1. [x] Add actual Section11 loop functions on this prescribed patch. Completed within function scope in [CURRENT_ORIGINAL_TRUE_RM_PATCH_FINITE_N_C5_2026_10_14.md](CURRENT_ORIGINAL_TRUE_RM_PATCH_FINITE_N_C5_2026_10_14.md); genuine roots, original N/eta/d_star/global phase and separate leading controls remain.
2. [x] Implement original log-coordinate dt/t transport, with x^-r incoming and exact log suffixes. The variable1/x measure is retained.
3. [x] Carry every nonzero `public_loop5.Rm_correction` row through the patch, separately from background history; pressure rate0 retains its inlet.
4. [x] Recover finite-N E/V, five signed densities, modified moments/pressure C5 and actual Rh endpoints with all cross/quadratic amplitude terms.
5. [x] Supply required native slow derivatives with original gamma first/second rows, Z-only implicit controls, full moment ODEs, moving-integral endpoint rule, x-Jacobian and explicit fast phase terms. Further mixed/second-radial sources and numerical norms remain open.
6. [x] Recover total U_y/mean C5, Q/Q_y C4, reduced S C3 and signed pressure contribution to Rh, with Rm/Rm^2*x measures, nonzero S0/S1 and original rminus boundary. Complete preRp C3 remains open.
7. [ ] Construct Rh_reference prescribed C6/generic C5 with the actual patch endpoint and original interpolation functions. Preserve all six leading background memories and independent P0.
8. [ ] Construct O2_slope C6/generic C5 from its original physical profile and exact endpoint data; retain original radius/phase origins.
9. [ ] Construct O2_axial C6/generic C5 with the original exponential physical Jacobian. Do not replace its logarithmic radial measure by a constant interval width.
10. [ ] Construct O2_buffer C6/generic C5 and its complete nonzero leading memories.
11. [ ] Construct O3_transition C6/generic C5 with original quiet-loop predicate and the full incoming finite-N correction history.
12. [ ] Construct O3_power C6/generic C5 through the original Rp. A quiet density does not imply zero correction memory.
13. [ ] Propagate loop/inverse/A/B C5, finite-N density/correction C5 and modified functions through each of these six charts. Keep one common original N/global phase and original support conditions.
14. [ ] Propagate total radial U_y C5, Q/Q_y C4 and signed reduced pressure C3 through all seven missing charts. Join retained actual lower-collar data at rminus before claiming complete preRp pressure C3.
15. [ ] Add genuine second total radial/mixed derivatives and direct Q_yy/N1p where diagnostics or required norms consume them; U_y does not supply U_yy.
16. [ ] Transfer original-scale leading control/background C6, comparison C7, loop/modified/history C5 and pressure C3 numerical bounds. Include complete actual core/density/body/tail errors, inverse conditioning, amplitude factors and physical log bases. This function receipt does not supply those norms.
17. [ ] Complete all17 normalized C4/C5 targets and source-derived D4/D5 bounds. Retain original compact J/N and Hessian factors4,3 and5,10; admit one N for all higher derivative/stress/heat conditions.
18. [ ] Compose all five complete n1 C3 defects and reserved inverse derivatives, with Astar derivatives and nonzero-mu divided differences. Recover repaired own mean_ZZ/Q/V Z1 with1+delta and retain the existing C1 terminal five zeros.
19. [ ] Complete required mixed/radial norms and global joins, compatible preheat pressure, finite-energy/exact heat matching, admissible stress-cone margins and flat remainder.
20. [ ] Implement genuine n>=2 coefficient equations and independent repairs, finite remainders/smooth sum, oscillatory averaged-stress cancellation and independent corrected Cartesian NS Linf/volume-L2 below1e-3.

Artifacts: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py`, same stem `.json.gz`, `_check.py`, `_check.json`, and the inspection script above. Checked function construction is distinct from numerical norm admission, global n1 closure, true recursion and corrected Navier-Stokes validation.
