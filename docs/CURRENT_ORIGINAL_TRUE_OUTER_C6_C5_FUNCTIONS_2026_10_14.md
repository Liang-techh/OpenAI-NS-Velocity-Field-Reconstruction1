# Six original outer C6/C5 sources and pressure C3 to Rc, 2026-10-14

The six remaining prescribed charts now have actual C6 fields and full signed generic C5. Together with the preceding ten charts and Rm patch, all17 prescribed source charts supply the required function orders. Actual selected-N loops, densities, nonzero-inlet correction histories and modified fields C5 continue through all six outer charts. Total radial C5, modified signed generic C4 and signed pressure C3 now reach **Rc**. The compact Rc-to-2Rc repair, quiet 2Rc-to-Rp continuation and complete preRp C3 remain subsequent work. Quantitative high-order/mixed bounds, normalized targets/compact controls, global n1, heat/stress, genuine scale recursion and corrected NS remain open. Full goal **ACTIVE / INCOMPLETE**.

This is the current entry point. Older queue paragraphs referring to these six charts as continuing directly to Rp are superseded: their actual endpoint is Rc. Do not skip the compact repair.

Inspect `python scripts/run_original_true_outer_C6_C5_functions.py` without rebuilding ancestors. On the original host reuse session55770/cache26, `public_Rm_finite5`, `public_Rm_patch6`, `outer65_owner`, `outer65_receipt`, `public_outer65` and earlier owners. On another host inspect the pinned archive and receipt; an in-memory Python session is host-local. The preceding [finite-N Rm patch milestone](CURRENT_ORIGINAL_TRUE_RM_PATCH_FINITE_N_C5_2026_10_14.md) remains valid.

Actual leading Rh fields come from `public_Rm_patch6.endpoints['Rh']`. The actual selected-N Rh correction and cumulative pressure integrals come separately from `public_Rm_finite5.Rh_correction` and `public_Rm_finite5.Rh['integrals']`. Every chart receives the preceding chart's complete endpoint by object identity. Leading repair N0, the current selected repair-only N and legacy N257 densities remain distinct. No old all-N normalized correction or numerical interval endpoint defines a field value.

The same original unique leading branch, exact full-weight map and complete-support Rh function identity are explicitly bound to the predecessor C6 receipt. The normalized Rh functions are

```text
E=exp(-1/2)/(1+Z^2), V=m=4Z/Pstar
h=5E/8, k=h*V, e=V^2-5E^2/12, p=5E^2/2
pressure=P0+p
```

The independent P0 remains. The original identity is an identity of smooth functions; the new C6 derivative rows come from the preceding actual implicit-function construction. An old C3 range or C1 correction is not promoted to a C6/C5 source. Equal E/V and equal incoming histories identify the propagated leading histories by uniqueness of the original normalized linear moment ODEs:

```text
m_y=V-m
h_y=E-3h/2
k_y=E*V-3k/2
e_y=V^2-E^2/2-e
p_y=E^2/2
```

The new representation uses integral kernels rather than dividing by tiny mu in closed formulas. Original profile recipes, all ordinary amplitude/control derivatives and the actual physical radial maps remain:

| Chart | Native domain | log(R/Rm) | dlogR/dx |
| --- | --- | --- | --- |
| Rh_reference | [-5,0] | 6+x | 1 |
| O2_slope | [0,1] | 6+x | 1 |
| O2_axial | [0,1] | 6+exp(40x) | 40exp(40x) |
| O2_buffer | [0,11] | 6+exp(40)+x | 1 |
| O3_transition | [0,1] | 17+exp(40)+x | 1 |
| O3_power | [0,2] | 18+exp(40)+x | 1 |

Thus `Rc=Rm*exp(exp(40)+20)`. The original canonical recipe is `logmu=log(1/1000)-4*(exp(40)+11)` and `Tw=-60logmu`; the original eta guard is retained. Toy parameters used in calculus checks do not replace these source parameters.

For leading histories and actual finite-N corrections, write y(x)=log(R/Rm), J=y'(x). Every rate uses

```text
H_j(x)=exp(-rate_j*(y(x)-y(a)))*H_j(a)
       +int_a^x exp(-rate_j*(y(x)-y(t)))*density_j(t)*J(t) dt
rates: m=1, h=k=3/2, e=1, p=0
```

The exponential O2 map is integrated exactly; its Jacobian occurs once. O3 quiet local density does not erase inherited correction: positive-rate memory decays and pressure rate-zero memory persists. All loops use the original global phase `frac(N*(log(R/Ra)-hbB*sc/2))` and the original full signed p2. Native differentiation holds phase fixed, then restores actual A_phi/B_phi terms. Standalone fractional_part smoothness is not claimed.

The pressure source retains the preceding rminus-to-Rh cumulative integrals and original rminus boundary. With S from actual physical modified U/mean and Q/Q_y,

```text
S0(x)=S0(a)+int_a^x R(t)*J(t)*S(t) dt
S1(x)=S1(a)+int_a^x R(t)^2*J(t)*S(t) dt
m13=-S0/2+(Q+Q_y)(R)-(Q+Q_y)(rminus)
m15=+S1/2-R*Q_y(R)+rminus*Q_y(rminus)
```

This is the rminus-to-Rc component, not the complete preRp pressure. Independent scoped checks cover source owners/rows/binders, original Rh theorem and normalized moment ODE, all maps/seams,168 actual profile/first-radial rows,30 quiet incoming correction rows,65 exponential transport rows and8 direct-N1p pressure rows. These checks establish formulas and provenance, not original-scale high-order norm admission.

Completed within function scope:

- [x] Rh_reference prescribed C6/generic C5 and actual selected-N correction C5.
- [x] O2_slope, O2_axial and O2_buffer C6/generic C5 and loop/density/history C5.
- [x] O3_transition and O3_power through their original terminal Rc.
- [x] Preserve all leading/correction/pressure inlet memories, original N/global phase and P0.
- [x] Retain the original exponential O2 Jacobian and canonical mu/eta guard.
- [x] Bind the original Rh leading function identity, preceding true C6 branch and normalized history ODE.
- [x] Construct all17 modified total radial C5/generic C4 and physical Q/Q_y C4 functions.
- [x] Carry signed rminus-to-Rc pressure C3 with every physical measure and boundary.
- [x] Save source-pinned archive/receipt and inspection command.

Execute the next tasks in dependency order. Mark a task complete only with its actual artifact, command and scoped evidence; retain open numerical/global gates.

Tasks1..9 are complete within the [actual compact C5/2Rc milestone](CURRENT_ORIGINAL_TRUE_COMPACT_C5_FUNCTIONS_2026_10_14.md). Numerical norms/global frequency and full preRp remain open; task10 is the next dependency.

1. [x] Build a true Rc compact-repair inlet adapter from `public_outer65.Rc`, preserving actual leading fields, five nonzero selected-N corrections, P0 and S0/S1. Bind leading seed equality by the original source theorem. Do not consume old C3 range records as field functions.
2. [x] Use `current_original_outer_to_repair_C3_bridge.py` (`compare_functions`, `correction_inlet_functions`) as the original convention reference. Extend the five-history interface from Z0..3 to genuine Z0..5. Account for one division by selected N only; current actual correction already includes its N factor.
3. [x] Produce genuine all17 normalized target C4/C5 and D4/D5 source functions needed by the compact five-control solve. Keep original signs, physical angular/axial scales and inherited pressure datum.
4. [x] Differentiate the original compact same-J implicit solve through C5, including all ordinary amplitude and inverse cross terms. Preserve original Hessian factors4,3 and5,10; a padded row or interval derivative is not an input.
5. [x] Lift `current_original_C3_repair_band.py` (`profiles_C3`, `signed_density_C3`, `leading_C3`, `build`) into a separate checked C5 provider on x=R/Rc in[1,2], using the actual Rc inlet and genuine control rows.
6. [x] Carry compact E/V and all five signed densities C5 with every linear/quadratic cross term. Preserve the original profile supports and flat endpoint behavior.
7. [x] Carry complete compact correction histories C5 with log decay, dt/t measure and pressure rate-zero memory. Expose an actual 2Rc endpoint; do not zero incoming memory in quiet subregions.
8. [x] Recover compact modified total radial U/U_y/mean C5 and Q/Q_y C4 with the actual native chain rule and fast terms where present. Preserve leading1-delta and positive-order own1+delta conventions separately.
9. [x] Carry signed compact pressure C3 from Rc to2Rc, inheriting this S0/S1 and the original rminus boundary, with exactly one R/R^2 measure.
10. [ ] Feed the genuine compact2Rc endpoint into the original `current_original_C3_power_to_Rp.py` power convention. Use s=log(R/(2Rc)), lengthTw-2-log2 and Rp=Rc*exp(Tw-2), retaining finite-mu memory and independent P0.
11. [ ] Construct the required quiet2Rc-to-Rp field/history C5 and pressure C3 projections. Differentiate actual functions; do not infer pressure zeros from small relative correction or flat local density.
12. [ ] Establish actual joins at Rc,2Rc andRp, including required radial/axial rows and all pressure endpoint terms. Join the retained lower-collar C3 at rminus; only then expose complete preRp C3.
13. [ ] Construct genuine second total radial and required mixed rows/direct Q_yy/N1p for diagnostics and norm admission. Current first total-radial functions do not close this requirement.
14. [ ] Transfer quantitative actual leading/control C6, comparison C7, loop/modified/history C5 and pressure C3 bounds, including all analytic core/density/body/tail errors, inverse conditioning, phase and original physical/log bases.
15. [ ] Admit one common N for all higher-derivative, compact, stress and heat conditions. The existing repair-only integer is not already a global frequency certificate.
16. [ ] Compose complete five n1 C3 defects from inner/lower/preRp/pulse components, with actual Astar derivatives and nonzero-mu divided differences. Retain all original signed sectors and finite endpoints.
17. [ ] Differentiate the actual selected reserved inverse through required Z2/Z3. Use complete source derivatives; preserve the proved C1 terminal five zeros and their pressure-tail consequences.
18. [ ] Recover repaired own mean_ZZ/Q/V Z1 and required higher/mixed norms using1+delta, then establish complete global n1 joins and finite energy.
19. [ ] Complete compatible preheat pressure, exact heat exterior, admissible stress-cone margins and flat remainder with actual derivative/tail evidence.
20. [ ] Implement genuine n>=2 coefficient equations and independent five-moment repairs. Coefficients must depend on earlier accepted orders; duplicated n1 or a shrinking visualization is not scale recursion.
21. [ ] Establish finite remainders and a smooth summed field, original two-family oscillatory corrections and averaged quadratic stress cancellation.
22. [ ] Evaluate corrected Cartesian [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] with restricted forcing; independently verify NS Linf and volume-L2 residual below1e-3. Keep the final goal open until this succeeds.

Artifacts: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py`, same stem `.json.gz`, `_check.py`, `_check.json`, and the inspection command above. No new PR, automation or additional agent is needed to resume this queue.
