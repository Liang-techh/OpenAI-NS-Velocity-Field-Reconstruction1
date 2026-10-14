# True modified Q C4, total radial C5 and signed pressure C3, 2026-10-14

Ten true modified source charts from original rminus to Rm now have complete Q/full signed inertial C4 functions, total E_y/V_y/mean_y C5 functions, physical U/U_y/mean C5 and Q/Q_y C4 functions, and signed m13/m15 pressure component C3 functions. The actual fast phase, finite intervals, independent pressure source, nonzero correction/background memories and both physical boundary terms remain. This is a function-construction result over ten charts. The actual Rm patch/six outer charts, complete pre-Rp C3 pressure, direct Q_yy/N1p, quantitative mixed bounds, global n1 and genuine recursion remain open. Full goal **ACTIVE / INCOMPLETE**.

Inspect `python scripts/run_original_modified_Q4_pressure_C3_functions.py` without rebuilding ancestors. On the original host reuse session55770/cache26, `public_loop5`, `modified4_owner`, `modified4_receipt`, `public_modified4` and earlier core/lower/pulse/terminal owners. The [preceding ten-chart loop C5 milestone](CURRENT_ORIGINAL_PRESCRIBED_LOOP_C5_FUNCTIONS_2026_10_14.md) remains valid.

Full modified generic recovery now uses genuine modified field/history/pressure rows0..5 and outputs rows0..4. It retains all original signed linear/quadratic/meridional sectors and the independent P0. In particular:

```text
d=1-Z^2, L=1-delta*Z^2
transport=(1-delta)*Z*m+d*m_Z
Q=(2*Z*V-transport)/L
Itheta=R*(theta_linear+Pstar*theta_quadratic)
Iz=R*(axial_linear+Pstar*axial_quadratic)
p2=Iz/E_N
```

The physical radius and Pstar quadratic weight occur once. Missing fifth mean/pressure input rows are rejected; they are not padded. These new modified generic outputs do not replace the prescribed inputs used to define the original loop.

Native differentiation of actual function bodies supplies the previously missing slow A_y/B_y C5. It applies the complete chain rule to function substitutions and at-inverse compositions, includes both moving endpoints of a definite integral, and differentiates the original monotone inverse with its positive lambda. The fast phase is explicitly held fixed in this operator. Accepted background E/V and all normalized leading moment ODE rows provide exact native derivative seeds, each multiplied by the original log-radius Jacobian. Analytic core derivative handles retain the unique twenty-term source/family and existing radial orders0..2; numerical higher-core containment is not asserted.

The total physical radial derivative then restores the original fast terms:

```text
E_N_y=exp(A/N)*(E_y+E*(A_y_slow/N+A_phi))
V_N_y=V_y+B_y_slow/N+B_phi
modified_mean_y=V_N-modified_mean
```

Here slow y is native differentiation divided by its original Jacobian, and physical y=logR. The same exact common N/eta/d_star/global phase and original lazy flat/collar extensions remain. Since the phase is N times log radius, its N factor cancels the 1/N in the displayed fast A_phi/B_phi terms. The new source is a derivative of the periodic composed loop functions; no standalone smoothness claim is made for fractional_part.

Physical U=Pstar*V_N and m=Pstar*modified_mean are recovered once. The leading auxiliary functions are:

```text
Q=(2*Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z)/L
Q_y=(2*Z*U_y-(1-delta)*Z*(U-m)-(1-Z^2)*(U_Z-m_Z))/L
```

This is the leading1-delta recovery. Positive-order n1 own-field recovery retains its separate1+delta coefficient; it is not changed here.

The exact reduced pressure expression is built through Z3 from Q/Q_y C4 and U C5:

```text
S=(Q+(1-delta)*Z*Q_Z/2+Q_y)/L
  +Q*(Q/2+Q_y)
  +U*(-2*Z*Q+(1-Z^2)*Q_Z-2*Z*Q_y)/L
N1p=-S/2+(Q_y+Q_yy)/R
```

The last line states the original identity used for integration by parts; this owner does not install a direct Q_yy/N1p function adapter. It constructs the equivalent original finite pressure integrals:

```text
S0=int_rminus^R S dR
S1=int_rminus^R R*S dR
m13=-S0/2+(Q+Q_y)(R)-(Q+Q_y)(rminus)
m15=+S1/2-R*Q_y(R)+rminus*Q_y(rminus)
```

Each chart carries the complete preceding signed S0/S1 functions and integrates to the actual query endpoint. Measures are R*Jacobian and R^2*Jacobian once. The lower boundary remains the original rminus at microscopic coordinate sc/2; the upper boundary is the actual chart query, reaching Rm. The m15 minus endpoint sign occurs once. Existing lower-collar, pulse and full-preRp C1 components remain valid; the new ten-chart component does not replace them or stand in for the seven missing charts.

Scoped source and independent calculus checks pass: actual background/loop/history identity and derivative seeds, unique native binders, no fast phase used as a substitution binder, moving-limit and lexical chain rules, an independently explicit implicit inverse, full total-y phase-chain calculus through C5, signed generic C4 with nonunit scales, and exact reduced-pressure equivalence to direct N1p integration through C3. Mathematical fixtures do not certify original-scale numerical accuracy. The public checked archive restoration and inspection command pass.

Completed:

- [x] Recover actual modified Q/full signed theta/axial inertial/p2 C4 from genuine modified C5 fields and histories.
- [x] Differentiate actual loop function bodies in native coordinates, with the phase held and original core/denominator/flat sources retained.
- [x] Include moving integral limits and the implicit inverse/at-inverse chain rules.
- [x] Restore total E_y/V_y/mean_y C5 with original A_phi/B_phi fast terms and Jacobian once.
- [x] Recover physical U/U_y/m/m_y C5 and leading Q/Q_y C4 with Pstar once.
- [x] Construct reduced S C3 and the two complete signed finite S0/S1 integral histories through all ten charts.
- [x] Keep original rminus/Rm boundary functions and expose signed m13/m15 pressure C3 endpoints.
- [x] Save scoped checks, source-pinned archive, public checked restoration and an inspection command.

Next tasks:

1. [x] Build versioned actual Rm-patch background C6 and generic C5 sources from the true Rm endpoint, original leading patch controls and support functions. Completed within function scope in [CURRENT_ORIGINAL_TRUE_RM_PATCH_C6_FUNCTIONS_2026_10_14.md](CURRENT_ORIGINAL_TRUE_RM_PATCH_C6_FUNCTIONS_2026_10_14.md). Old axial4 outputs and P0 remain; background and finite-N correction histories remain separate. New C6 numerical norms and patch correction transport remain open.
2. [ ] Extend the actual Rm patch plus Rh_reference/O2_slope/O2_axial/O2_buffer/O3_transition/O3_power sources, loops, finite-N densities and complete correction histories to the required genuine orders. Preserve one original N/global phase, the O2_axial exponential physical Jacobian and quiet O3 incoming memory.
3. [ ] Carry total U_y C5, Q/Q_y C4 and reduced signed pressure C3 through those seven charts to Rp; join with the retained actual lower-collar contribution at rminus. Preserve every signed source term, finite endpoint and physical Pstar/R/R^2 factor. Only then claim complete preRp C3 pressure.
4. [ ] Derive genuine second total radial derivatives and direct Q_yy/N1p functions where direct diagnostics or the required mixed/radial norms consume them. Retain slow-y, phase-y and mixed terms; never pad a missing radial source or treat U_y as U_yy.
5. [ ] Transfer quantitative background C6/C7, loop/modified/history C5 and pressure C3 enclosures with complete core/density/body/tail errors, inverse conditioning and original physical log bases. Use the correlated positive Poisson denominator; function identities alone do not close numerical accuracy.
6. [ ] Complete all17 normalized C4/C5 target functions and D4/D5 source bounds. Differentiate the same original compact implicit-control branch with its J/N and quadratic Hessian factors4,3 and5,10; admit a common N for all required higher derivative/stress/heat conditions.
7. [ ] Compose complete five n1 C3 defects and actual reserved inverse derivatives, including Astar derivatives and nonzero-mu divided differences. Recover repaired own mean_ZZ/Q/V Z1 with1+delta, retaining the existing C1 terminal zero functions.
8. [ ] Prove required mixed/radial norms and global joins; complete compatible preheat/finite-energy/exact heat matching, admissible stress-cone margins and flat remainder.
9. [ ] Implement genuine n>=2 coefficient equations and independent repairs, finite remainders/smooth sum, oscillatory averaged-stress cancellation and independent corrected Cartesian NS Linf/volume-L2 validation below1e-3.

Artifacts: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py`, same stem `.json.gz`, `_check.py`, `_check.json`, and the inspection script above. Full original goal and global scale recursion remain unfinished.
