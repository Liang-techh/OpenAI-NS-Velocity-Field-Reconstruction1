# Actual current finite-N prefix through O2 axial turnoff and the eleven-unit buffer

Checked source [2289af5f](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/2289af5fec977d4fbb6153f23d5a8036fe3dd652). Full reconstruction **ACTIVE / INCOMPLETE**. Previous actual input: [O2 slope exit](CURRENT_ORIGINAL_O2_SLOPE_FINITE_N_2026_10_09.md).

The same-current actual N257 correction now reaches Rd at Z0,.5 through the original O2 axial turnoff and all eleven buffer units. This closes the next production dependency: the genuine Rd incoming and original background seed are exported separately for O3. These are conservative source-function enclosures; sharp defects, whole-Z functional five-moment closure, cone/global N and n-dependent recursion remain open.

## Source and physical coordinates

The accepted original parameter definition is pinned to Md40, logPstar=exp(40)+11 and the same physical/implicit pressure family. Its checked parameter-frame receipt and entire source hashes are hydrated without constructing the old point owner. Radius is Rm_factor*exp(6+y), with the compact y and formal radius factor kept separate from enormous logCstar.

    axial: phase p in[0,1], y=exp(40*p), B(y)=sigma(1-log(y)/40),
    buffer: tau in[0,11], y=exp(40)+tau, B=0,
    dy=40*exp(40*p)*dp on axial; dy=dtau on buffer,
    Rd=Rref*exp(exp(40)+11).

The actual full source uses ordinary physical-log-radius cutoff derivatives, not phase derivatives. Axial width is exp(40)-1, not1; the buffer adds exactly11. Phase remains globally anchored as frac(N*(logRref+y-logRa-hb*s_c/2)), with the original positive microscopic source origin. A full-period cover encloses phase and its exact zero Z derivative; it is not a selected phase or sharp averaged cancellation.

## Original field and all five retained histories

Let u1 and m1,h1,k1,e1,p1 be the live original background at the accepted slope exit, t=y-1, and K_j(y)=integral_1^y exp(s-y)*B(s)^j ds. The actual formulas are:

    E=u1*exp(-t/2), V=4Z*B,
    m=m1*exp(-t)+4Z*K1,
    h=h1*exp(-3t/2)+u1*(exp(-t/2)-exp(-3t/2)),
    k=k1*exp(-3t/2)+4Z*u1*exp(-t/2)*K1,
    e=e1*exp(-t)+16Z^2/Pstar^2*K2-u1^2*t*exp(-t)/2,
    p=p1+u1^2*(1-exp(-t))/2.

The source is AST-bound to the original axial/buffer recipe. Only pure turnoff_kernels and turnoff_derivatives are reused. Directed endpoint-weighted rectangles retain the nonzero far kernel tail. In the buffer the exact identity K_j(exp(40)+tau)=exp(-tau)*K_j(exp(40)) is used, rather than resetting K to0 or recomputing an unrelated seed. All amplitude/history decays and the physical radius remain formal positive exponentials; no homogeneous incoming history underflows to zero.

The original shear a=1-2Ey/E=2 is exact throughout both charts. The full nonzero axial shear b=2Vy/(Pstar*E), its genuine midplane Z derivative, pressure, meridional terms and all signed inertia are retained. The original kappa-2 is collected as b^2/2 before interval addition/subtraction. Computing (2+b^2/2)-2 numerically could overwhelm the defining tiny positive term. Only this algebraically equivalent assignment is adapted in the accepted generic quotient function; all remaining signed equations and original q cutoff branches are unchanged.

At Z0, b0 and p2=0 but b_Z and p2_Z remain. In the buffer b0,a2 imply q^2=eta/2>0, so zero axial velocity does not make this region quiet. The signed all-u inverse/Poisson primitive bounds include these axis derivatives and all five nonlinear density C0/Z rows. Original analytic P0 stays the same live object and is not a finite-N correction.

## Real finite-N transport, evidence and limitation

Four phase quarter-cells per frame cover the whole axial source, followed by buffer cells[0,1],[1,5],[5,11]. Each integral uses its actual physical logarithmic width, positive own-rate mass enclosure and formal suffix to the chart exit. The real slope-exit incoming is restored from the current accepted report, guarded by source/frame/P0/N/radius. Axial output is the buffer incoming; buffer output is the genuine Rd correction. Pressure rate0 retains memory exactly1. The two source seams follow from exact zero kernel/decay and retained-kernel identities; 42 row overlaps per seam check consistency only.

Artifacts: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N.py, .json.gz, _check.py, _check.json. Producer 60.781s and checker 64.093s, terminal exit0. 1248 exact staged dependency hashes PASS. 14 source cells, 588 common coefficients, 140 signed density rows, 120 boundary/driver rows, 7 midplane cases and 6 positive-q buffer cases. Independent symbolic five-history/yZ/Jacobian/seam identities and nine explicit scalar cutoff/kernel/ordinary-derivative comparisons pass. 14 wrong typed inputs are rejected. Nonzero homogeneous history survives the huge axial width formally. One reused GPT-5.6 Luna/max read-only worker reviewed kernel/history pitfalls; no new agent or old source producer/checker runs.

The four axial cells and full-period phase cover are intentionally broad. They provide executable conservative original-source integration, not a small terminal defect, chosen velocity or whole-Z theorem. Sharper correlation/phase integration remains necessary for actual functional repair.

## Next production tasks

- [x] **W4c.5c.4c AXIAL SOURCE:** actual same-current original Md40 source, all five histories, full signed recovery and physical-coordinate derivatives at N257,Z0,.5.
- [x] **W4c.5c.4d AXIAL DRIVER:** actual physical widths/suffixes and original-rate C0/Z integrals with genuine slope-exit incoming.
- [x] **W4c.5c.4e BUFFER:** eleven-unit original history/finite-N transport, nonzero K tail, positive critical eta q and genuine Rd incoming/background export.
- [ ] **NEXT W4c.5c.4f O3 SOURCE:** consume exported actual_Rd_background_source and actual_current_Rd_correction_C0_Z. Guard same family,P0,N257,Z0,.5 and canonical exact_source_Rd_factor. Hydrate exact mu=.001*exp(-4*(exp40+11)); never select an enclosure endpoint or construct old source owners.
- [ ] **O3 FORMAL CRITICAL SHEAR:** implement E=E_d*exp(-tau/2-mu*J(tau)), V0 with original directed transition_kernels. Retain all old five histories. Express Delta=2mu*sigma(tau) directly before forming a=2+2mu*sigma; tiny mu must not disappear in ordinary addition to2. Preserve q/q_Z and true p2_Z.
- [ ] **O3 TRANSITION DRIVER:** original offset[0,1], physical dy=dtau, same all-u C0/Z signed primitives, real Rd incoming, pressure memory1 and genuine Rw output; exact seam/source checks.
- [ ] **W4c.5c.5 QUIET POWER:** bind actual Tw=-60logmu and Rc=Rw*exp2. Prove original q/A/B flatness on the actual power segments by mu/eta inequalities before claiming local driver0. Transport real incoming memory even on certified quiet segments. Do not use an old N1024 tail packet.
- [ ] **W4c.5c.6 ACTUAL Rc INPUT:** compose the complete same-current rminus-to-Rc route and export source-bound C0/Z functions with exact parameters/radii for independent five-control functional repair.
- [ ] **SHARPNESS / ALL-Z:** refine axial source cells and preserve local Duhamel/phase correlations, including a geometric physical-y atlas where useful. Use stable sharper own-rate masses rather than merely increasing checks. Extend all-Z and higher jets; a handful of axial frames is not functional closure.
- [ ] **RM-W5..W8:** actual five-function terminal repair, axis/analytic pressure/heat joins, admissible stress cone and one global admitted N.
- [ ] **RM-W9..W11:** exact heat/finite energy, stress/flat remainder, actual n-dependent coefficient recursion, mean/pulse stress cancellation, full corrected NS and physical core/material-winding diagnostics.

Mark DONE only with implemented source, bounded report/receipt and source commit. Continue with the actual O3 consumer. Full objective remains active and incomplete.
