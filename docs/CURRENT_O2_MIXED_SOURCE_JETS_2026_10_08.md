Current successor: [CURRENT_O2_REGULAR_MIXED_PHASE_2026_10_08.md](CURRENT_O2_REGULAR_MIXED_PHASE_2026_10_08.md), checked source 1e16263f. The three regular fixed-angle/fixed-phi/full-A-B tasks below are DONE only for the actual original O2 regular axis application. Whole nonzero-Z/signed domain coverage and changed five integrals remain open; follow the successor queue. The full reconstruction goal remains active.

# Original O2 slope whole-domain mixed source and varying parameters

Checked source [7f852778](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/7f852778afe78e633f426d62f58424e6b8a7caac). Predecessor: [CURRENT_REFERENCE_FULL_AXIAL_C1_2026_10_08.md](CURRENT_REFERENCE_FULL_AXIAL_C1_2026_10_08.md). Root implemented and computed; the existing GPT-5.6 Luna/max worker reviewed read-only. No child was spawned. The full original reconstruction goal remains active.

## Concrete advance

Original O2 slope now supplies issued, family-bound C0/y/Z/yZ source frames over the complete original y[0,1], Z[-1,1] domain. All64 continuous radial source cells cover that domain; exact-axis and signed subinterval frames are also supported. Full original E,V,p1,p2 coefficients, pressure prefix/suffix function relation and late pressure errors remain. This extends the reference constant-parameter interface to actual varying a,q and phase normalization nu. **nu here is the phase normalization 1+2q^2, not the physical Navier-Stokes viscosity.**

The defining integrals imply f_y=(1-a)f/2, H_y=f-3H/2, D_y=f^2/2-D and P_y=f^2/2. They differentiate the original coefficient templates and their ordinary Z rows, retaining the native radius derivative and applying L_Z once. The exact pressure datum has zero y derivative; differentiating alpha=P+W requires W_y=-P_y, which is retained before interval widening.

Original q is positive throughout, including y=1 where q^2=eta/2. It is stored as a formal original logq factor, never replaced by an exact zero. Actual rows use:

~~~
a_y=(6/5)*sigma_y
q_y=-(1+eta)*a_y/(2*a^2*q)
nu=2*(1+eta)/a
nu_y=-2*(1+eta)*a_y/a^2
~~~

The exact nu cancellation avoids differentiating a rounded q^2. On O2 slope only, b=t0=0 and a_Z=q_Z=nu_Z=a_yZ=q_yZ=nu_yZ=0. These facts must not be applied to axial or buffer routes. Source factors have the documented basis (logPstar, selected_logCstar, logL, original_logq, zero); radius powers are collected before large log arithmetic. Independent jet boxes enclose the original defined functions and are not selected as a new compatible field.

Four files with stem **lei_ren_part1_paper_compliant_current_original_O2_mixed_source_jets** are published under experiments/root_st073: producer, deterministic compressed report, checker and receipt. API: **OriginalO2MixedSources.source_frame(count,index,Z_lower=...,Z_upper=...)**, **describe(issued_frame)**, **OriginalO2MixedFrame**, **O2MixedAtlas**. Source levels must come from the accepted ordered mass archive; invalid levels, cells, axial domains and copied/unissued frames are rejected.

Focused evidence:4 independent defining-integral derivative identities;19 full original radial coefficient identities;5 independent ordinary-Z template identities;19 affine pressure/error and19 exact prefix/suffix derivative rows;4 exact q/nu identities and6 finite diagnostic derivatives;64 continuous source cells with2432 coefficient terms and1408 late pressure error terms;2304 mixed root rows;5 axis/signed frames and5 domain/frame rejections. **1232 Git-index dependency hashes PASS**. Producer 34.516s; checker 32.078s. Accepted ancestor suites were reused.

**Limit:** near y=1, the independent whole-cell q_y enclosure can be extremely wide because it divides by the microscopic positive q while separately enclosing sigma_y. This is a valid source enclosure, not a small derivative error claim. Mixed inverse/primitive evaluation, changed five integrals, all17/24 routes, actual controls, terminal closure, global N, joins/heat, stress/flat remainder, n-dependent recursion, two pulse families and corrected Cartesian NS remain open.

## Ordered next implementation

- [x] **O2-SLOPE-ACTUAL-PARAMETER-SOURCE-JETS:** complete original continuous C0/y/Z/yZ and actual a,q,nu with full pressure errors and positive endpoint q.
- [x] **O2-REGULAR-VARIABLE-q-FIXED-ANGLE (regular axis scope):** use the new issued frame and same atlas/ledger. Form u_i=(p2_i*q+p2*q_i)/dstar and all yZ terms, then r and hinv mixed jets. Reuse actual Fourier tails on the proved small-u predicate; do not divide by r or p2. Accept independent defining-integral comparisons including nonzero q_y.
- [x] **O2-FIXED-TRUE-PHASE-NU-TRANSPORT (regular axis scope):** use F=psi+T2-2*pi*phi*nu with phi fixed. Retain F_i=T2_i-2*pi*phi*nu_i, F_yZ and both t_i cross terms plus t_psi curvature. Keep source identities before caps; for slope 2*q*q_y=nu_y/2 is an exact useful cancellation. Do not use the old reference-only derivative adapter.
- [x] **O2-COMPLETE-A/B-PRODUCTS (regular axis scope):** recover A=a/2*(phi-psi/(2*pi)) and B=-a*E*totalT1/(4*pi) on slope, including all a/E first and mixed rows. Retain exact common symmetry traces without declaring a neighborhood exact-axis. Accept independent fixed-phi inverse derivatives and original T1/T2 quadratures.
- [ ] **O2-POSITIVE-SMALL-q-SIGNED-BRANCH:** derive a source-correlated weighted curvature bound valid for every positive q, including eta endpoint. The reference theorem requires q>=1/2. Prove actual p2 sign/carrier and regular/signed predicate coverage for O2; reference carrier evidence cannot be borrowed. No cap may become a selected field or a differentiated source.
- [ ] **O2-ACTUAL-C0/Z-FIVE-INTEGRALS:** connect actual common-N inverse phase, complete nonlinear densities, source y/Z transport, own-rate masses and inherited histories on whole source domains. Publish separate direct/averaged bounds, original large factors and global endpoints; local N candidates do not admit global frequency.
- [ ] **AXIAL/BUFFER-SOURCE-AND-SEAMS:** implement actual nonconstant a,b,t0,q,nu with their nonzero Z/yZ rows as required. Preserve original eleven-unit buffer, offsets, incoming moments and P0; typed seams must use the same source family.
- [ ] **RESTORE-ALL17/24-UNITS:** complete remaining original units using the predecessor route inventory; unsupported regions must remain explicit. Provide continuous C0/Z source/integral interfaces with actual error budgets and one common N.
- [ ] **SIGNED-NONLINEAR-MEANS/DRIFT:** integrate actual finite-N quadratic/exponential means with phase measure and y/Z drift. Keep signs and independent oracle errors; oscillatory reflection does not cancel nonlinear means.
- [ ] **FIVE-CONTROLS/TERMINAL/ALL-N:** use complete actual oracle constants in the centered all-N/unit-C1 bridge, derive global N inequalities with original large factors, and solve the five terminal identities as Z functions. Keep P0, pressure zero-rate memory and incoming corrections.
- [ ] **JOINS/EXACT-HEAT/ENERGY:** assemble the same-data core, annuli, flatten and collar with derivative joins, analytic preheat compatibility, exact heat exterior and finite physical energy/radial tails.
- [ ] **ADMISSIBLE-STRESS/FLAT-REMAINDER:** compute R_B=-div(T_B)+E_B and region-specific cone margins. Report stress and remainder separately; do not choose forcing after the residual.
- [ ] **ACTUAL-n-DEPENDENT-RECURSION:** implement distinct n=1 and n>=2 recovery equations, common core interval and independent moment repair, divergence-preserving cutoffs, truncation bounds and smooth sum. Rescaling one field is not coefficient recursion.
- [ ] **TWO-PULSES/CORRECTED-CARTESIAN-UVW:** realize both pulse families and averaged quadratic stress cancellation, then recover corrected u/v/w and validate forced Cartesian NS, physical contraction, relative axial elongation and material winding.

Mark DONE only with scoped code, evidence and commit. The full reconstruction objective remains active.
