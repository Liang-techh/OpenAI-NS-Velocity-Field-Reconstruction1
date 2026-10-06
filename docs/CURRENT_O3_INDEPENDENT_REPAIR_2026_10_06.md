# Current O3 independent source-correlated five-moment implicit repair

Current successor: [CURRENT_O3_REPAIRED_HISTORIES_2026_10_06.md](CURRENT_O3_REPAIRED_HISTORIES_2026_10_06.md), commit [90a587ac](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/90a587ac45de1a450af738c00b831543996c6d27). REPAIR5–8 and the local functional five-moment/velocity/pressure/radial q=2 exit (REPAIR9a) are now installed. Complete modified tensor, open continuation/physical dispatch/heat joins (REPAIR9b), common cone N and modified O2/O3 cones remain open. Earlier statuses below describe the implicit-map milestone.

Implementation and actual source receipts: commit [cb6b09cd](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/cb6b09cda6bbc05e3746969a457b3344ac8b2cd2).

## Constructed result

The actual modulated O2/O3 source now has a **new independent quiet-power five-bump map and a unique exact implicit constant repair vector**. The two axial and three swirl supports are fixed before sufficient-N selection. The complete map retains the same-center angular/axial cross terms, axial kinetic square and swirl energy/pressure squares.

The nearly coinciding axial rows and their actual source defects are both replaced by the same exact divided difference. Positive mu remains in every formula. The resulting scaled inverse and contraction proof are uniform despite the small mu. The actual sample N=10^12 satisfies the **repair-only** sufficient threshold 27,303,666, with contraction upper below 1.366e-5. This is not a sufficient common frequency for the modified cone or complete NS field.

The exact controls are constants in Z because the current defects share the C/C^2 factors and the common original repair amplitude cancels them. The source-dependent interval boxes enclose that unique implicit vector; they are not chosen coefficient values or midpoint substitutes. All five terminal increments vanish by the defining exact map. **This identity still needs connection to callable repaired partial-history/pressure/radial source functions before full repair installation is admitted.**

Original incoming moments and axis pressure datum are preserved. Original regional counts remain 15 strict nonzero whole regions plus separate exact zero exterior, with 17 original whole regions open. Inventory33/32/14 and primitive atlas14/8 are unchanged. Those counts do not admit the modified source. Actual coefficient recursion, waves, flat remainder, corrected Cartesian NS and finite-energy/dynamics goals remain open.

## Fixed independent supports and true radial units

Use the quiet O3 power band q in [1,2]. Let R0=Rw*exp(1), x=R/R0=exp(q-1), ell=1/40, and centers c0=1/5,c1=1/2,c2=4/5. These are **log-x centers**, not physical x centers. They lie strictly inside the quiet band. Distinct supports are disjoint; axial controls use the first and last swirl centers, so the two same-center cross terms remain.

The shared exact normalized even flat bump is

`beta(w)=raw_beta(w/ell)/(ell*B0)`, `B0=integral[-1,1] raw_beta(y)dy`,

where `raw_beta(y)=exp(-1/(1-y^2))` inside |y|<1 and zero outside. Define

`g_i(x)=beta(log(x)-c_i)/x`.

Then `integral x^p*g_i(x)dx=H(p)*exp(p*c_i)`, where `H(p)=integral beta(w)*exp(p*w)dw`. The normalization source and the change of variable bind **H(0)=1 exactly**. Fixed product weights are `integral x^p*g_i(x)^2dx=H2(p-1)*exp((p-1)*c_i)`. No old repair matrix or old-source defect certificate is reused.

The moment backend uses physical R, not absolute X. Any physical-R integral carries its R0 powers. The new implementation retains them through exact normalization instead of materializing enormous radii.

## Same actual amplitude and complete nonlinear map

At repair inlet t=log(R/Rd)=2, the original source has J(1)=1/2, U0=0 and

`A(Z)=E0(R0,Z)=Pstar*Ua*C(Z)*f2`, `f2=exp(-1-3mu/2)`, `C=1/(1+Z^2)`.

The source AST binds the actual transition endpoint and `pre.power` amplitude multiplication/exponent. Let alpha=1/2+mu, and set

`E=A*(x^(-alpha)+F)`, `Uz=A*G`.

The exact normalized repair changes in row order (M,J,I,S,Cp) are

| Row | Defining integral |
| --- | --- |
| M | integral G dx |
| J | integral sqrt(x)*(x^(-alpha)+F)*G dx |
| I | integral sqrt(x)*F dx |
| S | integral (G^2-x^(-alpha)*F-F^2/2) dx |
| Cp | integral (x^(-alpha)*F+F^2/2)/x dx |

Their physical primitive units are, respectively,

`(A*R0, sqrt2*A^2*R0^1.5, sqrt2*A*R0^1.5, A^2*R0, A^2)`.

The actual incoming modulated increments become `(Bm/f2,Bk/f2^2,Bh/f2,Be/f2^2,Bp/f2^2)`, independent of Z and Pstar. This concerns increments relative to the original histories; the original histories have not been set to zero.

At a later endpoint x_o, normalized history changes from these physical repair primitives are

`delta_m=A*M/x_o`, `delta_h=(A/Pstar)*I/x_o^1.5`,

`delta_k=(A^2/Pstar)*J/x_o^1.5`, `delta_e=(A^2/Pstar^2)*S/x_o`,

`delta_p=(A^2/Pstar^2)*Cp`.

These endpoint factors and Pstar sectors are required for the next installation.

## Correlated axial row and actual defect

The unscaled axial matrix has weights 1 and x^(-mu). Subtracting independently enclosed rows or defects and dividing by mu is unsound for useful bounds. Use the exact row

`D(c)=integral beta(w)*(exp(-mu*(c+w))-1)/mu dw`.

Its determinant is the stable defining expression

`det_A=-gap*H(-mu)*exp(-mu*c0)*integral_0^1 exp(-mu*gap*r)dr<0`, `gap=3/5`.

The three swirl powers are p=(1/2,-1/2-mu,-3/2-mu). Equal center spacing3/10 factors their block into common positive H rows and a Vandermonde matrix in exp(p_i*3/10). The S row sign is retained. The exact swirl determinant is strictly positive.

For the **actual modulated source defect**, at t=2 the normalized J-density divided by the M-density is

`exp(mu*(3/2-J(s))+A(s)/N)`.

With `A=mu*Abar`, `q(s)=3/2-J(s)+Abar(s)/N`, the scaled divided defect is

`dD=exp(-1+3mu/2)*integral exp(s/2-mu*J(s))*(beta(s)/sqrt(mu))*q(s)*integral_0^1 exp(mu*r*q(s))dr ds`.

This is formed **before** enclosing integrals. The proof replays the existing m/k scalar-density return programs, source transport factors exp(-2)/exp(-3), the f2 amplitude, shared cutoff/phase and actual A/Abar and beta/beta_scaled programs. It proves dD equals the scaled `(J-M)/mu` of the same existing source functions. A parallel plausible integral is insufficient.

## Scaled implicit solution and scope

Set `G=(sqrt(mu)/N)*(a0*g0+a2*g2)` and `F=(mu/N)*(e0*g0+e1*g1+e2*g2)`.

The scaled rows are `(M/(sqrtmu/N), (J-M)/(mu*sqrtmu/N), I/(mu/N), S/(mu/N), Cp/(mu/N))`. The scaled matrix B has the divided axial block and the three swirl rows. The nonlinear cross and kinetic coefficients are O(1/N); swirl-square coefficients are O(mu/N). No 1/mu blowup is introduced by enclosing the correlated defect.

The exact vector h=(a0,a2,e0,e1,e2) is defined by `B*h+Q_N(h,h)=-d_N`, within a fixed certified ball. For verified inverse norm CA, defect bound CD and quadratic norm CQ/N, set r=2*CA*CD. A sufficient repair-only integer bound is `N>4*CA*CQ*r`. The implementation proves strict ball inclusion and a strict contraction, and intersects directed fixed-point images to enclose the unique vector. It does not install box midpoints as fields.

API: `CurrentO3IndependentRepair(histories=checked_modified_histories).coefficients()`. Focused controller: `currento3independentrepair`. Reuse the checked warm graph; cold reconstruction is unnecessary for this stage. Files use the prefix `experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_independent_repair`: `_operator.py`, `.py`, `_check.py`, `.json`, `_check.json`.

Evidence: 39 exact source/normalization/physical/scaled-map identities;25 block-inverse identity enclosures; fixed support separation; complete overlapping cross/square terms; strict actual-N contraction; producer/checker, checked constructor/API, focused controller and compilation. Working/index source audit matched 834 dependency files. Read-only GPT-5.6 Luna/max review identified two source-binding gaps, both corrected before publication.

## Executable next tasks

- [x] **REPAIR1/MOD4:** fix independent log-translated compact supports strictly inside quiet power q[1,2], before common-N selection.
- [x] **REPAIR2/MOD5:** new divided axial determinant/inverse and source-dependent swirl Vandermonde block; preserve actual positive mu.
- [x] **REPAIR3/MOD6a:** full common-A normalized five-moment linear/quadratic map, including all overlap cross and both square terms.
- [x] **REPAIR4/MOD6b implicit controls:** actual signed source defect correlation, strict repair-only contraction and unique exact constant vector. Its defining map implies terminal moment cancellation; callable field connection remains next.
- [ ] **REPAIR5 — Bump source jets.** Implement g_i(x)=beta(logx-c_i)/x and ordinary logR derivatives through order4 with exact flat edges. Retain original A(Z)=Pstar*Ua*C*f2 axial jets through order5; constant controls have zero Z derivatives. Do not freeze A to an amplitude cap or midpoint.
- [ ] **REPAIR6 — Partial bump primitives.** Enclose continuous partial integrals for single weights p=0,-mu,1/2,-1/2-mu,-3/2-mu and product weights p=1/2,0,-1. Use the same exact H/H2 at full support and zero below support. Distinct support products vanish; retain both same-center axial/swirl products. Bind partial integral FTC to actual profile densities.
- [ ] **REPAIR7 — Own repaired histories.** Add original histories, transported modulation defects and partial repair primitives with all R0/current-R factors and Pstar sectors. Convert scaled control equations back to physical five-moment increments. No double counting or resetting incoming data.
- [ ] **REPAIR8 — Pressure and radial recovery.** Keep the original axis datum; integrate the repaired Cp in the same pressure. Recover radial velocity from own repaired M and Z derivative, with actual Uz restored from its Pstar sector. Then export callable ordinary logR/axial source rows.
- [ ] **REPAIR9 — Exit functional closure and flat joins.** At q=2 all bumps vanish with required derivatives and full primitives equal H/H2. Apply the exact implicit equation to prove equality of all five repaired histories with the original, for every Z[-1,1]. Connect that equality to pressure/radial/tensor source joins and inherited downstream analytic heat exterior. Only then promote the existing completed five-moment-repair gate.
- [ ] **REPAIR10/MOD7:** choose one finite N meeting derivative/shear, pressure, repair and full signed-cone inequalities. Repair-only threshold is insufficient. Preserve true frequency dependence in radial jets and any necessary Lipschitz constants.
- [ ] **REPAIR11/MOD9–MOD11:** whole affected O2/O3 source/tensor joins, independent O2 taper proof and modified closed O3 cone using complete energy, pressure, stress and remainder.
- [ ] **GLOBAL:** remaining core/O2/inner cones and smooth global admissible tensor/lift; actual n-dependent coefficient recovery with per-order moments, nonflat leading-origin cancellation and smooth sum; two actual oscillatory families with mean correction and quadratic cancellation; resolved physical u/v/w, independent corrected NS/energy and measured scale recursion/core widths/vorticity/material winding.

The full long-term objective remains active and incomplete. The immediate construction owner should implement REPAIR5–REPAIR9 on this exact source before attempting broad cone/NS validation.
