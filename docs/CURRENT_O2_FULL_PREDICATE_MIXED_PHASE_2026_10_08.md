# Full original O2 conditional varying-q mixed phase interface

Checked source [59e85ae7](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/59e85ae7d96f1b827a882ec8c431db870c989511). Predecessor: [CURRENT_O2_FULL_AXIAL_CARRIER_POSITIVE_Q_2026_10_08.md](CURRENT_O2_FULL_AXIAL_CARRIER_POSITIVE_Q_2026_10_08.md). The existing GPT-5.6 Luna/max worker reviewed read-only. Root implemented, computed and accepted; no child was spawned. The persistent full reconstruction goal remains active.

## Concrete implementation

The exact-axis restriction has now been extended by actual original-source predicates to all original O2 y[0,1], outer Z[-1,1], true phase[0,1]. Every one of64 radial source cells has three issued conditional branches: abs(u)<=1/4, u>=3/16, u<=-3/16. These predicates cover every real original u=p2*q/dstar, including zero; the overlap is3/16<=abs(u)<=1/4. There are192 source/phase frames, with original C0/y/Z/yZ A/B rows on each branch and18 real source inverse queries.

**This is a full conditional source/phase interface, not the completed O2 five integral contribution.** The separate angle/phase box is a conservative outer cover of the unique true inverse graph. No selected field or claim that an entire outer rectangle is regular/signed is made. All global-N, terminal, all-route, stress, recursion and corrected-NS admission flags stay false.

Eight published files under experiments/root_st073: producer, deterministic compressed report, focused checker and receipt for each stem:

- **lei_ren_part1_paper_compliant_current_original_O2_predicate_source_frames**: **OriginalO2PredicateSources.frame/describe**, **O2QUAtlas**, **O2PredicatePhase**.
- **lei_ren_part1_paper_compliant_current_original_O2_predicate_mixed_phase**: **OriginalO2PredicateMixed.primitive/whole_phase/evaluate**, **signed_fixed_phi_mixed**. The public adapter rejects copied or foreign frames and carries each frame's defining predicate.

## Source function and basis contract

Every root is rebuilt together onto one shared basis, context and ledger:

~~~
(logPstar, selected_logCstar, logL, original_logq, defined_original_logabsu_or_zero).
~~~

Actual radius powers were already collected into P/C before this basis. The legacy scale record's radius_power field means logabsu power only on this documented atlas; it is not compatible with a prior reserved-zero/radius slot. Addition uses maximum P/C/u and minimum L/q powers. Relative u exponents are nonpositive and logabsu has lower>=-2, so their positive exponentials remain bounded even when u powers are negative. Relative q exponents are nonnegative and original logq<=0. Formal q and u factors are not materialized or floored.

The accepted actual carrier is p2=Z*g, g/Lambda0<0 with Lambda0=Pstar^11*Cstar^10. The signed predicate implies a necessary lower q band from abs(p2)<=Lambda0*max(abs(g/Lambda0)), abs(Z)<=1. That restricts the same original q function on the named source domain; it does not prove the outer rectangle satisfies the signed predicate. The regular branch keeps the original q band. No double normalization of g occurs.

On every named predicate, p2=dstar*u/q, p2_y=p2*(g_y/g), and:

~~~
u_y = u*(g_y/g+q_y/q)
u_Z = p2_Z*q/dstar
u_yZ = (p2_yZ*q+p2_Z*q_y)/dstar.
~~~

Only strictly positive q/dstar and nonzero g are divided, so the regular identities extend through p2=u=0. Original p2_Z/yZ templates and their pressure errors remain. The saved g_y/g bound is an enclosure of the derivative of the same original carrier under its source identity; it is never selected or differentiated as an exact C1 function. Original pressure family, P0 datum, incoming histories, positive endpoint q and original varying a/q/nu y rows remain.

## Actual fixed true phase mixed formulas

The regular branch reuses the accepted regular derivative code unchanged. Its finite-offset arithmetic callback now allows large negative flat-source logs, which had rejected the last two source cells. Positive large offsets are still rejected until collected. Directed small exponential tails keep every formal q power; this is not a q floor or exact-flat replacement.

The signed branch constructs original r,s,hinv and fixed-angle chi jets, then T1=q*hinv/r*(chi-psi) and T2=q^2/r^2*((2-3*s)*chi+s*psi+2*r*sinchi) using the actual q MixedJet. K_y=r*(g_y/g+q_y/q) cancels the common u factor before bounding. Fixed true phi uses:

~~~
F = psi+T2-2*pi*phi*nu, nu=1+2*q^2
F_y = T2_y-2*pi*phi*nu_y; F_Z=T2_Z; F_yZ=T2_yZ
t_y = (K_y*(2*coschi-r)+q_y/q)*t-2*q*K_y*hinv.
~~~

The complete implicit mixed chain and total T1 mixed chain retain both cross terms and curvature. On the true inverse graph:

~~~
F_y=K_y*J+(q_y/q)*(4*pi*phi-2*psi), F_Z=K_Z*J
abs(4*pi*phi-2*psi)<=4*pi.
~~~

Both the main J^2 and the extra varying-q J curvature terms use the accepted all-positive-q theorem; the primitive J^2 bound retains the original formal q factor. The identity follows from nu_y=4*q*q_y and2*T2-8*pi*phi*q^2=4*pi*phi-2*psi on F=0. It is not asserted on arbitrary independent angle/phase pairs.

All a and E first/mixed product terms are retained in A=a/2*(phi-psi/(2*pi)) and B=-a*E*T1_total/(4*pi). Common0, half-period and full-period traces are exactly zero on all three branches. nu is phase normalization, not physical NS viscosity.

## Focused acceptance evidence

Source interface:4 source restriction identities;72 separate q/u arithmetic endpoint comparisons, including positive and negative powers;192 issued frames and6912 actual mixed root rows; shared basis/rebase/issued guards and original varying q/nu y rows. Producer 41.765s, checker 31.969s.

Mixed phase: exact same original rational direction/T2 derivative and period identities, fixed-phi implicit/total-T1 chain and varying-q cancellation; **228** independent rational defining-integral/implicit/A/B comparisons, **32** common-function overlap A/B comparisons and **12** independent finite inverse brackets. Both u signs, original flat a/q/nu variation at finite diagnostic eta and positive q<1e-6 are included. Actual original continuous coverage has192 frames, 1536 A/B rows,18 source inverse records and9 exact symmetry traces. Producer 78.953s, checker 149.937s. Finite parameters are diagnostics only, not replacements of the original data.

**1252 unique staged dependency hashes PASS.** Accepted ancestor receipts are reused; no ancestor suite was rerun. A read-only review of signed mixed chains and the source basis found no remaining concrete blocker after the nu_y factor accounting was resolved.

## Ordered tasks and acceptance

- [x] **O2-ISSUED-REGULAR/SIGNED-SOURCE-FRAMES:** original conditional abs(u)<=1/4 / u>=3/16 / u<=-3/16, both physical sides and zero, overlap, original mixed roots and histories; source59e85ae7.
- [x] **O2-q/LOGABS-u-COLLECTED-ATLAS:** separate original q/u logs, exact radius collection, shared rebase and relative-scale addition, explicitly documented fifth-slot semantics.
- [x] **O2-CORRELATED-Y-ROWS:** source p2/u identities, g_y/g+q_y/q, original transverse jets and K_y cancellation; no derivative of hulls or selected ratio.
- [x] **O2-SIGNED-VARIABLE-q-FIXED-PHI/FULL-A-B:** actual T1/T2, q_y/nu_y, both curvature pieces and all a/E products, exact common traces and source inverse queries.
- [x] **O2-REGULAR/SIGNED-COMMON-INVERSE-AND-DOMAIN-UNION (interface):** same rational integrands/primitive normalization and unique inverse, continuous named predicate cover. This checkbox does not admit changed five integrals.
- [ ] **NEXT O2-ACTUAL-N-C0/Z-DENSITIES:** implement stable F=E*A*exprel(A/N) and F_Z=E_Z*A*exprel(A/N)+E*A_Z*exp(A/N). Include all five complete nonlinear density products and cross terms from original Section11. Use one explicit candidate N and original source inputs. Accept only source-backed C0/Z rows and independent defining-density comparisons; no graph-bound derivative or constant-q substitution.
- [ ] **O2-PREDICATE-DENSITY-UNION:** evaluate all three source branches as restrictions of the same density function. Normalize/export each local q/u atlas before cross-frame sums; hull competing branch bounds where predicates overlap, never add duplicates. Preserve positivity and named domain labels. Acceptance: axis, both signs and all64 source cells, union cover receipts and foreign-frame rejection.
- [ ] **O2-SOURCE-COMMON-N-ENDPOINT-PHASE:** bind the actual original radius/common-N phase at each source endpoint. Retain literal source constants, original phase convention and global endpoint terms. Fixed-phi source derivatives support averaging; they are not complete derivatives of a spatial rapid phase. Acceptance: endpoint source identity and actual inverse, no independently chosen endpoint phases.
- [ ] **O2-FIVE-OWN-RATE-INTEGRALS-C0/Z:** combine complete densities, source y/yZ transport, directed own-rate masses, source/oracle errors and phase averaging/direct remainders across the full original window. Retain m/e rate1, h/k rate3/2 and pressure rate0 memory. Report signed contributions and bounds separately. Acceptance: actual full-window C0/Z contribution at an explicit candidate N, both signs and axis, not a finite axial sample terminal theorem.
- [ ] **O2-INCOMING-HISTORIES/P0-TRANSPORT:** consume original inherited five defects and common analytic pressure datum. Keep original large carrier/radius factors and all global boundary terms. A zero incoming test is not the actual driver. Acceptance: unchanged family/history hashes and exact pressure zero-rate memory.
- [ ] **AXIAL/BUFFER-MIXED-SOURCES:** recover remaining original charts with actual nonconstant a,b,t0,q,nu and C0/y/Z/yZ rows. Retain eleven-unit buffer offsets and joins. Do not reuse the O2 b=t0=0, q_Z=0 guards outside O2; explicit unsupported states are required.
- [ ] **ALL17/24-CONTINUOUS-ORACLE:** bind every original chart and cell to one pressure/history family, source-backed inverse/primitive/density/integral interface and source/oracle error. Continuous domain cover is required; no finite sample or reference profile substitutes the full route.
- [ ] **FINITE-N-MEANS/DRIFT:** integrate signed nonlinear means with actual phase measure and retain their finite-N biases. Attach all incoming/outgoing endpoint and quadrature errors to the centered all-N bridge. Acceptance: complete five function defects with C0/Z bounds, not local residual reduction.
- [ ] **FIVE-TERMINAL-CONTROLS/GLOBAL-N:** solve the five terminal identities as Z functions, with original large-scale inequalities admitting a common global N. Separate source, quadrature, oracle and Picard errors. Local N=7 or another candidate does not satisfy this task.
- [ ] **SAME-DATA-JOINS/ANALYTIC-PREHEAT/HEAT/ENERGY:** assemble inner, connecting annulus, flatten and collar with derivative joins, analytic preheat compatibility, exact heat exterior and finite physical energy/tails.
- [ ] **ADMISSIBLE-STRESS/FLAT-REMAINDER:** compute original region-specific cone margins and R_B=-div(T_B)+E_B with separate high-order/flat decay diagnostics and norms. Do not impose corrected residual tolerance on the unfinished background.
- [ ] **ACTUAL-n-DEPENDENT-RECURSION:** implement n=1 and n>=2 recovery equations, common core interval, separate moment repairs, divergence-preserving cutoffs, finite-order remainder and smooth sum. Scaling a single field is not recursion.
- [ ] **TWO-OSCILLATORY-FAMILIES/STRESS-CANCELLATION:** recover both pulse families and their averaged quadratic stress cancellation with the original admissible stress. Acceptance includes finite-frequency correction errors.
- [ ] **CORRECTED-CARTESIAN-UVW/PHYSICAL-DIAGNOSTICS:** assemble full corrected u(x,y,z,t),v,w and independent forced NS residual. Quantify radial contraction, relative axial elongation, amplitude/vorticity scaling and material winding; publish finite-energy and residual evidence.

The full goal remains active. Every DONE item must have scoped code, evidence and a commit. The next action is the actual O2 density/integral layer, not another prerequisite-only certificate.
