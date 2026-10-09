# Current original R110-to-Rsh finite-N source and local Duhamel

Checked source [2f240f3c](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/2f240f3ceb83d5a9ef155358db2f8134cd9b43ff). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [Rm weighted Z averaging and prefix starting points](CURRENT_ORIGINAL_RM_WEIGHTED_Z_AVERAGING_2026_10_09.md).

**RM-W4b has its first genuine conditional local implementation.** The full original long-reshape window now supplies common-unit source histories, original finite-N signed five-density C0/Z covers and a local inhomogeneous Duhamel operator at candidate N257, conditional on the current accepted Z=0,.5 source frames. This is production of a missing source window. It is not a complete finite-N prefix, an actual Rm inlet, five-moment closure or coefficient recursion.

The actual finite-N R110 incoming correction is still **unsupplied**. Its nonzero memory multiplier is retained explicitly, instead of substituting background histories or an earlier N1024 correction. Reference, restoration and postrestore candidate-background source drivers are still required.

## Current source and full finite kernels

The owner uses the current accepted fixed-radius R110 source, its actual anchored B axial tuple, frozen T=400A, logC, and the same live flow, five formal log bases and ledger as the current Rm parameter wrapper. Source family, frame, exact inlet radius110 and the B/T/logC tuples are guarded. Canonical analytic P0 coefficients agree exactly; recovery uses the same live P0 object as the reference owner.

For y=log(R/110), rho=y/T and R=110 exp(y), the actual common swirl amplitude is

    E=exp(B*(1-sigma(rho))-log(1+Z^2)+y/10-logC-logPstar)
    E_y=E*(.1-B*sigma'(rho)/T)
    a=.8+2*B*sigma'(rho)/T
    C=E*a.

Absolute large positive factors and exponentially small tails stay formal; exp(T) and the physical radius are not materialized. The correlated C identity is preserved before generic recovery. V is the actual current R110 physical axial source and V_y=V_yZ=0 in this window.

Each shape retains its true inherited leading-history term and a full finite source kernel

    K(k,m;y)=integral_0^y exp(-k*t+m*B*(sigma(y/T)-sigma((y-t)/T)))dt.

The pairs (k,m) are theta(1.6,1), pressure(.2,2), swirl(1.2,2). Positive rate bounds follow from the actual B/T theorem; ordinary axial Bell terms through order5 enclose the same integral. The finite upper limit is not replaced by infinity, shortened to a small interval or converted to a selected terminal kernel value. Both formal finite tails are retained. Mean/axial leading source terms and their inherited decay exp(-y) also remain.

The kernel C0 mass uses a cancellation-safe width*exp_average(-rate*width) cover for small arguments. Large original widths retain the formal tail expression. This allows rational microscopic subdivision without turning a positive source mass into zero.

## Original recovery and density units

Raw histories are built from the actual six shapes:

    m=mean
    h=E*theta
    k=E*theta_z
    e=axial/Pstar^2-E^2*swirl/2
    p=E^2*pressure/2.

Raw m,k,V are in physical axial units. The unchanged original generic recovery divides those three by Pstar exactly once. p is only the radial pressure increment; P0 is added separately in full signed inertia. Exact original first-y equations are

    m_y=V-m
    h_y=E-3*h/2
    k_y=E*V-3*k/2
    e_y=V^2/Pstar^2-e-E^2/2
    p_y=E^2/2.

The AST-bound recovery changes only the explicit long-window radius and correlated C assignments. All original signed meridional, inertial and pressure assignments are preserved. Pre-radius p2 uses the full signed pressure and inertia; one positive R factor is attached afterward. This window has b=t0=0 and Delta=a-2<0. q is computed from the same correlated source, with the current parameter eta and d_star. Parameter definitions do not establish the new source's stress-cone admission.

The original all-signed-u inverse theorem supplies A/B C0/Z/phi covers. The exact finite-N density expression keeps deltaE=E*expm1(A/N) and deltaV=B/N, including products and quadratic bias in all five densities. Actual radius phase is frac(N*(log110+y-logRa-hb*s_c/2)), with exact phase_Z=0. A full-period cover encloses it here; no narrow inverse or actual-phase cancellation is claimed.

## Local affine correction with nonzero source drivers

The full rho partition is 0,1/4,1/2,3/4,1. Physical-log-radius width and suffix are T*delta_rho and T*(1-rho_right). Each cell integrates source majorants with the original five own rates

    (m,h,k,e,p)=(1,3/2,3/2,1,0).

For rate lambda>0, mass=(1-exp(-lambda*width))/lambda, with a cancellation-safe bounded-argument implementation. For pressure, mass=width and all memory decays equal1. Each driver bound keeps its downstream suffix exp(-lambda*suffix), then the five signed covers are summed. The result is an enclosure of the true local inhomogeneous integral, not a selected numerical defect value.

    deltaH_j(Rsh)=exp(-lambda_j*T)*deltaH_j(R110)+local_driver_j.

The same relation holds for the first axial derivative because T, radius and phase are Z-independent. All ten incoming-memory rows are nonzero formal expressions. None, zero vectors and leading background histories are not an admitted finite-N incoming correction. A source-bound affine application and a real current N257 R110 vector remain tasks.

No small terminal defect, numerical convergence, functional closure or admitted global frequency follows from these conservative covers. No local weighted-IBP improvement over the long window has yet been measured; this layer installs its genuine direct baseline.

## Scoped evidence

Artifacts: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_reshape_finite_N.py, .json.gz, _check.py, _check.json. Final producer 12.719s and checker 21.094s, both terminal exit0. **1071 staged exact dependency hashes PASS**. Only this new producer/checker and staged byte audit ran; ancestor producers/checkers and old physical owner constructors were not executed. The existing read-only worker is **GPT-5.6 Luna / max**; root owns code, compute, commits and final scope.

Independent source checks include 36 full original kernel coefficient quadratures through Z order5, 37 common-unit amplitude/history/inertial/root comparisons, 20 original finite-N density C0/Z comparisons, 20 own-rate mass/decay/signed-forcing comparisons and 8 positive microscopic-width comparisons at width1e-1000. The scalar inverse Z difference is a diagnostic, not an exact derivative proof. Its auxiliary p1 is unused by A/B and does not assert actual source inertia or cone admission. Exact symbolic calculus plus live source enclosures establish the conditional bound.

Live/report replay covers 8 full-window source cells, 336 actual source coefficients, 80 density rows and 10 retained memory rows. Six wrong family/frame/radius source bindings are rejected. All full reconstruction gates remain false.

## Concrete next source map

[Current original reference/restoration functions](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reference_restore_functions.py) supply the current leading six shapes and shared P0. Their packet's actual_E_V110_minus4Z is the centered axial velocity forcing, **not** the generic shear amplitude. Generic E=Utheta/Pstar comes from log_Utheta_over_Pstar_axial5 and has E_y=E/10 on these windows.

Reference: g=gap*rho, gap=10(logC+logPstar)-T-8>0, logR=log110+T+g. Integrate with physical dy=gap*d_rho. Restoration: R=Rz*exp(t), t in[0,1], dy=dt. Postrestore: R=Rref*exp(x), x in[-7,-6], dy=dx, ending at Rm. Packet y rows already refer to physical log radius; do not divide them by gap.

Reference and postrestore have V_y=0. Restoration has V_y=(V_Rsh-4Z)*alpha_1, alpha_1=-sigma'(t); obtain the forcing and alpha from actual_E_V110_minus4Z and actual_alpha_ordinary_y_derivatives. Restore general b,t0,kappa and Delta from the full signed source. The long-window b=t0=0 / negative-Delta adapter cannot simply be reused throughout restoration.

Current packet methods are point evaluators. Endpoint values and selected caps cannot be promoted to whole-cell source functions. Build analytic interval cell providers with exact inherited memory, original full cutoff kernels, full axial5 rows and original raw first-y ODEs; then use the signed density service with source-bounded q branches. Accepted current Original* wrappers that hydrate receipts and pure function owners are permitted. Avoid ancestor recomputation and old physical assemblies/CompliantLongReshapeProfiles.

The six centered background rates (1,1.6,1.6,1,1.2,.2) are not the five finite-N correction rates. Keep the two sets separate. Every nonzero nonlinear candidate-background source must be integrated or proved identically zero; homogeneous tails alone do not complete any prefix window.

## Detailed production queue

- [x] **RM-W4b.1 CURRENT LONG SOURCE:** same R110 family/frame/radius/B/T/logC/P0, original six full-finite shapes, raw/common units and exact five first-y ODEs.
- [x] **RM-W4b.2 LONG SIGNED DENSITY:** source-owned a/q/E/p2 C0/Z, full signed pressure and inertia, one radius factor, original finite-N density with both N levels and exact phase-Z independence.
- [x] **RM-W4b.3 FIRST LOCAL LONG DUHAMEL:** all four full-T cells and five own-rate suffix weights, nonzero nonlinear drivers, microscopic stable masses and explicit incoming memory. Conditional Z=0,.5 / N257 only.
- [ ] **RM-W4a.1 ACTUAL R110 FUNCTION/FREQUENCY ADAPTER:** compare the native bridge/switch correction definitions with current endpoint functions, parameters, datum and ledger. A matching family string or arithmetic rebasing is insufficient. Saved N1024 is not current N257.
- [ ] **RM-W4a.2 PRODUCE REAL CURRENT N257 INLET:** recover actual candidate-background source corrections through every core/bridge/switch/power prefix window, integrate all original density sources and export five signed C0/Z rows at R110. Preserve each true transition endpoint and radius measure. Reject missing windows, wrong P0, reduced-basis transplants and cross-N data.
- [ ] **RM-W4b.4 TYPED AFFINE APPLICATION:** accept only a current-source/P0/frequency-bound five-row C0/Z correction, apply every own-rate memory and local driver, and export an Rsh correction with endpoint identity. Reject None, arbitrary zero placeholders and leading histories. No actual inlet admission before real input exists.
- [ ] **RM-W4c.1 REFERENCE CELL PROVIDER:** implement genuine interval cells over the full positive gap. Keep current Rsh inherited histories, centered source forcing and formal tails; compile E,V, five raw histories and current P0. Radius phase and dy=gap*d_rho must use the actual current source.
- [ ] **RM-W4c.2 REFERENCE SIGNED DRIVER:** recover same-owner generic source / q / original A/B, integrate all five finite-N C0/Z density drivers with physical-log-radius widths and suffixes. Export an affine Rz result retaining unknown incoming correction.
- [ ] **RM-W4c.3 RESTORATION CELL PROVIDER:** full original cutoff integrals on t in[0,1], actual nonzero V_y and general b/t0/kappa/Delta/q branches. Consume sufficient true axial orders; never differentiate selected bounds. Keep current P0 separate from the radial pressure increment.
- [ ] **RM-W4c.4 RESTORATION / POSTRESTORE DRIVERS:** integrate candidate-background sources with finite-N rates (1,1.5,1.5,1,0), carry all boundary memories across Rz and restoration exit, and compute the postrestore contribution to Rm. Prove each shared source endpoint identity.
- [ ] **RM-W4c.5 FULL PREFIX COMPOSITION:** compose true current R110 input plus long, reference, restoration and postrestore drivers. Export five incoming C0/Z correction functions at the current Rm owner. Admit the inlet only when every window, datum, frequency and endpoint matches; otherwise report the exact missing window.
- [ ] **RM-W5 SHARPER LONG/RM COVERS:** conditioned original inverse derivatives, valid correlation/subdivision and endpoint-retaining averaging. Quantify width gains relative to this genuine baseline; large upper-bound reductions are not convergence. Shared-source identity is required before endpoint cancellation.
- [ ] **RM-W6 COMPLETE Rc/ALL24:** actual outer-window contributions, complete Rc_E/Rc_E_Z and terminal defect functions using the genuine current prefix. Both original N levels and all suffixes remain.
- [ ] **RM-W7 WHOLE-Z / SPATIAL / HIGHER:** extend conditional source frames to functions, establish actual branch/cone/cutoff/pole/midplane guards, and add A_phiZ/B_phiZ, yy/ZZ and enough genuine source orders with fast N chains.
- [ ] **RM-W8 FUNCTIONAL FIVE-MOMENT / CONE / ONE N:** solve complete B*h+N*r+Q/N from real cumulative inputs, prove uniqueness and all five terminal identities throughout Z, actual stress margins and one consistent admitted frequency.
- [ ] **RM-W9 HEAT / PRESSURE / ENERGY:** compatible functional analytic preheat pressure, exact heat joins and controlled finite-energy tails.
- [ ] **RM-W10 REAL COEFFICIENT RECURSION:** admissible stress / flat remainder, distinct n=1 and n>=2 recovery, independent moment repair, curl-based cutoff and controlled smooth summation. Coordinate scaling alone is not recursion.
- [ ] **RM-W11 CORRECTED NS / DYNAMICS:** original mean/pulse families, averaged stress cancellation, smooth forcing and independent Cartesian residuals; measured contraction, relative elongation, accumulated material winding and finite energy.

Mark DONE only with code, scoped report/receipt and a commit. Prioritize genuine source production; no ancestor reruns are needed for this layer. Full goal remains active.
