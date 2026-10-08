# Complete original buffer period integration and new Rd/Rc histories

Checked source [a0e8c2a0](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a0e8c2a0c743a0856a22140a7962e7715394abb9); predecessor [CURRENT_COMPLETE_ORIGINAL_AXIAL_TRANSPORT_2026_10_08.md](CURRENT_COMPLETE_ORIGINAL_AXIAL_TRANSPORT_2026_10_08.md). Reused read-only worker metadata **GPT-5.6 Luna / max**. Root owns code, computation, acceptance and Git; the worker checked source/mean/error formulas and identified the next terminal target route. No new child was spawned.

## Implemented result

The entire original O2 buffer t in [0,11] now has actual nonlinear five-density C0/Z integral covers on the same Z[.36,.38] and [-.38,-.36] tiles at N1024. Each tile has 22 half-unit source cells, with 512 complete periods per cell and 16 full-phase bins. All **44 source cells**, **704 saved actual conditional phase/inverse density queries**, and **440 signed C0/Z local integral rows** execute and replay.

The new operator consumes the genuine new complete-axial exit correction from source 4033e3dc. It transports that incoming correction to Rd, then applies the already accepted complete O3 transition plus quiet power operator to obtain **updated Rc correction, own histories and absolute pressure in C0/Z**. The historical coarse tail Rd inlet is not reused. P0 remains separate and is added once; rate-zero pressure memory is retained. The expensive transition inverse/density calculation is not rerun because its source dependencies did not change.

This closes the previously missing buffer integration between complete axial transport and the checked Rd-to-Rc map on these two tiles. It does not establish five functional terminal identities, controls, axis/full Z, one global N, exact heat/admissible stress, true n-dependent recursion or corrected Cartesian NS. Those gates remain false. The integral estimates are enclosures of the original functions, not selected field values or final high-precision terminal solutions.

## Original buffer source and phase

With t measured forward from Y=exp(Md), the saved original endpoint background and original equations give

    E(t)=E_Y*exp(-t/2), V(t)=0,
    m(t)=mY*exp(-t),
    h(t)=hY*exp(-3t/2)+E_Y*exp(-3t/2)*expm1(t),
    k(t)=kY*exp(-3t/2),
    e(t)=(eY-E_Y^2*t/2)*exp(-t),
    p(t)=pY+E_Y^2*(1-exp(-t))/2.

The unchanged physical_mixed/raw_pre_stress_rows programs generate ordinary physical/primitive mixed4 and signed root y2/Z1 rows. The original radius is R=Rd*exp(t-11)=110*Cstar^10*Pstar^11*exp(t-11), retained in the same half-Pstar formal basis. All linear, cross, quadratic, energy and pressure terms remain, including full P0+p. The buffer has a=2, b=t0=Delta=0, q=sqrt(eta/2) and exact zero q slow jets; these buffer identities are not used on nonflat axial endpoints.

The saved Rd origin supplies phi(t)=Rd_phase+N*(t-11) modulo one. Every source cell contains exactly 512 cycles and each entire buffer **11264 cycles**. The geometry/phase origin have zero Z derivative. Actual origin intervals are retained; no phase midpoint is selected. Canonical phase bins cover every possible origin shift, including wrapped bin intervals.

## New nonlinear averaging with explicit error

For a fixed slow tuple, the original reflection A(1-phi)=-A(phi), B(1-phi)=-B(phi) gives the paired density at V=0, with x=A/N and z=B_over_Pstar/N,

    m_pair=0,
    h_pair=E*(cosh(x)-1),
    k_pair=E*z*sinh(x),
    p_pair=E^2*(cosh(2x)-1)/2,
    e_pair=z^2-p_pair.

Actual signed-u conditional inverse/first jets supply A/B before these nonlinear expressions are evaluated and branch covers are hulled. Eight half-period bin covers integrate this pair density into a frozen-slow mean enclosure. The nonzero factor in cosh(x)-1 stays formal through x^2*integral_0^1(1-s)*cosh(s*x)ds; sinh(x) keeps its x factor. The support theorem bounds only analytic factors and never defines A or B as a cap.

For period width Delta=1/N, let Kmax bound the original density at fixed phase, L_K its full fixed-phase slow-y derivative (including E_y=-E/2 and all source/implicit inverse effects), and M the positive original Duhamel mass for a source cell. Freeze the actual slow tuple at each period's left endpoint for the proof. Every such fixed-tuple mean lies in the computed whole-cell mean cover; no freeze point is selected numerically. The integrated C0 error satisfies

    |error| <= M*(Delta*L_K + 2*Kmax*(exp(lambda*Delta)-1)).

The first term controls the true slow envelope within each period. For the second, subtract the frozen unweighted mean and compare the varying period weight to its constant average; |K-mean|<=2Kmax. The weight term is exactly zero at the pressure rate lambda=0. Whole-cell derivative covers make this conservative, not a claim of a sharp left-endpoint expansion.

Z rows use direct nonlinear density derivatives in every full-phase bin. For any actual phase offset, each period spends exactly Delta/16 in a canonical bin, even if it wraps. The positive weight ratio obeys

    exp(-lambda*Delta) <= bin_mass/(period_mass/16)
                          <= exp(lambda*Delta).

Summing all complete periods gives the same uniform ratio for cell_mass/16. Original Z densities are integrated with these directed bin weights; no Z cancellation is claimed and no mixed yZ derivative is assumed available. Tightening these Z bounds with a justified higher-jet/reflection method is a precision task, not a reason to discard the complete original buffer source.

## Evidence and current limits

Six files under experiments/root_st073 use stem lei_ren_part1_paper_compliant_current_buffer_period_transport: producer/manifest, signed lossless gzip archives, checker/receipt. Archives total 22,530,932 compressed bytes. New focused checks include 60 original nonlinear pair comparisons (including a retained tiny nonzero cosh-minus-one factor), 72 independent analytic wrapped-bin mass checks, 6 weighted varying-envelope tests, 2430 original buffer physical/primitive comparisons and 216 signed stress y2/Z1/unit/radius comparisons.

The actual buffer checker recomputes original nonlinear C0/Z/y densities, reflected means, explicit error terms, positive per-rate masses, all signed local integrals and affine Rd/Rc propagation from lossless saved original inverse outputs. The expensive 704 phase inverse solves are not repeated; their accepted source program and all output/source bindings remain hashed. Complete ordered [0,11] coverage, exact 11-unit total coefficients, genuine new incoming correction and separate P0 are checked. **1296 Git-index dependency hashes PASS.** Producer 671.156s; focused checker 46.735s. No accepted upstream or transition producer/full legacy suite was rerun.

The half-unit slow hull and 16 phase-bin mean remain conservative. The Z integral deliberately retains phase and parameter correlation uncertainty rather than inventing cancellation. Actual terminal defect widths and control conditioning must determine whether any cell/phase refinements or mixed yZ calculations are required. Fixed-N two-tile coverage is not whole-Z/global construction admission.

## Authoritative next target interface

Use current_native_Rc_parameter_targets.py:target_rows and its original_joint_target_function_contract. The five original repair rows are M, D=(J-M)/mu, I, S, Cp; their distinct controls in current_generic_moment_repair_operator.py are axial0, axial2, swirl0, swirl1, swirl2. Let A=E_c, A_Z=E_c,Z and mu be the same original Rc amplitude/parameter, and let dm,dh,dk,de,dp denote the **correction** rows (not background-plus-correction own histories). The fixed-N normalization is

    M=dm/A, I=dh/A, S=de/A^2, Cp=dp/A^2,
    C=dk-A*dm, C_Z=dk_Z-A_Z*dm-A*dm_Z,
    D=C/(mu*A^2),
    D_Z=(C_Z-2*(A_Z/A)*C)/(mu*A^2).

The same contract has N^-1/N^-2 coefficient rows; the current exact finite-N integrals must be connected without assuming a coefficient split or replacing a function with a cap. Preserve joint C/C_Z at density/transport level before subtracting independently widened terminal hulls or dividing by tiny mu. Attach original_Rc_amplitude, the same CommonSourceCoordinates and source family. The control graph/Picard map is in current_native_Rc_functional_controls.py and consumes N_scaled_targets=N*r. The present result supplies continuous C0/Z function enclosures on two tiles, not a full-domain function-valued control oracle; that handle must be built or independently replayed over requested Z.

Raw preheat-pressure/exact-heat target constants are a separate bridge: current_pressure_terminal_closure.json retains the exact original P0 function, while angular/pressure terminal receipts still leave current_implicit_datum_to_native_raw_pressure_operator_identified and heat-pressure constant elimination false. Do not call those legacy constants a legal current-Rc target or redefine P0 to remove a residual.

## Executable next tasks

Only this leading handoff is active. Mark DONE only after committed implementation and necessary focused evidence; distinguish accepted from merged. Preserve source family, original P0 and normalized five units. Prior accepted producers should be consumed as saved data, not rerun for routine reassurance.

- [x] **BUFFER-ORIGINAL-BACKGROUND/STRESS/PHASE:** actual original endpoint histories, V=0 five ODEs, ordinary physical/primitive mixed4, signed stress roots y2/Z1, actual Rd phase and complete 11-unit geometry.
- [x] **BUFFER-ACTUAL-NONLINEAR-C0/Z-INTEGRALS:** all 22 cells per tile, real conditional inverse outputs, nonlinear reflected C0 means plus full slow/weight errors, original direct Z phase-bin integrals, all per-rate Duhamel masses.
- [x] **GENUINE-AXIAL-EXIT-TO-NEW-RD/Rc:** new complete-axial incoming correction, full buffer operator, accepted complete transition/power operator, updated Rc correction/own-history/absolute-pressure C0/Z, P0 once.
- [x] **ACTUAL-FIVE-RC-TERMINAL-COVER-ADAPTER (two tiles, N1024):** implement a source-linked RcTerminalDefect adapter from the new Rc correction C0/Z (not own histories), same A/A_Z/mu and original joint target contract above. Keep M,D,I,S,Cp units and all source errors. Attach the original N^-1/N^-2 coefficient contract only through a proved expansion/remainder interface; no relabelling of support caps or legacy target ranges.
- [ ] **CORRELATED-FIVE-DEFECTS:** form C=dk-A*dm and C_Z jointly from actual source/transport tuples before independent component hulls and mu division. Bind all five continuous-Z targets to the original control graph. Report signed covers, widths/conditioning and the dominant source contribution on both tiles. Add a function-valued C1 handle or independent source-function replay for requested Z; isolated samples/enclosures are not full-domain functional closure.
- [ ] **REFINE-WHAT-LIMITS-THE-DEFECT:** separate inherited slope/axial uncertainty, new buffer frozen mean/slow/weight/Z-bin errors, transition contributions and analytic datum/target uncertainty. Refine only the limiting cell/phase/error term. Preserve real source geometry and actual phase, not selected caps or midpoint values.
- [ ] **BUFFER-MIXED-yZ/SHARP-Z-AVERAGING (if material):** derive true mixed implicit inverse derivatives and full K_yZ bounds from original slow-root y2/Z1 data. Prove Z frozen-mean and weighted/slow errors before using cancellation. The current direct Z integral remains the accepted conservative baseline.
- [ ] **INDEPENDENT-ORIGINAL-FIVE-CONTROLS:** derive the five distinct control directions and their current continuous-Z derivatives against actual terminal defects. Solve them simultaneously, recording conditioning, uncertainty and effects on high derivatives/cone margin. Do not alter the shared P0 or add an arbitrary pressure tail.
- [ ] **CURRENT-DATUM-TO-RAW-PREHEAT/HEAT-BRIDGE:** identify the exact current implicit P0 operator with native raw-pressure normalization, then eliminate original heat/pressure constants using the same datum. Preserve P0/P0_Z and the original nonzero velocity source. Keep this false until the explicit operator identity, high-order joins and original target compatibility are proved.
- [ ] **AXIS/WHOLE-Z/GLOBAL-N:** construct the separate Z=0 source treatment, required edge/high-jet charts and one compatible finite frequency across the complete original domain. The current two strict-sign tiles at N1024 are only partial-domain evidence.
- [ ] **MATCHING/EXACT-HEAT/FINITE-ENERGY:** establish five functional terminal conditions, pressure compatibility, original high-order joins, axis regularity and physical finite-energy radial tails before admitting the exterior heat field.
- [ ] **ADMISSIBLE-STRESS/FLAT-REMAINDER:** recover -div(T_B)+E_B, all-chart cone margins and independent high-order/flat decay. Signed p1/p2 loop roots and period density estimates alone do not prove this gate.
- [ ] **TRUE-n-DEPENDENT-RECURSION:** implement n=1/n>=2 recovery equations, common core interval, independent moment repair per order, controlled truncation before curl and smooth sum. Ordinary source jets and coordinate rescaling are not recurrence.
- [ ] **TWO-FAMILY-CORRECTED-UVW:** implement both oscillatory families/mean corrections and quantified averaged stress cancellation, then export corrected Cartesian u/v/w, independent full NS/divergence diagnostics and measured contraction/elongation/swirl/material winding.

Complete original axial and buffer integration now feed a new Rc history map on the two actual tiles. The next main result must be actual five functional terminal defect/control construction, not another whole inherited validation pass.


Successor [CURRENT_RC_JOINT_TERMINAL_DEFECTS_2026_10_08.md](CURRENT_RC_JOINT_TERMINAL_DEFECTS_2026_10_08.md) implements the signed two-tile target cover adapter, all post-slope joint source transport and original linear control responses. The requested-Z function/control oracle and nonlinear terminal closure remain open. Its dominant-source ledger supersedes the former generic refinement instruction: prioritize slope source correlation/C1/pressure memory.
