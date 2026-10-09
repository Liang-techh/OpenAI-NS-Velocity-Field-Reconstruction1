# Successor: factored actual R100 endpoint and original pressure interface

Source [04dc4c9a](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/04dc4c9a4570c620f23562eb049b24583469357b) adds the exact R100 actual F/V, six histories, Q0..4, anchored amplitude derivative ratios and separate physical pressure interface. The accepted actual macro moment integrals below remain unchanged. Follow [CURRENT_ORIGINAL_R100_ENDPOINT_2026_10_08.md](CURRENT_ORIGINAL_R100_ENDPOINT_2026_10_08.md) for the actual first/second switches, true R110 functions and downstream tasks. Whole-Z, global N and actual coefficient recursion remain open.

---

# Original frozen-macro actual six moment functions

Checked source [d3eb501b](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d3eb501b7573fbe00009ccb79fb67becdf921fdb). Full reconstruction **ACTIVE / INCOMPLETE**. The actual frozen macro now has callable H/M/K/A/B/C finite Volterra integrals and their defining radial ODE rows, using the original nonlinear phi/V sources and the actual micro inlet errors. Native frames Z=0,.5 and exact rational macro fractions are admitted. Whole-Z source functions, both micro function providers, switches and the complete bridge remain open.

The accepted [original macro F/V functions](CURRENT_ORIGINAL_BRIDGE_MACRO_FUNCTIONS_2026_10_08.md), [finite terminal reshape kernels](CURRENT_ORIGINAL_RESHAPE_TERMINAL_KERNELS_2026_10_08.md), [general fixed-phase Z derivatives](CURRENT_GENERAL_CONDITIONED_SLOW_Z_2026_10_08.md) and terminal/reference/O2 partial histories remain reusable.

## Defining original histories

Use t=log(R/R0), R0=Ra*exp(2hb), S=q*(log(100/Ra)-2hb). Original radial geometry, h, Ra and hb are fixed with respect to Z; dR/R=dt. At q=1 the exact original endpoint is R=100. For each actual moment X,

    X(S)=exp(-rate*S)*X_micro_exit
         + integral_0^S exp(-rate*(S-t))*source_actual(t)dt.

| Moment | Rate | Actual source | Constant-source equilibrium |
| --- | ---: | --- | --- |
| H | 2 | 2 phi | phi |
| M | 1 | V | V |
| K | 2 | 2 phi V | phi V |
| A | 1 | V^2 | V^2 |
| B | 2 | phi^2 | phi^2/2 |
| C | 1 | phi^2 | phi^2 |

Inlets are the accepted second_exit.actual_own_six_moments_axial5 materialized ordinary Taylor coefficients of normalized moments, not comparison moments or hidden amplitude-factor packets. Nonzero incoming decay, core/micro source errors and the same original source family are retained. At q=0 the actual histories return the exact admitted inlet rows. No incoming moment is reset. Analytic P0 remains separate.

## Signed finite exponential polynomials and full error

The evaluator collects signed terms sum_(a,k) c_(a,k)(Z)*t^k*exp(a*t) before enclosure. Each c is an ordinary Z0..5 coefficient row in the same factored source basis and ledger. Products phi*V, V^2 and phi^2 are formed with this algebra before complete finite integration; generated integer exponents span[-6,8] and degrees0..6.

The zero-at-zero primitive of t^k*exp(a*t) is exact, including a=0. The complete Volterra mass uses b=a+rate. If b=0 it is exp(-rate*S)*S^(k+1)/(k+1). Otherwise it is

    exp(a*S)*sum_(j=0)^k (-1)^j*k!/(k-j)!*S^(k-j)/b^(j+1)
    - exp(-rate*S)*(-1)^k*k!/b^(k+1).

Formal radius factors are collected before log evaluation: exp(a*S)=R1^a/R0^a. No enormous exp(S) is evaluated, and no incoming decay is erased. Every exact integer resonance has an explicit branch.

For the original ell, the signed approximation is Epoly=1+ell+ell^2/2. It encloses the complete exponential with a nonzero remainder etaE<=B^3*exp(B)/6, where B is the accepted complete-prefix weighted ordinary Taylor norm and B<=1/2. This is a complete finite-integral enclosure, not a claim that exp(ell) equals Epoly. Actual phi and V are built with the original core-exit jet, micro ell/deltaV, known finite-width comparison modes, q0 and h/Pstar^2/F0^2 scales applied once. Native swirl derivative rows already carry normalized F0^2 derivatives; the base amplitude is not dressed twice.

Let etaPhi/etaV bound the complete field errors and PhiNorm/VNorm bound their signed polynomial parts. Source error norms are

    H:2etaPhi; M:etaV
    K:2(PhiNorm*etaV+VNorm*etaPhi+etaPhi*etaV)
    A:2VNorm*etaV+etaV^2
    B,C:2PhiNorm*etaPhi+etaPhi^2.

Multiply each by the positive complete Volterra mass (1-exp(-rate*S))/rate exactly once. Ordinary coefficient-n errors divide by the common positive Taylor weight to power n. The actual radial derivative rows are source_actual(S)-rate*X(S), using the defining ODE and complete field errors. They are not derivatives of a midpoint, interval selector or unrelated comparison trajectory.

## Source binding and evidence

OriginalBridgeMacroMoments uses the accepted OriginalBridgeMacroFunctions factory, source/cache/receipt hashes, same h/radii/basis/ledger and actual micro inlet packet. No original ancestor constructor, producer or full upstream checker is rerun. Generic ActualMacroMoments callers must supply valid same-owner inlet and defining source rows; algebra alone does not certify arbitrary caps as physical sources.

Producer/report/checker/receipt: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_bridge_macro_moments.py, .json, _check.py, _check.json. Producer 56.688s; checker 71.016s; terminal exit0. **308 staged exact dependency hashes PASS**.

Independent checks pass: 40 resonant/nonresonant Volterra quadratures, 12 primitive quadratures and 210 symbolic ODE/zero-inlet identities for the complete exponent/degree/rate set. Independent diagnostic finite fixtures use the full nonlinear exponential with both signed directions, nonzero swirl, pressure and incoming data: 144 actual history Taylor integrals, 144 radial ODE comparisons and 48 field comparisons. These fixtures are not native parameters.

Native evidence covers 216 history Taylor records, 216 radial ODE records, 24 nonzero complete nonlinear integral error norms, 36 retained incoming decays and 12 unchanged actual micro inlet joins. Ten invalid source/history/geometry/rate requests reject. Read-only mathematical review: **GPT-5.6 Luna / max**, scoped PASS; root owns code, compute and Git.

Full bridge/switch/R110/active patch/all-integral/five-control/global-N/stress/real coefficient recursion/corrected-UVW flags remain false. The result closes actual frozen-macro history functions conditional on the admitted actual micro inlet enclosures. It does not close the entire physical reconstruction or demonstrate scale recursion.

## Detailed production queue

- [x] **ORIGINAL MACRO F/V SOURCE FUNCTIONS.** Accepted finite-width known comparison modes, anchored original amplitude, ordinary Z0..5, signed nonlinear source transport and complete nonzero errors.
- [x] **ACTUAL MACRO SIX MOMENT FUNCTIONS.** All six complete finite Volterra integrals, signed quadratic products, nonzero inlet transport, ordinary Z0..5, defining radial ODEs and independent complete nonlinear checks.
- [ ] **NEXT: EXACT R100 ENDPOINT ADAPTER.** Export the admitted actual macro endpoint phi/V, H/M/K/A/B/C, Q, same original F0 amplitude ratios and separate P0 to the original first-switch source interface. Recover every normalization from defining source equations; retain factored radii, microscopic displacement, source family and derivative/error rows. Acceptance: no cap/selector substitution, same-source R100 inlet and no false R110 closure claim.
- [ ] **FIRST ORIGINAL SWITCH FUNCTIONS.** Recover 0<=t<=1 at R=100*exp(hb*t), using the original hb^2*(1-sigma(t)) drive, angular ratio, tiny log-phase displacement and original Q/P0. Integrate actual phi/V and all six own moments, preserving dy=hb*dt and fixed-Z derivatives. Acceptance: independently checked actual source functions and inlet/exit joins with nonzero errors.
- [ ] **SECOND SWITCH AND TRUE R110 HISTORIES.** Continue the second original switch; bind the exact terminal geometry and the true inherited theta_z/V/pressure histories. Acceptance: actual R110 phi/V and H/M/K/A/B/C functions, not symmetric110/12100 caps or comparison histories.
- [ ] **TRUE BOTH MICRO FUNCTION PROVIDERS.** Recover analytic core rows and original alpha=1-sigma((y-hb)/hb), chi=1-(1-hb)*sigma(y/hb). Preserve dy=hb*ds and original radial displacement. Integrate comparison and actual histories as coordinate-cell functions with functional joins to the macro; retain atom/width/pressure errors. Cached exits alone do not finish this task.
- [ ] **WHOLE-Z SOURCE/CELL PROVIDER.** Extend the native two-frame factory to original source functions on Z cells using the core recipe, anchored pole primitive and pressure atoms. Preserve source-correlated derivatives, sign/midplane splits and the same family; never extrapolate point rows to whole-Z certification.
- [ ] **RESHAPE/REFERENCE/RESTORATION BINDING.** Feed true R110 functions into the accepted three complete finite terminal kernels, preserving original mean/axial exponential transport, source-owned V110/E, Rsh/Rm tails and separate P0. Acceptance: source-functional downstream inlets and joins.
- [ ] **UNIQUE ACTIVE PATCH.** Solve the original implicit moment repair coefficient functions with Jacobian, uniqueness, derivatives and source-error bounds. Attach actual outer feedback, raw_patch_rows and original inertial compiler to the accepted general fixed-phase Z calculus. No midpoint roots or source-independent coefficient caps.
- [ ] **COMPLETE PATCH AND REMAINING SOURCES.** Integrate active[1,71/40] plus the accepted terminal piece using the original full integrating factor. Bind upstream reference/O2 histories, O2 axial/buffer/O3 and Rc_E/Rc_E_Z, exact original phase and both N-dependent coefficient levels. Acceptance: actual complete all24 source integrals.
- [ ] **FIVE TERMINAL FUNCTIONS AND ONE COMPATIBLE N.** Reconcile full source/repair/cone/interface frequency inequalities with complete-target N^-2, centered limit/tail/Picard contracts and finite N selection. Acceptance: five identities as functions of Z with one compatible finite N, not sampled zeros or rescaled N=1024.
- [ ] **MATCHED BACKGROUND.** Complete annular, flatten, analytic pressure and exact heat joins, axis regularity and finite-energy radial tails. Verify regional admissible stress and flat remainder separately.
- [ ] **REAL n-DEPENDENT RECURSION.** Implement n=1 and n>=2 recovery equations, shared core domain, independent moment repair, divergence-preserving cutoffs and smooth summation. Measure higher-order remainder; coordinate scaling alone is insufficient.
- [ ] **PULSES AND CORRECTED CARTESIAN FIELD.** Implement means and both pulse families, averaged quadratic stress cancellation and independent full corrected NS residual. Export u/v/w over physical x/y/z/t, then measure contraction, relative axial elongation, scale relationships and material winding separately.

Mark DONE only with defining code, scoped report/receipt and a commit. Perform changed-scope checks and then continue production. The whole objective remains active.
