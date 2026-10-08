# Original slope joint source and separated Rc error ledger

Checked source [79746360](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/79746360a9d55e907b76019e4bfd6a58426743aa). Predecessor: [CURRENT_RC_JOINT_TERMINAL_DEFECTS_2026_10_08.md](CURRENT_RC_JOINT_TERMINAL_DEFECTS_2026_10_08.md). Reused worker **GPT-5.6 Luna / max**, read-only formula review. Root implemented, computed, checked and published; no new child or child compute/Git task.

## Implemented result

At the original N1024 on Z[.36,.38] and [-.38,-.36], all **4096 slope source cells and 6148 saved phase queries** now have joint k-a_eff*m and ordinary Z derivative covers. The adapter restores each original native logPstar/delta/L/log|u|/R basis, signed roots, positive q and conditioned source branch. It evaluates the original primitive formula on the saved inverse coordinate; no inverse solver or ancestor producer is rerun. The original saved derivative inverse boxes provide the primitive Z covers.

Nonlinear joint density is formed before phase/source unions, integrated with the original rate-3/2 endpoint mass and transported through the original axial, buffer, transition and quiet power widths to Rc. The old composite slope inlet is now split into **pre-slope inherited memory** and **slope own source**, alongside the three accepted post-slope source contributions. The five correction-only M,D=(J-M)/mu,I,S,Cp C0/Z covers, N*r and first original linear control responses have been regenerated. Other four component histories retain their accepted original covers. P0 and rate-zero pressure memory are unchanged.

This completes the slope-own joint adapter. It does **not** recover the lost pre-slope root/phase correlation, original point/function controls, whole Z/axis, legal common N, nonlinear five-moment identities, exact heat/stress, true coefficient recursion or corrected uvw. These gates remain false.

## Exact slope amplitude and arithmetic

The original slope has

    E(y,Z)=qi*exp(y/10-(3/5)*J_sigma(y)), qi=1/(1+Z^2).
    d(y)=logPstar+3-y=exp(Md)+14-y.
    a_eff(y)=A_Rc*exp(d(y)/2)
            =qi*exp(3/10-y/2-5*mu/2).

The post-slope distance is logPstar+2=exp(Md)+13. It includes the original axial width exp(Md)-1 plus 11 buffer, one transition and two power units. The slope contributes one additional unit. Dropping either quiet power unit changes the effective amplitude and is forbidden.

The axial/buffer identity a_eff=E*exp(-5mu/2) is only valid at the slope endpoint, not over its interior. The slope uses

    dE=E*expm1(Aprimitive/N), dV=Bprimitive/N,
    C=V*dE+dE*dV+(E-a_eff)*dV,
    C_Z=V_Z*dE+V*dE_Z+dE_Z*dV+dE*dV_Z
        +(E_Z-a_eff_Z)*dV+(E-a_eff)*dV_Z.

The nonzero original expm1(-5mu/2) increment stays separate from the moderate slope mismatch. A bounded-log version of the same identity expm1(x)=x*integral_0^1 exp(t*x)dt retains the microscopic native log|u| factor near the slope exit. It changes arithmetic, not the source profile. It rejects source exponent covers outside [-1,1].

Native powers (p,d,l,u,r) collect exactly to half-Pstar power 2*(p-4d+10r), with offset

    original_offset-30d+l*logL+u*log|u|+r*(log110+10logC+y).

Neither huge positive source amplitudes nor tiny nonzero roots are replaced by point values. The original signed source records and true phase bindings remain hash-bound.

For inherited correction values at y=0,

    C_in=k_in-a_eff(0)*m_in,
    C_in,Z=k_in,Z-a_eff,Z(0)*m_in-a_eff(0)*m_in,Z.
    C_slope,out=exp(-3/2)*C_in+C_slope,own.

The inherited expression currently subtracts component covers. This is conservative and explicitly correlation-lost. The new slope own term replaces the old composite inlet contribution when composing Rc; it is not added to that composite a second time. Original axial/buffer/transition joint contributions are reused from the checked predecessor archives.

## What the new ledger establishes

The complete original pre-slope route retains twelve active charts and the source-defined zero correction at its initial flat inlet. Each chart's C0/Z contribution is carried to the O2 inlet and separately to Rc. The refined inner_reference pressure integral is retained; no original chart or pressure constant is discarded.

On both tiles, **active_first_bridge supplies the largest inherited Z cap in all five component rows**. The new separated slope-own divided-D Z cap has natural logarithm about 8.238e17; the inherited divided-D Z cap has natural logarithm about 5.937e+408906090034569676. The extreme total derivative bound is therefore inherited before the slope. This identifies a concrete next chart; it does not prove that its actual derivative has this size.

The total leading C0/Z caps are not materially smaller than the predecessor's. In particular, the slope-own divided-D C0 cap also has leading logarithm 8.238e17; splitting the source is not evidence that it is already small. Original N1024 plus current covers still cannot certify the nonlinear small repair. The pressure row also retains a large rate-zero memory cover. The main advance is a complete source-owned slope calculation and a correctly separated obstruction, not a claim of moment closure or new field accuracy.

The first bridge's active-kappa narrowing already uses the original identity Delta_Z=a_Z*(1-b^2/a^2)+2*b*b_Z/a. Its remaining wide bounds include independent products, phase inverse factors, branch hulls and retained linear q_Z/p2_Z terms. Do not rerun this already accepted narrowing and call it new work. The saved pre-slope serial contributions lack the phase-held tuples needed for a fully correlated joint evaluator; recover those from the original source owner, with a bounded new adapter if necessary.

## Evidence

Six files under experiments/root_st073 use stem **lei_ren_part1_paper_compliant_current_slope_joint_source**: producer/manifest, two signed gzip archives, checker/receipt. Archives total **46,886,686 compressed bytes** and retain complete lossless cell evidence.

Focused checks passed: 48 independent slope joint C0/Z and effective-amplitude comparisons; three independent exact different-rate integrals; three native-to-half power collection checks; tiny native expm1 and nonzero mu checks. All 4096 saved cell contributions are summed, 16 original source cells are independently replayed from the accepted inverse/source chunks, and both pre-slope attribution/Rc target compositions are replayed. **1308 Git-index dependency hashes PASS.** Producer 98.203s; checker 15.110s. Accepted upstream producers and inverse solvers were not rerun.

## Executable next tasks

This handoff supersedes the predecessor's primary slope-own task. Mark DONE only after committed implementation and necessary focused evidence. Keep function enclosures, actual function evaluations and certified identities distinct. Reuse accepted receipts; do not repeat whole ancestor producers or suites merely to reconfirm unchanged files.

- [x] **SLOPE-OWN-JOINT-SOURCE:** complete 2048 cells per tile, original roots/q/phase/inverse, exact slope a_eff, joint C/C_Z before hulls, original masses and native basis conversion.
- [x] **SEPARATED-RC-SOURCE-LEDGER:** pre-slope inherited, slope own, axial, buffer and transition; corrected total targets and original linear responses; unchanged pressure memory/P0.
- [x] **FIRST-INHERITED-C1-CHART-ATTRIBUTION:** all twelve source contributions; active_first_bridge is the largest inherited derivative cap on both tiles. Root-level dependency classification remains open below.
- [ ] **FIRSTBRIDGE-PHASE-HELD-SOURCE-ADAPTER (primary):** own only active_first_bridge and its export adapter. Consume the existing source family, true geometry, signed root jets, same original a/b/p2/q and cutoff branch support. Save same-source phase-held primitive C0/Z tuples and actual inverse provenance on the two strict-sign tiles. Do not select branch endpoints as field values or replace q_Z/p2_Z by zero. If live original source queries are required because old archives lack tuples, execute only this chart, without recomputing the other eleven charts or full downstream route.
- [ ] **FIRSTBRIDGE-C1-TERM-LEDGER:** record the individual a_Z, b_Z, p2_Z, q_Z, inverse-coordinate and primitive product terms before unions. Find the first operation producing the extreme cap. Compare source root bounds, reciprocal/branch conditioning and unit conversion separately. Distinguish dependency loss from a proven original derivative bound. Retain all linear terms and derivative inverse contracts.
- [ ] **FIRSTBRIDGE-JOINT-K-M-INTEGRAL:** derive local a_eff from exact source-to-Rc log-radius distance and common A_Rc. Form k-a_eff*m and its Z derivative before products/phase hulls; use the original rate-3/2 Duhamel mass. Export a joint source contribution and replace only this chart's inherited component approximation. Replay all twelve chart transport contributions without discarding any. Compare old/new divided-D widths, with unchanged m,h,e,p/P0.
- [ ] **SLOPE-J-MISMATCH-CORRELATION:** bind the original cached J_sigma cell source to E and a_eff. Where helpful use exp(.3-.6*y+.6*J_sigma(y)-2.5*mu)-1 as the relative mismatch instead of independent E/base interval subtraction. Enclose expm1 on the actual moderate slope range; keep the tiny mu increment separate. This must remain a source identity, not a fit or endpoint cancellation assumption.
- [ ] **SLOPE-PERIOD-CANCELLATION:** derive the same-source frozen mean plus slow/Duhamel-weight errors for nonconstant E,V,q. Include epsilon_y=(-.6+.6*sigma)*(1+epsilon) and its extra -E*epsilon_y*dV term in a slow joint derivative. Do not copy buffer C_y unchanged. Derive mixed yZ from the original implicit inverse before claiming sharp Z averaging. Start with the limiting slope cells; direct signed integral covers remain the baseline.
- [ ] **PRESSURE-MEMORY-CANCELLATION:** separate first-order oscillatory slope pressure terms from their quadratic mean; integrate with rate zero. Tighten the original firstbridge/other inherited pressure accumulation if it dominates. Keep the shared P0 and preheat datum operator separate; no artificial pressure tail or reset of the inherited constant.
- [ ] **N-DEPENDENT-ORIGINAL-ORACLE:** expose original r(N,Z) through function-valued source/integral callbacks or valid uniform-N coefficient/remainder bounds. Actual phase, inverses, controls and all integrals must consume the identical N. The two N1024 receipts are not a legal all-N admission.
- [ ] **COMMON-N-AND-NONLINEAR-CONTROLS:** after source C1 bounds are improved, use the exact original divided axial/swirl matrix and Q(mu,h)/N to derive one compatible N and a controlled nonlinear solve for axial0/axial2/swirl0/swirl1/swirl2. Keep the mu^-1 cross term, energy/pressure quadratic terms and quotient derivatives. Current linear responses are not repaired control functions.
- [ ] **REPAIRED-VELOCITY-AND-TERMINAL-FUNCTIONS:** install the original F/G bumps on Rc..2Rc with derivative controls; compute actual five cumulative moments and terminal identities as Z functions. Check joins/cone changes caused by the actual repair. An invertible matrix or interval target recipe alone is insufficient.
- [ ] **AXIS/WHOLE-Z/DATUM/EXACT-HEAT:** supply the Z=0 degeneration, edge/high-jet charts and common N/controls; bridge the exact current P0 operator to raw analytic preheat normalization. Match the legal finite-energy exact heat exterior with the same five identities, not arbitrary constants.
- [ ] **ADMISSIBLE-STRESS/FLAT-REMAINDER:** reconstruct residual=-div(T_B)+E_B, signed cone margins over all joins/collars and independently bounded high-order flat decay. Preserve the distinction between background residual and corrected full NS residual.
- [ ] **TRUE-n-DEPENDENT-RECURSION:** implement distinct n=1 and n>=2 recovery equations, independent five-moment repair per order, common core interval and a controlled smooth sum. Apply truncation before curl. Source jets and coordinate rescaling are not coefficient recursion.
- [ ] **TWO-FAMILY-PULSES/CORRECTED-UVW:** against admitted stress, implement both oscillatory pulse families and mean corrections, test averaged quadratic cancellation, export physical u,v,w and independently measure NS/divergence/energy, core contraction/elongation and material winding. Animation remains secondary.

The next production pass should recover and integrate the **first bridge's correlated phase-held derivative source**, then replace that chart's current inherited approximation in this reusable terminal composition. Additional buffer precision does not address the present derivative obstruction.
