Current successor: [CURRENT_REFERENCE_NEAR_MIDPLANE_MIXED_C1_2026_10_08.md](CURRENT_REFERENCE_NEAR_MIDPLANE_MIXED_C1_2026_10_08.md), checked source 989880da. The four regular mixed-source/primitive/product/averaging tasks below are now DONE on the whole signed central window. This direct-integration evidence remains valid. Read the successor for the next source-cell/overlap/full-Z dependencies.

# Original near-midplane reference: whole signed axial intervals

Checked source [28498a53](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/28498a5314e4ad951a6f73d376c191460d692370). Predecessor: [CURRENT_REFERENCE_MIDPLANE_2026_10_08.md](CURRENT_REFERENCE_MIDPLANE_2026_10_08.md). The existing GPT-5.6 Luna/max worker reviewed the carrier and original pressure provenance read-only. Root implemented, computed and published. No child was spawned; the full reconstruction goal remains active.

## Implemented source and integration

**OriginalReferenceNearMidplane** now encloses the original reference source as a function over a whole signed neighborhood, rather than only at exact Z0. The two computed windows are zeta in[-1e-8,1e-8] and[-1e-6,1e-6]. The second is 100 times wider. Both cover the whole original radial reference interval y in[-5,0].

~~~
Lambda0=Pstar^11*Cstar^10
Z=zeta/Lambda0
Q=1+Z^2
L=1-delta*Z^2 > 0
~~~

Each original Z^k is collected with the radius/Pstar/delta/L powers before conversion. Its exact additional powers are(-11*k,-10*k,0,0,0). The physical Z map remains a nonzero signed formal source; neither tiny Z nor the large derivative carrier is replaced by zero or a native floating-point value. Q and log L use directed positive whole-window hulls. These hulls enclose the defined original functions; they do not preserve exact arithmetic correlation in zeta.

The normalized pressure uses the actual defining alpha enclosure and bounds on the original remainder:

~~~
P0/Pstar^2 = -alpha/Q^2 + R0
(P0/Pstar^2)_Z = Z*(4*alpha/Q^3 + H_R)
H_R = R0_Z/Z, smoothly extended at zero

|R0| <= 5*exp(3/5)/(2*Pstar*Q^2)
|H_R| <= 10*exp(3/5)/(Pstar*Q^3)
|R0_ZZ| <= 44*exp(3/5)/(2*Pstar*Q^2)
~~~

The fourteen actual pressure density derivative factors establish the odd first-jet carrier, including the accepted late and finite-end source budgets. H_R is not identified with R0_ZZ. Finite alpha quadrature error remains in the same defining coefficient enclosure. Every late-error source term retains its one original Pstar^-1 factor.

The jet boxes enclose one defined original remainder and its derivatives, but use separate interval enlargements. They do **not** retain joint arithmetic correlation among R0,H_R,R0_ZZ. Downstream cancellations must use source identities or a stronger shared representation; they cannot assume independent boxes are a correlated field. Ordinary Z rows come from the original source templates, not derivatives of a selected interval bound.

The whole p2 C0 enclosure crosses zero without being identically zero. The genuine p2_Z remains strictly negative and nonzero. The original u=p2*q/dstar is formed directly; the regular Fourier branch handles both signs and zero without p2_Z/p2 or r^-1. On every 16-cell expanded-window cell, |u| <=0.146716, below the original1/4 guard. Exact-Z0 stays compatible with the separate accepted midplane service. The strict signed service's guards remain unchanged.

The unchanged original finite-N coefficient program computes exprel(A/N), exp(A/N), both coefficient orders and all ordinary-Z product rules. Actual phase unions, inverse primitives, positive own-rate masses and radius Jacobian remain. Incoming histories and the separate P0 datum are never reset; pressure zero-rate memory totals5. Source and contribution exports retain the working precision.

Four files with stem **lei_ren_part1_paper_compliant_current_original_reference_near_midplane_integrals** are published under experiments/root_st073: producer, deterministic compressed manifest, focused checker and receipt.

## Computed evidence and scope

Five genuine source levels were computed: narrow window(4 cells,N160),(16,N160),(16,N16384); expanded window(16,N160),(16,N16384). At N16384/16 cells, rounded outward bounds are:

| Row | C0, narrow | C0, 100 times wider | Z/Lambda, narrow | Z/Lambda, wider |
| --- | ---: | ---: | ---: | ---: |
| m |0.000006188|0.000006528|5.672|9.233|
| h |0.000010029|0.000010660|1.601|1.814|
| k |0.000004076|0.000004334|3.999|6.852|
| e |0.000013545|0.000014315|2.069|2.324|
| p |0.000051574|0.000053100|5.966|6.433|

**Lambda=Pstar^11*Cstar^10*L^-2**. The last columns retain this astronomical original positive factor formally; they are not small unscaled ordinary-Z errors. Lambda0 is the coordinate carrier without L^-2. N16384 is a local candidate, not an admitted global frequency. These are whole-window contribution bounds, not terminal five-moment defects or Cartesian momentum residuals.

Producer54.594s, checker29.484s. Focused evidence:82 exact collected source-power identities;14 original pressure derivative factors;30 independent derivatives of original C0 functions with a consistent nonzero even pressure remainder;60 independent polynomial replays and60 corresponding pressure remainder covers;42 independent defining-integral primitive comparisons at negative/zero/positive u;14 accepted midplane root containment comparisons;5 frame/domain rejections;68 actual whole y/zeta cells and25 own-rate mass comparisons. **1208 Git-index dependency hashes PASS**. No ancestor producer or unrelated test suite was rebuilt.

The implementation removes the missing open neighborhood at the midplane for this source window. It does not yet connect that neighborhood to a complete source-signed interval atlas. The next dependency is actual y/yZ regular source and phase averaging, followed by overlap and whole-Z coverage. Full Z[-1,1], all17/24 units, actual five controls, global N, recursion, pulses and corrected NS remain open.

## Ordered executable queue

Read this before historical handoffs. Execute the first unfinished dependency, preserve the source family and actual error contracts, and mark DONE only with code, scoped receipt and commit. Do not repeat accepted checks unless inputs or relevant behavior change. Keep the full original objective and the predecessor detailed queues.

- [x] **FACTORED-NEAR-MIDPLANE-SOURCE:** typed exact zeta window, original physical-Z carrier, polynomial source-power collection, actual C0/Z roots and positive L/Q hulls.
- [x] **ODD-PRESSURE/p2-CARRIERS:** actual fourteen-stage derivative identities, true odd first-jet carrier, finite alpha enclosure, full even errors and ordinary-Z template rows.
- [x] **REGULAR-WHOLE-NEIGHBORHOOD:** original u guard, Fourier primitives and inverse Z derivatives across zero, whole-reference actual-N own-rate contributions.
- [x] **HUNDREDFOLD-COVERAGE:** original same-family source and complete contribution evidence on zeta[-1e-6,1e-6], retaining the larger Z bounds and full branch guard.

- [x] **NEXT REGULAR-y/yZ-SOURCE:** create a separate near-midplane mixed service. For each original term R^r times coefficient, use d_y coefficient = r*coefficient + sum(rate*x*d_x coefficient), with rates(f,H,D,P)=(1/10,1/10,1/5,1/5). Apply this to the original C0 and original Z rows; do not differentiate an interval box, apply L_Z twice, or multiply by Lambda0 twice. Compile the resulting original polynomials with the same exact zeta power collector and pressure error provenance. Acceptance: whole-cell C0,y,Z,yZ roots, independent finite-unit C0 differentiation including full even pressure remainder, and correct exact-midplane u_y=0 while u_yZ remains nonzero.
- [x] **REGULAR-FOURIER-MIXED-PRIMITIVES:** differentiate the original M48 defining Fourier series and directed tails using actual u_y,u_Z,u_yZ. Keep the true-phase implicit denominator1+t^2 and both inverse cross terms; do not call the strict signed r^-1 routine through zero. Acceptance: independent original T1/T2 defining-integral yZ comparisons at both signs/zero, finite tails and source factors, and exact periodic/halfperiod traces.
- [x] **COMPLETE-E/V-MIXED-PRODUCTS:** form B_yZ with E_yZ*T1, E_y*T1_total_Z, E_Z*T1_total_y and E*T1_total_yZ. Keep actual V/E derivatives and original pressure rows in all five densities. Acceptance: complete original product-rule jets with no envelope derivative substituted as a field.
- [x] **WHOLE-NEIGHBORHOOD-C1-AVERAGING:** attach the above mixed primitives to the accepted actual-N coefficient and own-rate IBP programs. Preserve Q_Z=E_Z*A^2*R2(A/N)+E*A*A_Z*exprel(A/N), all nonlinear terms, physical global endpoints, source-valid internal telescoping and pressure zero-rate memory. Compare true direct and averaged enclosures and retain the tighter bound. Acceptance: whole signed-zeta C0/Z contribution receipts retaining formal Lambda, actual N recomputation, separate oracle/error scope and genuine mixed source evidence.
- [ ] **REGULAR/SIGNED-SOURCE-OVERLAP:** extend the source service to positive/negative zeta cells without zero crossing where appropriate. Refine radial/axial cells until the original signed |u| lower guard and the regular upper guard share an actual overlap. Use the same family, pressure error provenance, original phase and positive rho/s/hinv. Acceptance: independent same-source values and ordinary-Z/mixed derivatives contained in both chart representations on an open overlap, no sign inferred from point samples.
- [ ] **WHOLE-Z-REFERENCE-ATLAS:** combine the factored central carrier with source-signed physical-Z intervals to cover[-1,1] as functions. The current physical neighborhood is[-1e-6/Lambda0,1e-6/Lambda0], not[-1e-6,1e-6] in physical Z. Cover intermediate amplification ranges and branch guards; fail closed on uncovered cells. Acceptance: one typed reference-window C1 source/integral API with full interval coverage, actual seam contracts and no uniform-Z sampling substitute.
- [ ] **SOURCE-CORRELATION-WHERE-NEEDED:** if overlap, cancellation or nonlinear mean computation requires joint pressure/L correlation, retain the defined shared remainder/jet relationships and exact L function map in a stronger representation. Acceptance: source-backed correlated operations with independent checks; do not reinterpret the current separate hulls as exact functions.
- [ ] **SYMBOLIC-FREQUENCY/ORACLE-ERROR-CONTRACT:** retain Lambda and other large factors in all-N whole-route C1 constants. Derive admissible frequency conditions from actual complete source bounds, compatible with the accepted centered all-N/unit-C1-ball bridge. Keep original radius phase, endpoint provenance, numerical/oracle error and Picard tail separate. Acceptance: executable typed contract without materializing an impossible native exponential; bounded160/16384 do not establish global N.
- [ ] **SIGNED-NONLINEAR-MEANS/DRIFT:** integrate actual quadratic/exponential means and their controlled y/Z drift when signed terminal closure requires them. Reflection of the leading density does not cancel nonlinear means. Preserve original phase measure, Jacobian and finite-N error.
- [ ] **O2/AXIAL/11-UNIT-BUFFER:** implement actual C0,y,Z,yZ source rows with original inherited pressure/moment histories, nonconstant parameters and all offsets. Integrate genuine densities with explicit seam and initial-data contracts; do not reset incoming corrections.
- [ ] **O3/RESTORE/ALL17-24-SOURCE-ORACLE:** follow the predecessor unit queue, recover every remaining definition and history, and publish complete same-N C1 function/integral outputs. Unsupported units fail closed.
- [ ] **ACTUAL-FIVE-CONTROLS/TERMINAL-CLOSURE/PICARD:** connect the complete actual oracle to the accepted centered bridge, solve the five terminal identities as Z functions, admit N and replay corrections with separate oracle/Picard errors.
- [ ] **GLOBAL-JOINS/EXACT-HEAT/FINITE-ENERGY:** assemble consistent core, annulus, flatten, collar and exact heat exterior with analytic pressure compatibility, all required joins and radial tails.
- [ ] **ADMISSIBLE-STRESS/FLAT-REMAINDER:** compute the actual background residual as -div(T_B)+E_B, establish cone margins across every original region and genuine flat decay separately from stress.
- [ ] **ACTUAL-n-RECURSION:** execute distinct n=1 and n>=2 recovery equations, original common core interval, separate moment repairs, divergence-preserving cutoffs and controlled summation. A coordinate-scaled old field does not complete this task.
- [ ] **TWO-PULSE-FAMILIES/CORRECTED-UVW:** implement both oscillatory families and averaged quadratic stress cancellation, recover the corrected Cartesian u/v/w, independently validate forced NS residual, and measure physical core contraction, axial aspect ratio and material winding.

Whole-Z terminal closure and the full reconstruction objective remain active.
