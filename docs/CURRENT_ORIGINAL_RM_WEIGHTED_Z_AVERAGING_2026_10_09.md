# Successor: genuine long-window finite-N source production (2026-10-09)

[Current full R110-to-Rsh source and local Duhamel](CURRENT_ORIGINAL_LONG_RESHAPE_FINITE_N_2026_10_09.md) implement the first conditional RM-W4b.1/.2/.3 layer, including nonzero nonlinear density drivers and explicit incoming memory. Real current R110/Rm inlet, downstream source windows, functional closure, cone/global N, whole-axis control and coefficient recursion remain OPEN. Full reconstruction **ACTIVE / INCOMPLETE**. Historical Rm averaging evidence below is unchanged.

---

# Actual Rm genuine axial phase primitives and weighted Z integral bounds

Checked source [43f79b63](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/43f79b638a2a131a8baa487aa2a93aa054b4e899). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [genuine mixed source/primitive yZ](CURRENT_ORIGINAL_RM_MIXED_YZ_PRIMITIVES_2026_10_09.md).

RM-W3 now has its first complete conditional local implementation: genuine f1_Z/f1_yZ, integral-defined G_Z/G_yZ, full nonlinear second-remainder Z derivatives, actual radius endpoint phases, and own-rate endpoint-retaining N^-2 averaging. All ten local Z absolute integral bounds over x in [1,e] at Z=0,.5 strictly improve relative to the accepted complete direct baseline at candidate N=257. Existing ten C0 bounds are bitwise unchanged. **These are bounds, not measured defect values, five-moment closure or scale-recursion progress percentages.**

## Actual reduction routes

| Density | Z=0 whole bound tighter | Z=0 cells using IBP | Z=.5 whole bound tighter | Z=.5 cells using IBP |
| --- | --- | ---: | --- | ---: |
| m | yes | 4 / 10 | yes | 0 / 10 |
| h | yes | 4 / 10 | yes | 0 / 10 |
| k | yes | 4 / 10 | yes | 0 / 10 |
| e | yes | 0 / 10 | yes | 0 / 10 |
| p | yes | 4 / 10 | yes | 0 / 10 |

At Z=0, 16 cell/density rows select the new IBP bound. At Z=.5, the direct bounds on the new genuine-source partition are tighter than conservative mixed IBP, so the direct route is used. The Z=.5 improvement must not be attributed to oscillatory cancellation. Every cell chooses a valid direct or averaged cover; the sum then competes with the accepted whole-patch baseline.

The old and new absolute logarithmic upper bounds are still astronomically large. Their full exact interval tuples are retained in the report; no floating-point conversion of these logarithms to infinity or a rounded100% reduction is used as a scientific result. Tighter covers remove overestimation, but do not establish small defects, numerical convergence, a positive stress-cone margin or an admitted frequency.

## Original axial phase primitive service

Use the same actual mixed4/generic/quotient source, P0, formal Rm/radius, bases and ledger. Original A/B reflection gives zero uniform phase mean of the signed leading vector. Differentiating these original functions at fixed phi gives

    f1=(B, E*A, E*(V*A+B), 2*V*B-E^2*A, E^2*A)
    G_Z(y,Z,phi)=integral_0^phi f1_Z(y,Z,s) ds
    G_yZ(y,Z,phi)=integral_0^phi f1_yZ_fixed_phi(y,Z,s) ds
    G_Z(0)=G_Z(1)=G_yZ(0)=G_yZ(1)=0.

B denotes B/Pstar normalized exactly once. The genuine source yZ retains correlated shear/q, full signed pressure and inertial terms, and p2_yZ=R*(p2bar_yZ+p2bar_Z). First y and Z rows alone are not mixed data. The complementary original integrals yield

    |G_Z(phi)| <=min(phi,1-phi)*sup_phase|f1_Z|
    |G_yZ(phi)|<=min(phi,1-phi)*sup_phase|f1_yZ|.

phase_primitive_Z_enclosure returns these integral-defined function enclosures, not selected antiderivative values or derivatives of selected caps. At closed endpoints they vanish exactly. The accepted general conditioned first graph can tighten its eight original C0/y/Z/phi rows; genuine yZ remains supplied by the checked all-signed-u mixed theorem. Choosing a tighter valid derivative enclosure does not differentiate that choice.

## Retained nonlinear axial bias

Retain deltaE=E*A/N+R_E/N^2, R_E=E*A^2*R2(A/N), deltaV=B/N and F_N=N*deltaE. For absolute source caps and expcap=exp(|A|/N),

    |R_E|<=|E|*|A|^2*expcap/2
    |R_E_Z|<=(|E_Z|*|A|^2+2|E|*|A|*|A_Z|)*expcap/2
               +|E|*|A|^2*|A_Z|*expcap/(6N)
    |F_N|<=|E|*|A|*expcap
    F_N_Z=N*E_Z*expm1(A/N)+E*exp(A/N)*A_Z
    |F_N_Z|<=(|E_Z|*|A|+|E|*|A_Z|)*expcap.

The masses1/2 and1/6 are integral(1-s) and integral s(1-s). The F_N formula includes the exprel derivative through exprel(x)+x*exprel'(x)=exp(x). The exact signed second remainder vector is

    (0, R_E, V*R_E+F_N*B,
     -E*R_E+B^2-F_N^2/2, E*R_E+F_N^2/2).

Every Z product term is included before taking a positive majorant, including V_Z*R_E, E_Z*R_E, F_N_Z*B and F_N*B_Z. Neither the nonlinear mean nor its Z derivative is zeroed. Only the proved bounded |A|/N scalar is exponentiated; actual physical positive factors remain formal.

## Original weights, endpoints and candidate binding

For y=log(x), dy=dx/x, rates m1,h3/2,k3/2,e1,p0 and K=exp(-lambda*(1-y)), the actual Rm map has phi_y=N and phi_Z=0. Thus

    integral K*f1_Z/N dy
      =([K*G_Z]_left^right-integral K*(G_yZ+lambda*G_Z)dy)/N^2.

The implemented bound keeps both actual endpoint phase boxes, every own-rate positive mass and downstream suffix, and the nonlinear remainder derivative. Pressure has rate0 and exact decays1. No inter-cell endpoint cancellation is claimed. Source identity, P0, Rm and exact phase-Z independence are guarded. The baseline comparison rejects a candidate N different from its checked N257; this is frequency consistency, not global N admission.

The affine terminal operator still requires

    deltaH_j_Z(Rh)=exp(-lambda_j)*deltaH_j_Z(Rm)+local_integral_Z([1,e]).

True finite-N incoming corrections at Rm are **unsupplied**. Leading histories, None and arbitrary zeros must never replace them. Total spatial mixed yZ needs A_phiZ/B_phiZ and remains OPEN; it is not required for this fixed-phase slow mixed IBP because phi_Z=0.

## Scoped evidence

Artifacts: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_weighted_Z_averaging.py, .json.gz, _check.py, _check.json. Final producer 83.875s; checker 142.844s; both terminal exit0. **1211 staged exact dependency hashes PASS**. Existing GPT-5.6 Luna/max reviewed calculus, ownership and scope read-only. Root owns code, compute and Git. No ancestor producer/checker or old-owner contribution was run.

Independent symbolic checks establish five exact signed density-Z splits, all second-remainder products, F_N_Z/exprel identity, both integral masses and the endpoint-retaining Z IBP. 27 original exponential/remainder comparisons use negative/zero/positive A and nonzero E_Z,A_Z,B_Z,V_Z. 30 original GenericShearLoop comparisons use both signed p2, nonzero slow/mixed E/V, original angle-to-phase Jacobian and nonvertex endpoints. Quadrature precision90/allowance1e-65 is a diagnostic; exact calculus plus live source bounds establish conditional covers.

Live/report agreement covers 20 cells, 400 leading C0/y/Z/yZ rows, 400 endpoint primitive rows and all ten strict Z reductions. All ten C0 bounds are bitwise preserved; wrong phase/source/endpoints and cross-N comparisons are rejected. All full reconstruction gates remain false.

## Concrete prefix starting points, still unbound

The read-only inventory identifies the existing [native bridge/switch C1 history operator](../experiments/root_st073/lei_ren_part1_paper_compliant_current_native_bridge_switch_C1_histories.py) as a reference for original finite-N correction histories at R110. Its saved candidate is N1024, not the current N257. Its owner/function provenance must be compared against the current original owner; matching family metadata is insufficient. Do not scale or rebase its saved values into a claimed current inlet.

The next source windows are implemented in [original long-reshape endpoint](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.py) and [original reference/restoration functions](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reference_restore_functions.py). These provide background affine history kernels. Their homogeneous decay parts alone do not supply an entire finite-N correction: every nonzero candidate-background density source must also be integrated, or proved identically zero for that window. The pressure/phi normalization uses the genuine background endpoint fields, not a phi reconstructed from correction histories.

The [original Rm generic inputs](../experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs.py) provide the current exact point/cell E,V, first ordinary radial rows, five histories and P0. A typed adapter must prove current function/parameter/P0/ledger identity across native reduced and full five-log bases before using earlier arithmetic rows. Pure rebasing changes representation and is not that proof. Avoid cascading old physical-assembly constructors; consume accepted raw source data and explicit kernels under the current owner.

## Detailed next production tasks

- [x] **RM-W3a SOURCE yZ:** same-owner correlated shear/Delta/q/pressure recovery with one radius factor and correct source-order consumption.
- [x] **RM-W3b FIXED-PHI A/B/LEADING yZ:** general original all-signed-u mixed implicit/Fourier enclosures, complete products and endpoint terms.
- [x] **RM-W3c CONDITIONAL G_Z/G_yZ / Z IBP:** genuine original integral functions, endpoints/own rates/nonlinear Z bias, valid direct fallback, same-N comparison; ten local Z bounds strictly improve. Sharper caps and whole-axis extension remain tasks below.
- [ ] **RM-W4 PRIORITY: REAL FINITE-N PREFIX AND Rm INLET.** Start from an actual source window, compute common E,V and five histories/P0 with original recovery, and bind actual radius phase. Restore the loop's q branches/parameters for that same owner; compile both N levels of original candidate densities. Integrate every source window with original dy=dr/r and own-rate tails, including true transition endpoints. Carry corrections continuously into Rm and export C0/Z incoming vectors. Do not reset history or transplant an old-owner receipt. Add a source-bound affine application that refuses unsupplied corrections.
- [ ] **RM-W4a CURRENT SOURCE / FREQUENCY ADAPTER.** Read the saved bridge/switch operator as a definition reference, inventory exact current R110 endpoint source objects, and prove datum/parameter/function identity before conversion into the current full basis. Build fresh current-owner N257 C0/Z source corrections; N1024 saved data are not the N257 inlet. Reject missing identity, mismatched N, wrong P0 or reduced-basis transplantation.
- [ ] **RM-W4b FIRST DOWNSTREAM WINDOW.** Implement R110-to-Rsh correction transport with exact history conversion, all own-rate decays, and genuine candidate-background forcing. Background phi and phi_Z must come from the current endpoint; do not normalize correction rows as total histories. Prove any zero driver rather than assuming homogeneous transport.
- [ ] **RM-W4c REFERENCE / RESTORATION / RM BOUNDARY.** Continue the affine correction through each actual source window, compute retained driver integrals and endpoint matching, and bind the final vector to this current Rm source. Expose values plus source-order/provenance/matching evidence; leave complete inlet admission false until every prefix window is included.
- [ ] **RM-W5 SHARPER MIXED/PRIMITIVE COVERS.** Improve all-u mixed caps through valid conditioned original inverse derivatives, slow-source correlation and justified subdivisions. Demonstrate integral-width reductions; do not treat astronomical upper-bound tightening as convergence. Prove exact shared-source endpoint identity before cancelling cell boundaries.
- [ ] **RM-W6 COMPLETE Rc / ALL24.** Feed the genuine prefix and current local C0/Z contribution through actual outer windows with all suffixes and both original N levels. Compute Rc_E/Rc_E_Z and terminal defects as functions under the current owner.
- [ ] **RM-W7 WHOLE-Z / SPATIAL MIXED / HIGHER.** Generalize beyond conditional0,.5, prove cone/cutoff/pole/midplane branches, and add A_phiZ/B_phiZ, yy/ZZ and further genuine raw orders with all fast N chains.
- [ ] **RM-W8 CONE / FUNCTIONAL FIVE-MOMENT REPAIR / ONE N.** Prove actual owner margins, solve complete B*h+N*r+Q/N from genuine cumulative inputs, uniqueness and all five terminal identities throughout Z. Candidate257 is not frequency admission.
- [ ] **RM-W9 PRESSURE / HEAT / FINITE ENERGY.** Compatible functional preheat pressure, exact heat exterior joins and controlled finite-energy radial tails.
- [ ] **RM-W10 STRESS / FLAT / REAL n-DEPENDENT RECURSION.** Admissible divergence-form stress and flat remainder, distinct n=1/n>=2 recovery, independent repairs, curl-based truncation and controlled smooth summation. Scale recursion is still OPEN.
- [ ] **RM-W11 PULSES / FULL NS / DYNAMICS.** Original mean/oscillatory correction and averaged stress cancellation, smooth forcing and independent Cartesian residuals; quantitative contraction, relative elongation, accumulated material winding and finite energy.

Mark DONE with code, scoped report/receipt and commit. Full goal remains active. Proceed to source production; no ancestor reruns are required for this local layer.
