# Actual Rm phase primitives and weighted C0 integral reduction

Checked source [85a24450](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/85a24450c1ca121c8d0f66224bcd9db9ec4261b6). Full reconstruction **ACTIVE / INCOMPLETE**. Predecessor: [genuine first spatial source/density jets](CURRENT_ORIGINAL_RM_FIRST_SPATIAL_JETS_2026_10_09.md).

All five conditional Z=0 local C0 integral absolute bounds over x in [1,e] are now strictly smaller than the accepted complete direct bounds. The reductions are approximately40%–67%. The conditional Z=.5 baseline remains the tighter cover, so it is retained. Ordinary Z integral bounds remain bitwise unchanged: genuine mixed yZ rows are still needed for their oscillatory reduction. Real finite-N incoming histories, functional five-moment closure, actual cone/global N, heat/stress and true coefficient recursion remain OPEN.

## Measured local C0 improvement

These ratios compare proved absolute bounds of the same original local source integral at candidate N=257, not numerical defect values or whole-project progress percentages.

| Density | New / previous absolute bound | Bound reduction | Previous / new |
| --- | ---: | ---: | ---: |
| m | 0.4122 | 58.78% | 2.426 |
| h | 0.4366 | 56.34% | 2.290 |
| k | 0.3304 | 66.96% | 3.027 |
| e | 0.4839 | 51.61% | 2.066 |
| p | 0.6035 | 39.65% | 1.657 |

Twenty source cells cover both conditional frames. The new operator chooses a valid IBP or direct cap cell by cell, then compares the sum against the accepted complete local bound. 20 cell/density rows use the IBP cap. Bounds depend on the original source and can remain large; improvement does not prove small defects or exact terminal moments.

## Original phase primitive service

For the same actual leading Rm source, original A and B/Pstar are odd under phi ->1-phi. The general Section11 reflection identity uses t0=-b/a and is valid when t0 is nonzero. E and V are fixed with respect to phi. Thus the original signed leading vector has exact zero uniform phase mean:

    f1=(B, E*A, E*(V*A+B), 2*V*B-E^2*A, E^2*A).

Here B denotes the original normalized B/Pstar; no second unit conversion is applied. Source-bound original phase primitives are

    G_j(y,Z,phi)=integral_0^phi f1_j(y,Z,s) ds
    G_j_y_slow=integral_0^phi f1_j_y_slow(y,Z,s) ds
    G_j(0)=G_j(1)=G_j_y_slow(0)=G_j_y_slow(1)=0.

The new phase_primitive_enclosure API encloses these actual integral-defined functions and their genuine fixed-phi y derivatives. Complementary integrals give

    |G_j(phi)| <= min(phi,1-phi)*sup_phase|f1_j|
    |G_j_y_slow(phi)| <= min(phi,1-phi)*sup_phase|f1_j_y_slow|.

The actual endpoint radius phase boxes are used to bound the distance to0/1, including modulo seam unions. Exact endpoints return zero. These outputs are function enclosures, not selected G values or derivatives of an absolute cap. No separate narrow antiderivative value or cancellation between independent cell endpoint boxes is claimed.

Full-period primitive caps use the same live mixed4/generic/quotient objects, real first-y source roots and the original conditioned first graph or all-u theorem. Original E,V,P0,Rm, geometry, canonical bases and ledger are retained. The full-phase y derivative is the slow derivative; the total N*phi spatial derivative must not be used in the IBP slow term.

## Signed nonlinear density split

No zero-mean assumption is made for a complete finite-N density. Retain the exact exponential remainder and both original N levels:

    deltaE=E*A/N+R_E/N^2
    R_E=E*A^2*R2(A/N), R2(x)=integral_0^1 (1-s)*exp(s*x) ds
    deltaV=B/N
    F_N=N*deltaE=E*A*integral_0^1 exp(s*A/N) ds
    deltaDensity_j=f1_j/N+remainder_j/N^2.

The exact signed second vector is

    (0, R_E, V*R_E+F_N*B,
     -E*R_E+B^2-F_N^2/2, E*R_E+F_N^2/2).

Its caps keep all linear, cross, quadratic and pressure terms. R2 is bounded by exp(|A|/N)/2; F_N is bounded by |E*A|*exp(|A|/N). Only the proved bounded A/N scalar is exponentiated. The actual enormous positive radius/Pstar factors remain formal. The independent original frozen fixture confirms that nonlinear pressure phase mean is strictly positive while the leading vector mean is zero.

## Endpoint-retaining own-rate integration

Use true y=log(x), dy=dx/x, original own rates m1,h3/2,k3/2,e1,p0. For a cell a<=y<=b and downstream target log(x)=1,

    K(y)=exp(-lambda*(1-y)), K_y=lambda*K
    phase=frac(N*y+phi0), phase_y=N, phase_Z=0
    integral_a^b K*f1/N dy
      =([K*G]_a^b-integral_a^b K*(G_y_slow+lambda*G)dy)/N^2.

Each cell keeps both actual endpoint terms. For width w=b-a, local decay d=exp(-lambda*w), suffix decay d_s and positive mass M=integral_0^w exp(-lambda*(w-s))ds, the implemented cap is

    d_s/N^2 * (sup|f1|*(distance_right+d*distance_left)
               + M*(sup|f1_y_slow|/2+lambda*sup|f1|/2+sup|remainder|)).

Pressure rate0 has M=w and both decays exactly1. No extra R,Pstar,N or cell-width multiplier is applied. The exact partition preserves the six active support edges and terminal suffix cells. Actual incoming correction is still separate:

    deltaH_j(Rh)=exp(-lambda_j)*deltaH_j(Rm)+local_integral_j([1,e]).

Nothing in the averaging operator supplies or resets deltaH_j(Rm).

## Scoped evidence

Artifacts: experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_weighted_averaging.py, .json.gz, _check.py, _check.json. Producer78.625s; checker127.594s; both terminal exit0. **1203 staged exact dependency hashes PASS**. Existing GPT-5.6 Luna/max reviewed the weighting, endpoint, signed split and scope read-only; root owns code, compute and Git.

The checker independently proves all five exact signed leading/remainder algebraic splits, the general nonzero-t0 reflection and weighted IBP sign/N^-2 identity. Two original frozen signed-p2 GenericShearLoop cases give30 independent weighted-integral, phase-antiderivative and zero-mean comparisons through the original angle-to-phase Jacobian. They retain nonvertex endpoints and a nonzero nonlinear pressure mean. Quadrature precision90/allowance1e-70 is a finite diagnostic; exact identities and live source bounds establish the conditional local enclosures.

Live/report agreement covers20 original source cells, 200 genuine leading C0/slow-y rows and400 endpoint primitive rows. All five whole-patch C0 reductions at Z=0 are checked; Z=.5 preserves the tighter baseline. All10 ordinary Z rows are bitwise unchanged. No missing global or recursive gate is admitted.

## Detailed next production tasks

- [x] **RM-W1 ACTUAL ORIGINAL PHASE PRIMITIVES.** Bind the general A/B reflection to the same actual source; implement G and genuine G_y_slow enclosures through exact integral definitions, full-phase leading pairs and complementary endpoint distances. Keep phase unions and exact endpoint zeros.
- [x] **RM-W2 OWN-RATE C0 AVERAGING / RETAINED NONLINEAR BIAS.** Apply exact endpoint-retaining IBP with true dx/x, all own-rate masses and suffixes, both original N levels, all cross/quadratic/pressure terms. Choose valid direct/averaged covers; achieve five conditional Z=0 reductions without overwriting tighter Z=.5 or Z bounds.
- [ ] **RM-W3 PRIORITY: GENUINE MIXED yZ AND Z AVERAGING.** Extend actual source recovery, correlated a/q/p2 and original implicit inverse to genuine mixed yZ; include the radius derivative once and account for consumed axial order. Compute f1_Z and f1_yZ, then G_Z/G_yZ and the Z IBP. Preserve all existing Z covers until improved ones are independently valid; separate first y/Z rows do not supply yZ.
- [ ] **RM-W4 PRIORITY: REAL FINITE-N Rm INLET.** Reuse actual core/micro/macro/switch/long-reshape/reference/restoration functions and attach the original phase/density compiler to each physical source window. Carry actual corrections and every own-rate tail into Rm; keep P0 separate. Supply the true C0/Z incoming vector to the local affine operator. None, leading histories and arbitrary zeros are not this vector.
- [ ] **RM-W5 SHARPER PHASE ANTIDERIVATIVE / CORRELATED CELLS.** Improve G/G_y endpoint enclosures using the actual original phase antiderivative or justified integral subdivisions; keep common flat beta carriers and source correlation. Prove and implement internal endpoint cancellation only for the identical source function/phase across adjoining cells. Report actual integral width reductions rather than just smaller primitive caps.
- [ ] **RM-W6 COMPLETE Rc / ALL24 / FINITE-N REPAIR INPUTS.** Combine the true finite-N prefix, local Rm contribution and actual outer windows with every original suffix and both levels. Produce Rc_E/Rc_E_Z and terminal defects as functions. Avoid reset at Rm/Rh or inheritance from an older patch owner merely because source-family metadata matches.
- [ ] **RM-W7 HIGHER MIXED JETS / CUTOFF / WHOLE-Z.** Restore yy/yZ/ZZ and further source orders using exact raw data and recovery equations; prove active/flat/cutoff branches and midplane/pole behavior. Retain N*phi and N^2*phiphi chains for actual spatial derivatives. Conditional0,.5 frames are not a whole-axis provider.
- [ ] **RM-W8 ACTUAL CONE / UNIQUE REPAIR / ONE GLOBAL N.** Prove same-owner H0,D,J relaxed-cone margins and selected parameter admission; solve B*h+N*r+Q/N from complete cumulative sources, then prove all five terminal moment identities throughout Z. Candidate257 and improved local ranges do not constitute global frequency admission.
- [ ] **RM-W9 ANALYTIC PRESSURE / HEAT / ENERGY.** Construct the compatible functional whole-axis preheat pressure datum, exact heat joins and controlled radial tails with finite physical energy. Keep pressure data coupled to the core and matching; no independent fitted residual-canceling tail.
- [ ] **RM-W10 ADMISSIBLE STRESS / FLAT / REAL n-DEPENDENT RECOVERY.** With completed matching, build divergence-form stress and flat remainder, prove cone/norm/scale bounds and implement distinct n=1/n>=2 coefficient equations, independent moment repairs, curl-based truncation and controlled summation. This averaging layer does not complete scale recursion.
- [ ] **RM-W11 PULSES / FULL CORRECTED NS / PHYSICAL DYNAMICS.** Implement original mean/oscillatory corrections and averaged quadratic stress cancellation, restricted forcing and independent Cartesian residual checks. Quantify physical core contraction, relative axial elongation, material winding, scale exponents and finite energy.

Mark DONE with source, scoped report/receipt and commit. Full goal remains active. Continue actual prefix and mixed-jet production; no ancestor producer/checker or old-owner/pressure fallback is needed for this layer.
