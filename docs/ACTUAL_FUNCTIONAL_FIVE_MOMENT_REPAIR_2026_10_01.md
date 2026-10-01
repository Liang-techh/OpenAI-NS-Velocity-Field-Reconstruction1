# F30: actual functional five-moment repair

The F27–F29 family now admits the Section 10 five-bump repair as functions of the entire real axial variable Z in [-1,1]. The source core, selected logarithmic radius, preheat pressure P0, and fixed cutoff/bump data remain the same. This is an implicit smooth coefficient family with callable directed enclosures, not a finite pointwise reconstruction of the whole physical field.

## Source defects and admission

`experiments/root_st073/lei_ren_part1_paper_shared_five_defect_admission.py` derives all five centered defects from the actual inherited cumulative primitives. It includes the R <= 110 contributions, long reshaping, axial restoration, and the reference contributions down to the axis. The apparently singular reference pressure integral is integrable and is included exactly. The scalar restoration integrals are enclosed by 2048 directed closed cells.

The signed defects have the form d1 = k1 E1 + tau1, d2 = k2 E1 + tau2, d3 = tau3, d4 = k4 Am^-2 E1^2 + tau4, and d5 = tau5, where E1 = v1 - 4Z and Am = exp(-.6) Pstar/(1+Z^2). The source identity and whole-axis C1 bound are retained, including all positive nonzero tail caps. Source functions and tails are not replaced by interval midpoints or zeros.

The complete summed C1 defect upper bound is 4.9017390830378e-23, or approximately 0.0350558 percent of e_star. This is the five-moment admission quantity, not the Navier–Stokes momentum residual. The complete sufficient Section 10.19 test passes without changing Cstar or P0.

This norm bound is for the correlated true source functions. Independent value/derivative boxes enlarge that source class: their arbitrary summed box norm can be larger (j+2rho versus j+rho, or twice a tail cap). The analytic Banach proof uses the true summed source norm; the separate interval inverse is certified on the enlarged boxes. These two bounds are not interchanged. The callable patch also fails closed when a denominator cannot be proved positive and checks the selected-family Rm >= 16 provenance.

## Functional correction and coupled fields

`lei_ren_part1_paper_shared_five_moment_repair.py` uses the original directed bump integrals and fixed CA = 344, CQ = 1081, CS = 57749553. The actual physical delta is recovered from the pressure datum; the analytic tube's delta upper bound is not substituted as a parameter.

The scalar C1 Banach bound gives ||h||C1 <= 3.37239648913e-20. A separate anisotropic source gate covers the angular coefficient box, including the complete C1 Am^-2 E1^2 source. The directed inverse proves strict self-inclusion, contraction, and derivative recovery for the same actual source family. The original five moments are recovered from the invertible centered system.

The API `SharedFiveMomentRepair.coefficients(Z)` encloses the unique five coefficient functions and their axial derivatives. `evaluate_patch(x,Z)` computes swirl, axial velocity, partial cumulative moments, pressure P0+Mp, radial velocity, and normalized stress quantities from these same data, for x = R/Rm in [1,e]. Correction support is strictly inside (1,2), namely [49/40,71/40]. Partial bump integrals use directed closed cells; completed supports use the original full integral enclosures.

The corrected inner relaxed-cone bound remains positive through Rh: |bw| <= 8.9742822813401e-6, kappa < 1, and the recorded conservative cone margin exceeds 3.7999. The original core and inner admissible collar are preserved.

## Independent check and limits

`lei_ren_part1_paper_shared_five_moment_repair_check.py` derives five correction integrands directly from the physical moment definitions, differentiates five cumulative primitives, and independently recovers Q, N and Ur. All 13 symbolic identities pass. Five callable interval examples pass their field/cone diagnostics. These finite examples are diagnostics; the whole-axis closure claim rests on the functional contraction and source bounds.

Read the producer and independent checker receipts together. The producer's pending-independent-check flag describes its own generation boundary; the companion receipt provides the separate completed check, without circular hashes.

Remaining work: build the outer candidate and its own angular/axial corrections at the selected Rref, recompute the exact heat profile and its correction coefficients, close pressure and radial velocity across the outer/heat match, then construct the global admissible stress and flat remainder. Actual n-dependent temporal coefficient recursion, oscillatory corrections, finite-energy integration over the complete domain, and full Cartesian residual validation remain unfinished. Physical Cstar/Rm are not materialized, exact coefficient point functions are not evaluated, and Ur_Z is not available from this C1 patch API.
