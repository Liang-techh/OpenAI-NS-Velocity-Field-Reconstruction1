# Rsh physical interface and correlated finite-N primitive shear

Rsh implementation and focused receipt: [d962130e](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/d962130e08ff98ea27035ac921ac203f5378ad10). Joint shear implementation and focused receipt: [c2296cde](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/c2296cde076b6c01559d3e3fc02d59785ac96f1e).

**CONT4f2a and the primitive parts of BOUND1a5a–b are implemented.** The already proved reshape/reference join now has a two-sided physical velocity/absolute-pressure query provider. A new exact correlated finite-N inequality preserves a positive primitive shear reserve through both flat tapers. The independent moment repair and completed signed tensor/direction conditions remain the next cone work; global velocity-interface composition remains open.

## Reshape/reference source to physical trace

CurrentRshPhysicalVelocityInterface in experiments/root_st073/lei_ren_part1_paper_compliant_current_Rsh_physical_velocity_interface.py consumes actual_Rsh_source_join_check and the actual_long_reshape_mixed_C4 / actual_reference_restore_mixed_C4 endpoint packets. It checks the same feedback family, implicit source and analytic pressure datum, rebinds the original source call sites, and retains the two distinct saved source germs.

interface('reshape_reference',Z=(-1,1),log_tau=(-3,-1),theta=None,viscosity=1) returns independent left/right physical traces and216 common spatial4/time1 contribution bounds. The common upper encloses both sides only after the accepted source theorem identifies their functions. Frozen original swirl, sqrt(R/2) radial and Pstar-squared pressure units are consumed after physical differentiation. The mapper does not differentiate those normalization factors again. The finite left long-reshape neighborhood is preserved.

Both exact original radius recipes give log(Rsh)=log(110)+T with the same selected original T. The left packet's T is checked against its defining parent. No gigantic Rsh exponential is materialized. physical_interface(Z,log_tau,theta=None,viscosity=1) uses actual q=log(lambda)=(log_tau-log(1-Z^2))/2 and signed physical z. It accepts finite |Z|<1 and one positive constant viscosity; the source enclosure API permits the Z closure endpoints. Primitive cells still enclose whole Z when the operator Z request is narrower.

Evidence:217 source/radius/unit/linear-map identities;216 common contribution groups;744 independent side-to-common nonzero upper comparisons across two sectors; a fresh log_tau[-12,-11], Z[.1,.4], nu=.8 sector; independent physical lambda and log-radius identities;8 invalid query guards;908 working/index dependencies. The predecessor's135 function-level mixed4 rows are consumed, not relabeled as a new theorem.

Files use prefix lei_ren_part1_paper_compliant_current_Rsh_physical_velocity_interface: .py, _check.py, .json, _check.json, _views.json.gz. Run its producer and checker directly. This provider closes one of the32 adjacent interfaces. The six named patch supports remain available from the previous atlas; a complete single global interface export still needs the remaining source providers and current modified germs.

## Correlated primitive shear through the flat support edges

CurrentCorrelatedFiniteNShear in experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_correlated_finite_N_shear.py binds the actual finite-frequency profile assignments before deriving the inequality. query(offset,N) supports the actual shared O2/O3 logR offset[-2,1], with the current positive mu interval and any finite integer N>=1. It returns chi, chi_y, sigma, original log slope L, actual translated phase N*offset, both finite-N normalized shears, phase-specific absolute error caps and a source-certified joint margin lower envelope.

For these original finite-N profiles, define M_N=a_minus2+b^2/(2+3mu). The exact source-correlated result is:

**M_N >= 2mu*sigma + (mu/4)*chi^2, for every real phase and every finite integer N>=1.**

In the negative O2 buffer sigma=0 and L=-1/2. In O3, L=-1/2-mu*sigma and |L|<=1/2+mu. The actual squared-exponent cutoff has 0<=chi<=1 and mu<=.001. Set S=sin(2pi*phase), C=cos(2pi*phase), d=chi_y/(pi*N*chi), ell=L/(pi*N), and q=exp(-2A/N)/(2+3mu). The exact joint expression after removing2mu*sigma is mu*chi^2 times:

`C^2-S^2+d*S*C+q*(2*S-(d+ell)*C)^2`.

For chi>0, completing the square gives:

`C^2+S^2-ell*S*C-S^2/(4q)+q*(d*C+ell*C-2*S+S/(2q))^2`.

The elementary exponential lower bound and pi>=3 give q>=11999/24036>2/5. Since S^2+C^2=1 and |S*C|<=1/2, the joint coefficient is at least3/8-|ell|/2>=583/2000>1/4. This estimate is uniform in the cutoff logarithmic derivative d, including its divergence near a flat edge. The square retains the coupling between theta and axial shear that independent absolute error/margin ratios discard.

At the exact edges, the numerical API and source endpoint proof use the undivided formula. chi=chi_y=0 yields M_N=0 at O2_buffer9 (offset-2) and M_N=mu at O3 offset1/2. The original zero endpoint is preserved. The lower envelope is strictly positive as a source function wherever chi or sigma is nonzero. A very small cutoff may still have a numerical cap with zero lower endpoint; the directed_scalar_margin_lower_bound explicitly reports that numerical limitation while the source inequality retains its exact positive function.

Phase-specific error caps retain chi, chi_y and L. They use |exp(-x)-1|<=|x|exp(|x|), without microscopic exponential subtraction. Both caps vanish exactly at the two flat source endpoints. They are diagnostic/error inputs; the joint primitive floor is proved with the correlated square, not by requiring separate uniform error ratios.

Evidence:19 exact source/square/exponential/rational/endpoint identities;128 independent evaluations of the original two shear expressions, including either sign and magnitude1e60 of chi_y/chi;8 actual current-source taper queries; exact left-zero/right-mu endpoints and a positive directed floor at the old O2/O3 seam;5 invalid query guards;919 working/index dependencies. Read-only mathematical review: **GPT-5.6 Luna / max**, no worker edits or constructors/tests.

Files use prefix lei_ren_part1_paper_compliant_current_O3_correlated_finite_N_shear: .py, _check.py, .json, _check.json. Run producer and checker directly. No ancestor graph is rebuilt. This is the primitive shear inequality before independent repair. It does not prove alignment, directional pressure/stress conditions, the second signed vector cone, repaired shear, completed tensor admissibility, sufficient common N, energy or coefficient recursion. N>=1 in this inequality does not select N=1 for the final construction.

## Next executable tasks

- [x] **CONT4f2a:** reshape/reference physical spatial4/time1 provider; two source germs, original radius and frozen units, common bounds and actual-lambda/nu queries.
- [x] **BOUND1a5a (primitive):** retain actual chi, chi_y, L and translated phase in local shear/error queries; preserve both exact flat endpoints without dividing by chi.
- [x] **BOUND1a5b (primitive):** complete the correlated square and prove M_N>=2mu*sigma+mu*chi^2/4 uniformly in phase and every finite integer N>=1.
- [ ] **CONT4f2b:** reuse the actual_reference_restore_mixed_C4 packets for reference_restore and restore_buffer. Bind both radius endpoints, current histories and analytic P0 before physical transport; preserve the frozen original units and source branch differences.
- [ ] **CONT4f2c–4 / CONT4f6a4b/6b/7:** add remaining original inner/core adapters, compose existing pulse/outer providers and8 support interfaces, and replace the affected seams with all12 current modified germs. Publish a missing-provider inventory and admit the global interface gate only after complete source coverage.
- [ ] **BOUND1a5c / BOUND2:** retain this primitive correlation in the full signed stress and pressure/directional expressions. Bound alignment and the two-vector cone margins on each taper; do not substitute the primitive inequality for the completed admissible cone.
- [ ] **BOUND1b1:** use the actual five physical density majorants to bound signed cumulative defects, the independent implicit inverse/Jacobian and repair controls versus the same N. Preserve original pressure datum and coefficient uncertainty.
- [ ] **BOUND1b2:** include actual bump derivatives, shrinking partial strips, full histories, terminal implicit closure, radial velocity and pressure recovery. Establish the repaired shear perturbation relative to the correlated reserve.
- [ ] **BOUND2/3 / COMMONN:** combine all repair, radial, diagonal, pressure, alignment and signed-vector inequalities and select one sufficient finite integer N. Primitive N>=1 and the previous repair-only N=10^12 are not final certificates.
- [ ] **ENERGY / REC:** include kinetic cross terms and physical volume/time tails, then actual n=1/n>=2 recovery equations, common core domain, independent per-order moment repairs, divergence-preserving truncation and smooth summation.
- [ ] **WAVE / PHYS / DYNAMICS:** oscillatory/mean correction and averaged stress cancellation, resolved physical u/v/w/p and corrected Cartesian NS residual, measured core contraction, relative axial elongation and material winding.

The controller's last full runtime stage remains currentmodifiedphysicalvelocity. These focused source/physical stages are standalone. Preserve the full persistent reconstruction goal and unrelated working files. Previous source/parameter handoff: [CURRENT_PATCH_ATLAS_AND_ACTUAL_N_BOUNDS_2026_10_06.md](CURRENT_PATCH_ATLAS_AND_ACTUAL_N_BOUNDS_2026_10_06.md).
