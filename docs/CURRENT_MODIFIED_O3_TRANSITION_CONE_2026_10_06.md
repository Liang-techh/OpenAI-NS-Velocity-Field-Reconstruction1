# Modified closed O3 signed stress cone, all finite integer N >= 22

Latest successor: [CURRENT_MODIFIED_O2_TAPER_CONE_2026_10_06.md](CURRENT_MODIFIED_O2_TAPER_CONE_2026_10_06.md), implementation [f60f89f8](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/f60f89f8747166009d6e8daac3479e5842c88e19). Modified O2 open taper is now complete for N>=22; the left support edge remains degenerate and quiet/common-N/global gates remain open.

Checked implementation: [bdbe19b1](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/bdbe19b13b3024c1a42bdd38514cfa258d06763f). **The actual modified source now satisfies the signed two-vector stress cone on the entire closed O3 offset interval [0,1], all Z in [-1,1], all oscillation phases and every finite integer N>=22.** This includes the old zero-shear seam at offset=0 and the right flat taper. The complete original pressure, energy, incoming moments and all stress cross terms are retained.

This is a scoped O3 source theorem. The left O2 taper, the quiet repair cone and a sufficient common frequency remain open. N=22 is sufficient for this region only; it does not satisfy the separate quiet repair threshold. The original unmodified registry counts remain unchanged. No finite-energy, actual coefficient-recursion or full corrected NS claim is made.

## Callable implementation

CurrentModifiedO3TransitionCone in experiments/root_st073/lei_ren_part1_paper_compliant_current_O3_modified_transition_cone.py consumes checked completed tensor error majorants, correlated finite-N shear and the checked whole original O3 direction bounds. query(offset,N) accepts finite interval subsets of [0,1] and integer N>=22. It returns signed g/h interval diagnostics, the actual phase N*offset and the source theorem. The proof is uniform in phase and does not depend on those diagnostic phase boxes. No ancestor constructor, fixed-N modified source value or saved control midpoint is used.

Scoped completed gate: current_modified_closed_O3_signed_two_vector_cone_certified. Remaining O2/quiet/common-N/global/energy/recursion/remainder/full-NS gates are false. Producer/checker .py and hash-bound .json/_check.json use the common prefix lei_ren_part1_paper_compliant_current_O3_modified_transition_cone.

## Signed algebra and units

Let K=Pstar*sqrt(R/2)>0, theta=Ttheta/K and zeta=sqrt(mu)*Tz/K. Actual shears are a=2+mu*g and bs=+sqrt(mu)*h=+2*Uz_y/Utheta. The common physical Pstar cancels with its original positive sign.

- H=a*theta-h*zeta=a*D/K, where D=Ttheta-(bs/a)*Tz.
- Kc=a*zeta+mu*h*theta=sqrt(mu)*a*J/K, where J=Tz+(bs/a)*Ttheta.
- mhat=2*g+h^2+mu*g^2=a*(vs-2)/mu.
- W=(2*a-mu*h^2)*theta^2-2*a*h*theta*zeta-a*g*zeta^2.
- a^3*(Q/K^2)=2*a*H^2-mhat*Kc^2=(a^2+mu*h^2)*W, where Q=2*D^2-(vs-2)*J^2.

These source identities avoid division by a cutoff or axial coordinate. The weighted axial component keeps microscopic shear factors rather than discarding their cancellations in a large unweighted axial cap.

## Whole-region inequalities

The exact source and genuine cutoff cap |chi_y|<=128 give |g|<3, |h|<4 and 1.997<=a<=2.003 for all N>=22. The accepted correlated primitive lower transfers to the actual vs-2 because a<=2+3*mu. For offset<=1/4, chi=1; for offset>=1/4, sigma>=sigma(1/4). Thus vs-2>=mu*min(1/4,2*sigma(1/4))>0 on the closed O3 interval. The old seam offset=0 is strictly positive after modulation. The left O2 edge offset=-2 has exact zero primitive excess and is not included.

The original whole O3 baseline bound gives a signed positive theta lower and a bound on 2*mu*(Tz0/Ttheta0)^2, with its canonical U factor and C(Z) normalization preserved. All 14 native stress error pieces (7 theta and 7 axial) from the complete source lie below exp(-1010) relative to that baseline lower at N=22. Summing each component gives at most exp(-1000). Formal logarithms are retained; astronomically small profile amplitudes and giant radii are never materialized.

Leading stress order zero uses local velocity through y1 and cumulative histories at y0 only. The source replay verifies this dependency and excludes higher y derivatives that can grow with N. Those local/history absolute caps are positive increasing functions of q=1/N, and the stress error bound is a sum of positive monomials. This proves the N=22 ledger is valid for every N>=22; checks at N=37 and 10^12 are supplementary comparisons.

The modified weighted axial ratio is below exp(-390). Resulting uniform normalized reserves are H/theta>1.9, W/theta^2>3.9, D/Ttheta>0.9 and Q/Ttheta^2>1.9. These are estimates for the complete modified signed source, including the retained pressure and moments.

## Focused evidence

The checker passes 37 exact source/weighted-cone/all-N identities, 25 strict whole-source inequalities, 64 independent original signed-cone fixtures with both admitted and rejected cases, 14 independent scale identities, 28 larger-N native error comparisons, 6 current queries and 7 invalid guards. Working and staged hashes agree for 927 bound files. Read-only source review: GPT-5.6 Luna / max.

## Next tasks for GitHub agents

- [x] **BOUND2b1-O3:** derive actual signed a/bs/vs/D/J/Q with physical normalization, both signs of bs and correlated primitive reserve.
- [x] **BOUND2b2-O3:** transfer the whole original positive stress/direction to the complete modified source; prove all-N error monotonicity and strict D/Q margins.
- [x] **BOUND2b3-right:** include the closed old seam, right taper and right flat edge without dividing by chi.
- [ ] **BOUND2b-left1:** derive whole original O2 baseline stress bounds on offset[-2,0] from its actual pure-buffer source. Preserve canonical U, pressure, energy and incoming moment factors. Do not apply the original positive-offset O3 certificate outside its domain.
- [ ] **BOUND2b-left2:** transfer the complete O2 native stress error ledger to the same weighted cone. Prove signed direction/alignment on (-2,0], all phases and all sufficiently large N using actual chi/chi_y. Keep arbitrarily small positive cutoff factors correlated until cancellation is complete.
- [ ] **BOUND2b-left3:** prove the exact chi=0 flat-edge identities at offset=-2 and reconcile that degenerate source with the accepted adjacent buffer cone notion. Do not label vs-2=0 as strict or silently remove the endpoint from a global domain.
- [ ] **BOUND2b-buffer:** audit the unchanged earlier O2 buffer offset[-11,-2], its a=2/bs=0 geometry and exact transported stress; establish its correct admission/closure rule and joins.
- [ ] **BOUND2c1:** derive quiet bs=+2*Uz_y/Utheta directly from original source assignments. The earlier uniform repair majorant's negative b label was only used for absolute caps; it is not a signed quiet-cone convention. Bind the corrected signed theorem without rewriting ancestor receipts.
- [ ] **BOUND2c2:** specialize the completed error sectors on each disjoint quiet bump and all gaps using the unique finite-N implicit controls. Retain partial/complement histories, absolute cumulative pressure, original radial cross terms and C versus C^2 factors. Obtain local signed baseline theta/axial bounds rather than a single coarse whole-strip cap.
- [ ] **BOUND2c3:** prove quiet D>0/Q>0 with its a-2 reserve, every entry/internal/terminal join and all N above an explicit threshold. Exact terminal source equality alone does not prove the preceding strip.
- [ ] **BOUND3 / COMMONN:** combine repair, history, pressure/radial, diagonal, primitive, direction and alignment constraints into one sufficient finite N. Current O3 threshold22, repair27,303,666 and quiet primitive68,533,403 are distinct partial bounds; publish the complete inequality ledger.
- [ ] **CONT4f2b:** implement accepted reference_restore/restore_buffer physical interface providers, remaining adjacent/internal providers and global composition with each mapper's documented lambda convention. Do not relabel obsolete full33 views.
- [ ] **ENERGY:** bound physical kinetic cross terms using the actual Jacobian and radial/time tails. Moment closure does not erase positive kinetic density.
- [ ] **REC:** actual n=1 and n>=2 recovery equations, common inner domain, independent moment repairs, divergence-preserving truncation, finite-order remainder and smooth summation. Leading source cone certification is not coefficient recursion.
- [ ] **WAVE / PHYS / DYNAMICS:** two oscillatory/mean families, averaged quadratic stress cancellation, resolved u/v/w/p, full corrected Cartesian residual and measured contraction/aspect ratios/material winding.

Mark only the completed scoped task, publish checked source hashes and remaining gates, and update the top-level handoff. Preserve unrelated dirty files and historical sections. Previous report: [CURRENT_COMPLETED_TENSOR_ERROR_BOUNDS_2026_10_06.md](CURRENT_COMPLETED_TENSOR_ERROR_BOUNDS_2026_10_06.md). Full goal remains active.
