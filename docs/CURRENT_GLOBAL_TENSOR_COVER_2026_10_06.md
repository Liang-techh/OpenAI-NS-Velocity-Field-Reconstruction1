Historical global physical T/E cover milestone. Continuation: [CURRENT_ORIGINAL_CONE_2026_10_06.md](CURRENT_ORIGINAL_CONE_2026_10_06.md); the original cone map and whole current pulse-end signed two-vector cone are now available. Actual global covariance waves and genuine coefficient recovery remain open.

# Current source-bound global physical tensor cover (2026-10-06)

Implementation and scoped independent global physical T/E cover receipt: commit [a4898929](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/a4898929aa4158348dc228d8dd6de74cacc6e3fb).

## Result

The same checked background T/E graph now has an **independent global physical chart cover**: for every finite physical position, positive time-to-terminal tau and finite positive constant viscosity nu, each fixed admitted construction source has a finite implicit scale and belongs to its original radial atlas, including the Cartesian axis and unbounded exact Gamma exterior.

This is a source-function cover. Parameter uncertainty can still produce several candidate regions; the output retains their complete T/E union. Unique numeric region/seam selection, converged physical u/v/w/p, admissible stress cone/lift, corrected NS, prescribed-domain energy, completed flat remainder and genuine coefficient recursion are not certified by this milestone.

The full graph remains 33 regions / 32 adjacent interfaces / 14 internal supports, with separate same-core analytic axis extension. Primitive velocity/pressure atlas remains 14 / 8. Earlier receipts retain their original narrower scopes and flags; this layer adds its own cover admission without rewriting them.

## Implementation and API

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_global_tensor_cover_operator.py`: independent explicit 33-chart radius atlas, original-source identity comparison, positive derivatives/contiguous endpoints, axial onto-map theorem and source-bound directed candidate argument.
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_global_tensor_cover.py`: actual hypothesis audit, immutable exact common-width witness, checked graph wrapper and source-location replay.
- Matching `.json`, `_check.py` and `_check.json`: source manifest, focused independent admission and fresh physical interface evidence.
- Focused controller: `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentglobalcover`.
- Warm API: `CurrentGlobalTensorCover(field=checked_correlated_tensor_locator)`.

```python
cover = CurrentGlobalTensorCover(field=checked_correlated_tensor_locator)
point = cover.cartesian(x_phys, y_phys, z_phys, time, viscosity=nu, terminal_time=1)
location = cover.locate_log_radius(log_r_phys, z_phys, log_tau, theta, nu)
tensor = cover.tensor_log_radius(log_r_phys, z_phys, log_tau, theta, nu)
source_family_tensor = cover.correlated_tensor(request, z_phys, log_tau, theta, nu)
```

The returned `original_physical_source_result` preserves the complete earlier result and its source flags. The outer wrapper supplies the separately admitted cover scope and common-width witness. Every location is replayed through the source locator; family/source/datum and candidate union cannot be forged to acquire a cover certificate.

## Independent axial onto map

For fixed finite tau>0, nu>0 and the admitted 0<=delta<1, introduce the strictly positive complement C=1-Z^2 on -1<Z<1:

\[
\lambda=\frac{\sqrt\tau}{\sqrt C},\qquad
z_{\rm phys}=\sqrt\nu\,\tau^{(1-\delta)/2}ZC^{-(1-\delta)/2}.
\]

The total derivative, with C'=-2Z, is

\[
\frac{dz_{\rm phys}}{dZ}
=\sqrt\nu\,\tau^{(1-\delta)/2}
  (1-\delta Z^2)C^{(\delta-3)/2}>0,
\quad 1-\delta Z^2=(1-\delta)+\delta C>0.
\]

The two endpoints tend to negative/positive infinity because (1-delta)/2>0. Thus this continuous map is one-to-one and onto all finite z. The positive-real-power identities are checked with C>0, and the same lambda satisfies the original viscosity-normalized implicit equation. Every finite point has |Z|<1 and finite positive lambda. The radial map R=r_phys^2/(2 nu lambda^2) then supplies every R>=0, and periodic angle supplies the x/y plane. Axis uses the analytic Cartesian source and does not require a polar angle.

The numerical log-root program retains a finite directed bracket under its global positive derivative bound. Sign tests change an endpoint only with a sign for the complete parameter box; interval Newton retains every parameter root by the mean value theorem. `step_limit` and `parameter_or_precision_limited` report accuracy limits, not a lost root. A narrow bracket is unnecessary for conservative coverage. The Z enclosure can touch +/-1, while its exact finite physical source remains strictly inside; no finite field value at an infinity-limit endpoint is asserted.

## Independent radial chain and exact common width

The new operator specifies all 33 radius formulas directly, then compares each expression/domain/derivative against the original AST-derived atlas. Every interior radius derivative is strictly positive under actual audited source margins: h>0, both macro widths>0, T>0, loggap-8>0, Md>0, Tw>0, 0<mu<1/4, Lrel-4>0, Ts>0 and wait>0. Core/patch derivatives are 1/native_radius on their positive interiors; the O2 axial derivative is Md exp(Md phase). All other derivatives reduce to positive margins, positive reciprocal mu or 1.

All 32 adjacent endpoints agree exactly. The core log radius tends to negative infinity as rho tends to zero; the final heat log radius tends to positive infinity. Continuity, strict monotonicity and these contiguous ordered endpoints prove every finite logR lies in the atlas. This is a function-level argument, not an inference from sample locations or interval overlap. Shared source traces supply the same T/E at chart joins; uncertainty about which numeric seam was requested does not invalidate union coverage.

The exact shared width is `h_b=epsilon_b=exp(log(cstar)-100 log(K))=cstar K^-100`. A frozen witness binds the original formal recipe, positive cstar, actual physical norm definition, selected C-star/analytic-core/inner-parameter families and ledger hash. Exact transport assignments and the one-upstream bridge/adapter/switch object graph carry that witness to both regions. Equal interval endpoints are only a consistency check. The original finite log(h) datum establishes h>0 even when its numerical cap includes zero.

Both absolute and correlated source paths are AST-bound. Unresolved positive widths retain the entire legal phase domain; every possible region is retained. Core inverse boxes touching rho=0 use the constrained Cartesian same-core source. Heat exterior uses `registry.exterior()` and its identically zero unbounded Gamma T/E, rather than an infinite finite-chart coordinate or radial truncation.

## Focused evidence

- Independent proofs: all 33 original source formulas/domains/strict derivatives, 32 contiguous joins, radial endpoint limits, axial Jacobian/onto limits and original implicit relation, with actual source hypotheses audited.
- Four independent scalar roots with only one numerical solver step, including delta near 1, delta=0 and small time; 16 corner roots in a wide parameter box; an extreme finite log(tau)=-10^1000 bracket. These complement the analytic argument rather than replacing it.
- Four fresh physical/source requests yield nine complete T/E candidate views: **917 component groups / 10,907 retained signed contributions**. They exercise the exact physical axis, an extreme-time off-axis union of five regions, the joint core/bridge boundary, and an unmaterializable exact Gamma radius.
- Absolute outer source uncertainty still retains 25 candidates. Global coverage does not relabel that union as a unique coordinate or resolved field value.
- Eight invalid source hypotheses and four foreign/changed location/request cases are rejected.
- GPT-5.6 Luna / max read-only source/math review accepted the universal cover scope after common-width provenance and source replay were bound.
- Focused controller and checked cover API pass; the fresh receipt is reused only after complete source-hash and manifest matching. The warm graph remains intact.

New gates: `current_original_33_chart_source_radius_ordering_certified`, `current_global_tensor_physical_cover_certified`. Unique coordinate/seam selection and all cone/lift/point/NS/energy/flatness/recursion gates remain false.

## Next construction tasks

- [x] F57C-global-cover1f: independent global physical chart/T/E cover on the fixed admitted source family, every finite point at tau>0/constant nu>0; actual hypotheses, exact axis, unbounded Gamma and conservative candidate unions.
- [ ] **Next F57C-cone1a:** extract the original admissible stress cone and stress-lift equations from the supplied papers. Cite exact version/page/equation; map each physical stress component, normalization and admissibility variable onto this same signed T/E graph. Distinguish exact cone inequalities from existing absolute-value bounds and from relaxed historical comparison cones.
- [ ] **F57C-cone1b:** construct source-correlated regional cone margins and the original lift. Include core/axis, inner exit, macro matching, patch/O2/O3, pulse/gap/end, flatten/angular/steep/waiting/collar and exact zero exterior. Report failing or inconclusive signs and address their source construction; do not set a global cone gate from isolated samples.
- [ ] F57D-recursion1a/b: implement actual n=1 recovery functions and independent moment repair on the common inner interval. Genuine corrections must cancel or absorb the nonflat leading origin remainder: the checked unchanged E_z grows as tau^((-3+delta)/2).
- [ ] F57D-recursion2a/b,3: distinct n>=2 equations, two nontrivial orders, finite-order remainder bounds and final smooth summation. Coordinate-map coverage and repeated spatial scaling do not establish this recursion.
- [ ] F57C-points1/energy1: converged physical u/v/w/p with independent error control, prescribed-domain finite energy and tail integration.
- [ ] F57E/F: both oscillatory pulse families, mean corrections, averaged quadratic stress cancellation, independent corrected Cartesian NS and measured contraction/slenderness/material winding.

Keep the complete long-term objective active. Reuse the admitted graph; move to the original cone/lift construction.
