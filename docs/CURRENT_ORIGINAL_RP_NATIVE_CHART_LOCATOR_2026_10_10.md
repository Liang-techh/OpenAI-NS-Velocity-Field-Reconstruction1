# Update: whole-box source callbacks installed (2026-10-10)

[Source e9060cd4](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/e9060cd43351da8f95053a4633b55c5143fc216a) completes the bounded fifteen-chart native/Z interval callback layer. Read the [new source API and next tasks](CURRENT_ORIGINAL_RP_NATIVE_BOX_SOURCE_2026_10_10.md). Correlated original-scale physical inversion, successful physical source admission/error delivery, seam calls and u/v/w/p materialization remain open.

---

# Current original Rp native chart locator (2026-10-10)

Full reconstruction **ACTIVE / INCOMPLETE**. [Source ca93c2e7](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/ca93c2e7ea3a7e4cb7ed67b34664330d6f4dbb0a) installs the exact affine inverse and directed radial locator for all fifteen current Rp-to-heat charts. It consumes unchanged live physical inverse results or supported current exact logR function nodes, returns every possibly overlapping chart and its native-coordinate interval, and identifies radii definitely before Rp. A source-bound directed logRp evaluator is now available. Native interval source callbacks, seam source calls and numerical physical u/v/w/p remain open.

## Current callable interface

```python
owner = CurrentOriginalRpNativeChartLocator(accepted_physical_inverse)
physical = owner.before.cartesian('0.8', '-0.3', '0.2', '0.95')
located = owner.locate_inverse(physical)
located['classification']             # before_current_Rp for this input
located['possible_charts']
located['chart_rows']['heat_exterior']['directed_native_coordinate']

logR = owner.source_log_radius('heat_exterior', '9/2')
native = owner.locate_function(logR)
native['contained_charts']            # ['heat_exterior']
native['chart_rows']['heat_exterior']['exact_native_coordinate_function']
```

Implementation and receipt quartet: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_chart_locator.py`, `.json.gz`, `_check.py`, `_check.json`. Reuse accepted current runtime objects. No source selection/repair/future is replayed, and no old tensor registry or 33-chart owner is constructed.

The public function path accepts only exact source logR handles registered by this locator through source_log_radius or source_endpoint_log_radius. It does not accept a supplied interval or numeric evaluator. The inverse path requires a live solve witness from this exact accepted inverse owner, reassembles the view and compares the complete packed result. A mutated view or imported numerical row is rejected.

## Current source bounds and cancellation

The accepted current absolute source identity is

```
logRp = 10*logC + 11*logP + log(110) + 1 - 60*log(mu)
logP = exp(40)+11
mu = (1/1000)*exp(-4logP)
```

The current heat-closure source proof already identifies this formula with the full original Rc/quiet-to-Rp graph. The new evaluator binds logC to the existing repair heat owner's exact selected dyadic singleton; it validates its source leaf/family/receipt. Mu, Lrel, Ts and integral-defined waiting W have directed bounds from the same accepted current pulse/postpulse owners. These bounds enclose exact functions and never redefine them.

An exact rational polynomial normalizer expands graph sums/products/negations, canonicalizes commuting monomials and removes zero coefficients before directed arithmetic. The microscopic hbB*sc origin cancels exactly without needing its numeric value. Common logRp, pulse length, Lrel/Ts/W terms likewise cancel before native differences are evaluated. Quotients and analytic/integral functions remain exact graph atoms; no waiting-box endpoint or arbitrary pulse coordinate is selected.

Absolute **log**-radius bounds can be held as MPF exponent/mantissa values. The physical R=exp(logR) is not materialized. Only the original exp(40) and mu parameter nodes are whitelisted for exponential evaluation. Unregistered constants or moderate exponentials, huge absolute scales and unknown source variables are rejected. Native domains and the exponential whitelist are rebound in assert_graph.

## Fifteen actual affine maps

Let B=logRp, P=13/mu, L=Lrel, T=Ts and W=waiting. Each inverse retains the exact function `(logR-base)/J`, where base is the current original radius map at native coordinate zero.

| Chart | Native domain | base | J=dlogR/dv |
|---|---|---|---|
| pulse_entrance | [0, (1/50)/mu] | B | 1 |
| pulse_main | [1/50, 10] | B | 1/mu |
| pulse_exit | [10, 11] | B | 1/mu |
| pulse_gap | [11, 12] | B | 1/mu |
| pulse_gap_end | [-1/mu, -4] | B+P | 1 |
| pulse_end | [-4, 0] | B+P | 1 |
| flatten | [0, 100] | B+P | 1 |
| outer_power | [0, 1] | B+P+100 | L-4 |
| outer_angular | [-4, 0] | B+P+100+L | 1 |
| steep_entry | [0, 1] | B+P+100+L | 1 |
| steep_power | [0, 1] | B+P+101+L | T |
| steep_exit | [0, 1] | B+P+101+L+T | 1 |
| waiting | [0, 1] | B+P+102+L+T | W |
| heat_collar | [0, 3] | B+P+102+L+T+W | 1 |
| heat_exterior | [3, infinity) | B+P+102+L+T+W | 1 |

These expressions and Jacobians come directly from the current `_map` definitions. The heat exterior is explicitly classified as an unbounded native tail.

`possible_charts` retains every chart whose domain may intersect the directed native interval. `contained_charts` requires the entire native enclosure to lie inside the conservative inner bounds of that chart domain. Each result retains both the unconstrained exact native function/enclosure and a **conditional** chart intersection. The clipped interval applies only if the point lies in that chart; it is not the selected native coordinate. A boundary overlap is not rounded into one provider. The locator does not call a source seam or claim source point accuracy. The physical inverse's numeric |Z|<1 resolution flag is preserved separately.

## Focused evidence

- Fifteen actual current interior source-function anchors: exact affine inverse identities and directed native containment pass. No velocity packet is used as a coordinate value.
- Fourteen actual current inherited radial boundaries: exact same-radius identities pass and both adjacent chart candidates are retained. A private interval cover fixture also retains flatten/power overlap; it is labelled a cover fixture, not a physical source value.
- Nine unchanged live Cartesian/log input inverses are routed. These particular test points are all definitely before the current Rp radius, so no outer source provider is selected for them. This does not supply the missing inner/global field.
- Independent directed evaluation of the original logRp formula overlaps the polynomial graph evaluation. The selected logC dyadic and exact original N definition are retained.
- Twelve inadmissible-source rejections: numeric radius boxes, unregistered constants/moderate exponentials, unknown graph variable, absolute-radius exponential, foreign graph, foreign-owner/mutated inverse results, substituted logC bound, changed native domain, expanded exponential whitelist and an unadmitted exponential oracle.
- Accepted constructor, fresh exterior native function and fresh physical inverse calls pass. Source dependency hashes are audited against the Git index before upload. The sole read-only reviewer remains GPT-5.6 Luna / max.

## Completed bounded tasks

- [x] **CURRENT-RP-NATIVE-AFFINE-LOCATION** Current fifteen affine inverse functions, directed domains, all possible/contained charts, before-Rp classification and unbounded heat-tail label.
- [x] **CURRENT-RP-SOURCE-BOUND-LOGRP** Directed current logRp evaluation using accepted exact source identity and live selected logC datum.
- [x] **CURRENT-RP-NATIVE-COMMON-FACTOR-CANCELLATION** Exact rational graph cancellation before numeric local differences, retaining original source scales and Jacobians.
- [x] **CURRENT-RP-RADIAL-BOUNDARY-CANDIDATES** Retain both neighboring charts at all fourteen inherited radial seams. Source seam calls remain separate.

## Next implementation tasks and acceptance

Priority is interval source delivery, followed by grouped physical values. The predecessor contains the full stress/recursion/oscillatory task backlog; those tasks remain open.

- [ ] **CURRENT-RP-PHYSICAL-INVERSE-CONTRACT** Finish source admission/interval consumption. The physical inverse and radial locator are now installed, but a possible radial chart plus exact |Z|<1 is not a sufficiently accurate source point. Require interval-native source delivery or meaningful refinement, and propagate coordinate errors.
- [ ] **CURRENT-RP-NATIVE-CHART-LOCATION** Finish the broader caller contract: the affine radial subtask is complete, while accepted seam provider calls and interval-native evaluation remain open. Do not repeat the fifteen affine inverse implementation.
- [ ] **CURRENT-RP-CORRELATED-LOG-PHYSICAL-INPUT** Extend current inverse inputs to source-bound exact graph log_r/log_tau/log_abs_z where needed. Current scalar rational/binary inputs cannot conveniently encode the original correlated astronomical logC/time scales. Preserve same-graph provenance and cancellations; no forged interval evaluator or practical substitute radius. Exercise an actual outer chart with a correlated near-critical time and known source coordinates.
- [ ] **CURRENT-RP-INTERVAL-SOURCE-INPUT** Partial: actual fifteen-chart native/Z whole-box callbacks with exact function scales and directed mixed rows are installed in the new native box source module. Independent source cells pass; successful correlated original-scale physical inverse consumption and requested physical error remain open. Do not repeat the completed wrapper adapters.
- [ ] **CURRENT-RP-NATIVE-BOUNDARY-CALLERS** Invoke the accepted fourteen same-function seam callers with explicit left/right unit conversions when an exact boundary is established. For a crossing interval enclose both sides or refine/partition into admitted exact rational cells. Preserve chart intersections as conditional.
- [ ] **CURRENT-RP-NATIVE-LOCATOR-REFINEMENT** Report when native or axial coordinate widths prevent source accuracy. Refine the true inverse/source parameter enclosures or use correlated exact expressions; preserve every possible chart until exclusion is established. Do not declare point admission from a finite-width location alone.
- [ ] **CURRENT-RP-HEAT-TAIL-INTERVAL-DELIVERY** Use the explicit unbounded heat-exterior classification to bind actual Gamma/tail evaluation and error for native intervals. Do not treat the collar endpoint or a truncated practical radius as the exterior.
- [ ] **CURRENT-RP-RAW-POINT-ERROR-UNITS** Give real directed phase/quadrature/coefficient/Gamma errors in each source amplitude unit, including coordinate uncertainties. Convert these through each physical scale group and state achieved absolute/relative accuracy.
- [ ] **CURRENT-RP-EXACT-POINT-GRAPH-EVALUATOR** Evaluate admitted physical inverse/atan2, sin/cos and signed integer-power nodes with source binding and correlated log arithmetic. Native log-geometry evaluation does not implement the full physical-scale oracle.
- [ ] **CURRENT-RP-PHYSICAL-ROW-MATERIALIZER** Produce directed finite u/v/w/p and spatial/time rows from admitted source packets, exact grouped scales and all errors. Retain signs, pressure radial/baseline unit split, independent P0 and original N.
- [ ] **CURRENT-RP-PHYSICAL-POINT-API** Combine current inverse, radial/axial admission, interval source/error callback and materializer. Demonstrate actual pulse, angular/flatten, waiting/collar/exterior points and inherited seams, not only coordinate anchors or synthetic fixtures.
- [ ] **CURRENT-RP-PHYSICAL-POINT-INDEPENDENT-CHECK** Independently check completed physical queries and fixed-x/spatial derivatives with real source accuracy and a correlated original-scale input.
- [ ] **CURRENT-GLOBAL-MIXED-DERIVATIVES** Add actual current core, transition, Rh/O2/O3, compact/quiet and leading-to-Rp rows and seams to the same physical graph. Ordinary physical inputs which route before Rp need these providers.
- [ ] **CURRENT-PHYSICAL-AXIS-LIMIT** Construct actual current core parity/limits and Cartesian derivatives/pressure/divergence at r=0.
- [ ] **CURRENT-CARTESIAN-VELOCITY** Deliver the full callable [u(x,y,z,t),v(x,y,z,t),w(x,y,z,t)] after point/error/global/axis admission. Publish actual supported domain and error.
- [ ] **CURRENT-RP-POSTPULSE-QUANTITATIVE-SEAMS** / **CURRENT-RP-PULSE-UNIFORM-BOUNDS**: establish genuine full-domain mixed derivative/support bounds beyond finite diagnostic anchors.
- [ ] **CURRENT-RP-STRESS-DERIVATIVE-ORDER-LEDGER** / **CURRENT-STRESS-CONE-AND-REMAINDER**: construct original signed stress and independent flat remainder with needed derivatives, cone margins and max/physical-volume L2 bounds.
- [ ] **CURRENT-PHYSICAL-ENERGY**: integrate actual physical-volume kinetic energy and full tail with the lambda-dependent Jacobian.
- [ ] **CURRENT-TEMPORAL-RECURSION-N1** / **CURRENT-TEMPORAL-RECURSION-HIGHER**: actual recovery equations, independent five-moment repairs, curl-preserving truncation and smooth summation.
- [ ] **CURRENT-OSCILLATORY-FAMILIES**: both original families/mean corrections, averaged quadratic stress cancellation and independent remainder.
- [ ] **CURRENT-DYNAMICS-DIAGNOSTICS** / **CURRENT-FULL-CARTESIAN-RESIDUAL**: actual widths/aspect ratio/growth/material winding/recurrence/energy, then full corrected forced NS max and physical-volume L2 residual.

Predecessor: [directed physical inverse and complete detailed backlog](CURRENT_ORIGINAL_RP_PHYSICAL_INVERSE_2026_10_10.md). Start with interval-native source input and error delivery. Keep full reconstruction ACTIVE.
