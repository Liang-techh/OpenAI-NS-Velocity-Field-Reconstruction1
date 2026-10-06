# Actual implicit physical coordinates and original radius inverses (2026-10-06)

Implementation and scoped actual physical coordinate/radius inverse candidate receipt: commit [95ca039d](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/95ca039db83139964205bcef785728fe8b30079b).

The checked 33-region full tensor graph now has a physical Cartesian entrance and a logarithmic radius entrance. The new layer solves the actual implicit similarity scale, evaluates every original regional radius inverse as a directed candidate enclosure, and connects those candidates to the complete original T/E views. It keeps the actual time parameter and all signed source sectors.

This completes F57C-global-cover1e-1 and the candidate-formula part of 1e-2. It does not complete exact physical seam selection, correlated outer-radius location, global physical cover, converged u/v/w/p, admissible stress/lift, coefficient recursion, completed flatness, energy or corrected NS. No percentage of the full reconstruction is inferred from the coordinate milestone.

## Source and API

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator.py`: original radius AST before numerical caps, 32 symbolic endpoint identities, exact implicit physical identities, stable directed root map and constrained Cartesian core source.
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_physical_tensor_locator.py`: one admitted registry, actual source parameter bindings/positive slope guards, 33 inverse candidates and full T/E routing.
- Matching `.json`, `_check.py` and `_check.json`: source manifest and independent scoped admission.
- Focused controller: `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentphysicallocator`.
- Warm construction: `CurrentPhysicalTensorLocator(registry=checked_current_tensor_registry)`.

```python
field.cartesian(x_phys, y_phys, z_phys, t, viscosity=nu)
# For a locator without full tensor evaluation:
location = field.cartesian(x_phys, y_phys, z_phys, t, viscosity=nu, tensors=False)

# Avoid materializing astronomical/small absolute radii:
location = field.locate_log_radius(log_r_phys, z_phys, log_tau, theta, nu)
view = field.tensor_from_location(location)
```

`terminal_time=1` by default; time-to-terminal is tau=terminal_time-t and must be strictly positive. Inputs and constant nu must be finite, nu positive. `log_r_phys=-inf` represents the exact axis. Input boxes crossing the axis retain all possible coordinates and regions. The output is a union of complete directed source-function tensor enclosures, not a selected midpoint field or a converged physical velocity value.

## Actual scale and time semantics

With q=log(lambda), the original constant-viscosity map is

\[
\lambda^2-\frac{z_{\rm phys}^2}{\nu}\lambda^{2\delta}=\tau,
\qquad Z=\frac{z_{\rm phys}}{\sqrt\nu\lambda^{1-\delta}},
\qquad R=\frac{r_{\rm phys}^2}{2\nu\lambda^2}.
\]

The production root is enclosed through

\[
F(q)=2q-\operatorname{logaddexp}\left(\log\tau,
2\log|z_{\rm phys}|-\log\nu+2\delta q\right),
\qquad 2(1-\delta)\le F'(q)\le2.
\]

The lower bracket is directed log(tau)/2. The upper bracket makes each positive term no more than half exp(2q). Bisection and interval Newton retain the complete parameter box; trial midpoints never replace a source parameter. z=0 and sign-crossing boxes are handled separately. Very negative exponentials receive conservative nonnegative enclosures, while their exact logs remain attached to the source. Astronomical scale factors are not materialized.

The true finite-point relation log(1-Z^2)=log(tau)-2q is retained. Its complement also narrows the directed Z enclosure. That enclosure may touch Z=+-1 when precision cannot resolve the strictly positive true complement; the true finite physical point still satisfies |Z|<1.

Native regional APIs retain actual `log_tau`. Their `gamma*log_tau/2` was a conservative sqrt(tau) factor for nonpositive gamma. It is not the off-axis actual scale. Each retained signed row now has separate `actual_log_lambda`, `actual_lambda_log_prefactor=gamma*q` and `actual_lambda_log_absolute_upper`, alongside its unchanged original native row and conservative bound.

## Original radius sources and ambiguity

All 33 inverse recipes come from the existing original radius program. Only its two numerical logR clamps and two rounded singleton shortcuts are removed in the symbolic replay; every exact source expression and original branch remains. Native adapters retain the original core rho, pulse xi, reciprocal gap phase and distinct O2 axial/buffer coordinates. Symbolic simplification removes common large bases before deriving affine widths. The 32 source endpoint identities are exact algebraic equalities, not interval overlap tests.

The width h is the original positive source exp(logh). Its numerical `[0,cap]` bound participates in candidate range/inverse calculations only as an enclosure. Candidate metadata explicitly reports this use and retains the exact logh source. Unresolved positive widths use their full original legal phase domain. Runtime guards enforce positive Md/T/Tw/Ts/wait, 0<mu<1/4, loggap>8, Lrel>4 and the positive bridge/switch macro widths.

At the core, the original full Cartesian source is replayed with an independently enclosed rho under the exact joint polar constraint X=sqrt(2rho)cos(theta), Y=sqrt(2rho)sin(theta). This avoids the artificial outside-disk radius from independently squaring X/Y interval boxes at rho=4. All original physical derivative operators and six remainder sectors remain unchanged. The exact Gamma exterior uses its admitted unbounded identically zero T/E source, so a wide tail offset cannot trigger astronomical exponentiation.

The current source data exposes a real limitation: at ordinary moderate physical radii the locator selects one macro region; at absolute outer log radii, source parameter uncertainty can overlap **25 candidate regions**. The returned union is conservative and explicitly ambiguous. It does not resolve short local offsets relative to an astronomically large common radius. Increasing ordinary precision or selecting a midpoint would not complete that construction.

## Focused evidence

- Seven independent direct-exponential root fixtures, including zero/negative/tiny/large z, q below zero and delta zero; two parameter boxes with 16 corner roots/Z each; one astronomical negative-exp enclosure case.
- Original native coordinate enclosed by the inverse for all 33 regions; all 32 original symbolic radius endpoint identities replayed.
- Fresh physical axis, ordinary point and joint core-boundary routing: six complete T/E views, 704 component groups and 15,608 retained signed contributions. The boundary union includes the core and original bridge candidates.
- Actual lambda factors are checked separately from unchanged native source rows; actual off-axis q differs from log(tau)/2.
- Invalid terminal/after-terminal time, nonpositive nu, nonfinite x and forged candidate lists rejected.
- GPT-5.6 Luna / max read-only mathematical/source review; width-cap provenance and source positivity feedback incorporated.
- Focused controller reuses the real fresh receipt only after full source hashes and manifest match; no historical producer rebuild or second numeric scan.
- True gates: `current_actual_implicit_physical_coordinates_directed`, `current_33_region_original_radius_inverse_candidate_locator_available`.
- Exact physical seam selection/global locator/global physical cover and all cone/lift/NS/point/energy/flatness/temporal recursion gates remain false.

## Next construction tasks

- [x] F57C-global-cover1e-1: recover the actual implicit lambda relation and constant-viscosity map; direct finite-coordinate inversion with real log(tau), true |Z|<1 relation and locator width/status.
- [x] F57C-global-cover1e-2a: all 33 original inverse candidate formulas, positive source guards, axis path and full T/E candidate union with separate actual lambda factors.
- [ ] **F57C-global-cover1e-2b:** represent physical log radius with an exact shared original anchor and signed local offsets, retaining correlations in logRref/logRp/logRv and the microscopic exp(logh) source. Do not collapse these to one absolute MPF interval. Bind the anchor to the same admitted family/source/datum and original radius AST.
- [ ] **F57C-global-cover1e-2c:** use the correlated representation to distinguish axial_restore/buffer/patch/O2/O3/pulse/flatten/angular/steep/waiting/collar domains and exact microscopic phases. Return sound inverse coordinate errors and explicit unresolved cases; no midpoint/cap source or arbitrary single candidate.
- [ ] **F57C-global-cover1e-2d:** bind named adjacent and support boundary source recipes to physical requests. Use the already admitted 32 adjacent / 14 support traces for exact common physical rows. Axis remains a separate analytic extension, not a counted support trace.
- [ ] **F57C-global-cover1f:** independent finite-positive-time physical cover using actual implicit source, all exact radius endpoint identities and audited positive-width hypotheses; include the axis and full unbounded Gamma. Publish its own source-bound scoped receipt. Candidate inverse formulas alone must not set this gate.
- [ ] **F57C-cone1a/b:** recover the original cone/lift conditions on this same global signed T/E; compute correlated regional margins and retain failures as construction targets.
- [ ] **F57D-recursion1a/b, 2a/b, 3:** actual n=1 and distinct n>=2 operators/functions/moment repair, at least two nontrivial orders, leading origin E cancellation/absorption, finite-order remainder estimates and smooth sum. The existing strictly positive origin E_z grows as tau^((-3+delta)/2); the unchanged leading E is not flat.
- [ ] **F57C-points1/energy1; F57E/F:** convergent physical u/v/w/p and prescribed-domain energy, oscillatory families/mean corrections, independent corrected Cartesian NS and measured contraction/slenderness/material winding.

The full long-term reconstruction goal remains active. Reuse the admitted graph and build the correlated locator next.
