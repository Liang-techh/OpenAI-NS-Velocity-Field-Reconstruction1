# Current correlated source radius tensor locator (2026-10-06)

## Result and scope

The admitted full background T/E graph now accepts an exact original radius anchor plus a signed local offset. Common logRref/logRp/logRv terms are cancelled symbolically before interval inversion. This distinguishes source locations which were mixed into a union of 25 regions when represented by one absolute radius interval. All 33 independently specified interior source targets select their expected region; all 32 adjacent and 14 internal support radii select their original complete tensor trace route.

These requests describe **source-dependent physical point families**. A request denotes the physical point at the construction's exact radius, for each admissible source parameter. This does not resolve one fixed absolute physical coordinate or produce converged physical u/v/w/p. The earlier absolute-coordinate API still returns its conservative ambiguous union. Arbitrary physical boundary selection, independent global physical cover, original admissible cone/lift, energy, corrected NS, completed flatness and genuine temporal recursion remain open.

The graph and counts remain 33 regions / 32 adjacent interfaces / 14 internal supports. Primitive velocity/pressure atlas remains 14 / 8. Axis is a separate same-core analytic extension, not an added region or support count.

## Implementation

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_correlated_radius_operator.py`: original exact source equalities, positive microscopic width source, canonical rational logarithms, original radius derivatives/inverses and exact 46 boundary recipes.
- `experiments/root_st073/lei_ren_part1_paper_compliant_current_correlated_tensor_locator.py`: source-family requests, all-region inverse candidate scan, complete T/E views and original boundary trace dispatch.
- Matching `.json`, `_check.py` and `_check.json`: manifest and focused source-bound admission.
- Controller stage: `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentcorrelatedlocator`.
- Warm construction: `CurrentCorrelatedTensorLocator(locator=checked_actual_physical_locator)`.

```python
from sympy import Rational
from lei_ren_part1_paper_compliant_current_correlated_radius_operator import SYMBOLS

field = CurrentCorrelatedTensorLocator(locator=checked_actual_physical_locator)
request = field.radius('Rref', '-7.5')
location = field.locate(request, z_phys='.7', log_tau='-2.6', theta='.41', viscosity='.8')
tensor = field.tensor(request, '.7', '-2.6', '.41', '.8')

# Exact microscopic phase, even though a numerical width enclosure touches 0:
micro = field.radius_expression(SYMBOLS['log_Ra'] + SYMBOLS['h_bridge']*Rational(1417,1000))

# A source offset box crossing a seam retains both legal regions:
crossing = field.radius('Rref', ('-7.1', '-6.9'))

# Only an exact source identity can select a common boundary trace:
seam = field.radius_expression(field.boundaries['gap_end']['expression'])
traces = field.traces(seam, field.locate(seam, '.7', '-2.6', '.41', '.8'))
```

Use exact SymPy expressions, integers or decimal strings. Float source expressions, unrelated symbols, unbound/nonfinite offsets, invalid log/reciprocal domains, foreign source requests and forged trace selectors are rejected. Candidate selection scans the entire original inventory; the expected region is never supplied to the locator.

## Exact radius and physical coordinate semantics

The original physical map remains

\[
\lambda^2-(z_{\rm phys}^2/\nu)\lambda^{2\delta}=\tau,
\qquad
\log r_{\rm phys}=(\log2+\log\nu+2q+\log R_{\rm source})/2,
\quad q=\log\lambda.
\]

Actual q is separately enclosed using the already admitted directed implicit map. Actual log(tau) is retained throughout. Every complete native signed contribution keeps its original source row and obtains the separate actual-lambda prefactor; native sqrt(tau) bounds are not silently reinterpreted as off-axis lambda.

The exact width is `h_bridge=h_switch=exp(exact_source_log_hb)>0`, bound to the same checked finite log datum. This defining expression establishes strict positivity. Its numerical nonnegative cap remains an enclosure only. Shared source aliases and local widths cancel before inverse evaluation. Exact rational logarithms from the source are normalized algebraically, so `log(10)+log(11)-log(110)=0` cannot leave an interval rounding residue larger than h. Larger user constants may remain unsimplified and conservatively ambiguous.

Crossing boxes preserve the candidate union. O2 axial logarithms retain their legal positive argument; the reciprocal pulse-gap slope is exactly `(1-4*mu)/mu>0` under the admitted `0<mu<1/4`. Axis uses the same analytic core. Finite but unmaterializable exterior source radii use the admitted unbounded identically zero Gamma T/E; no astronomical exponent or nominal midpoint field is needed.

Exact boundary selectors compare canonical source expressions, not interval overlap. Common trace rows are supplied by the checked full registry. Their native sqrt(tau) triangle bounds remain conservative; the actual q is reported separately, rather than treating the aggregate bound as an individual signed contribution.

## Focused evidence

- All 33 independently specified source targets select one region; source expressions are independently written from construction anchors, not native forward/inverse round trips.
- The same axial-restoration source represented as an absolute enclosure selects 25 candidates; its exact correlated request selects one. This comparison changes the input semantics explicitly.
- All 46 exact source radius selectors retain the expected neighboring regions and complete trace route. Each of the 32 adjacent boundaries is probed on both sides by a strictly positive displacement `h/10^30` (64 unique one-sided results).
- Five source offset boxes cross microscopic, restoration, pulse and both O2 axial seams; all expected unions are retained without asserting an exact boundary identity.
- Four fresh complete T/E views: 284 component groups and 4,225 retained signed contributions, including the microscopic second bridge, axial restoration, reciprocal pulse-gap end and an unmaterializable exact Gamma radius. Gamma contributes zero nonzero signed rows.
- Three fresh full trace evaluations: gap_end, patch_support_49 and outer_angular_-3_-1, each with 71 common component groups.
- Exact axis extension and invalid/foreign source guards are exercised.
- Focused controller reuses the fresh receipt only after matching complete source hashes and manifest, and preserves the admitted graph.
- Read-only mathematical/source review uses GPT-5.6 Luna / max.

Only `current_correlated_source_radius_inverse_locator_available` and `current_46_exact_source_radius_boundary_trace_selectors_available` are newly admitted. Earlier native-registry and actual-coordinate gates are inherited through checked constructors. All broader physical cover/cone/lift/point/NS/energy/flatness/temporal-recursion gates remain false.

## Next tasks and completion criteria

- [x] F57C-global-cover1e-2b: retain exact common radius anchors/local offsets on the same family, datum and original radius program; keep exp(logh) provenance explicit.
- [x] F57C-global-cover1e-2c: original 33 source-region inverses with source cancellation before enclosure; exact small phases, conservative crossing unions and complete signed T/E routing.
- [x] F57C-global-cover1e-2d, **source-dependent point-family scope**: bind all 32 adjacent / 14 support source recipes to the admitted full trace registry. General fixed-coordinate boundary certification stays open.
- [ ] **Next F57C-global-cover1f:** construct an independent finite-positive-time physical cover argument. Prove unique actual implicit scale for finite z and positive tau/nu; prove strict ordering/monotonicity of every original source radius under the actual audited source hypotheses; include the exact axis and unbounded Gamma. Cover arbitrary physical inputs and document unresolved interval unions. Publish a new scoped source-bound receipt; do not promote locator fixtures to a global theorem.
- [ ] F57C-cone1a: extract the original admissibility inequalities and exact stress lift from the supplied papers. Bind them to this graph's signed T/E and pressure history, with a clear comparison of physical/construction conventions.
- [ ] F57C-cone1b: compute correlated regional margins on the same source family; retain failing signs as construction targets. An interval upper bound on absolute stress is not a positive cone margin.
- [ ] F57D-recursion1a/b: implement actual n=1 recovery equations, common inner domain and independent moment repair; cancel or absorb the certified leading axis E term through genuine corrections. Its unchanged value grows as tau^((-3+delta)/2), so it is not flat.
- [ ] F57D-recursion2a/b,3: distinct n>=2 operators and at least two nontrivial orders, finite-order remainder estimates and final smooth summation. Coordinate rescaling alone is not coefficient recursion.
- [ ] F57C-points1/energy1: converged physical u/v/w/p, independent point error controls and prescribed-domain finite energy.
- [ ] F57E/F: both oscillatory pulse families/mean corrections, averaged quadratic stress cancellation, corrected Cartesian NS and measured contraction/slenderness/material winding.

The long-term goal remains active. Reuse the warm admitted graph and proceed to the independent physical cover construction.
