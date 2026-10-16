# Complete first-micro cumulative m15 - 2026-10-16

`OriginalFirstMicroFullM15` extends the cumulative signed leading pressure moment from the previously accepted s=1/4 endpoint through the entire original first_micro chart, ending at s=1. It preserves the six core/cutoff/lower inlet sectors and the original nonzero singular pressure history, and carries complete factored selected-N correction errors over every traversed source cell. The earlier quarter-slice implementations and reports remain unchanged.

The full goal is **ACTIVE / INCOMPLETE**; tasks 14/15 remain IN PROGRESS. Full first-micro leading pressure transport with correction errors is available. Resolved finite-N signed centers, useful total relative cancellation, second_micro and the full active/compact/Rp join remain open.

## Source domain and unchanged leading equations

The actual prescribed first-micro source is defined for s in[0,1], with R=4epsilon*exp(hb*s), dy=hb*ds and chi=1-sigma(s)+hb*sigma(s). Its formulas do not change at s=1/2. The narrower half-chart admission in the original n1 velocity transport module is not invoked to extrapolate a raw n1 solution. Only its leading-field preparation arithmetic is reused, with original signed first-chart field increments and complete errors; the cut n1 velocities have already stopped, and their cumulative pressure history is retained.

Source evidence is in `current_original_micro_functions`, `current_original_prescribed_bridge_C6_functions`, `compliant_actual_bridge_integrals` and the first-chart [0,1] admission of `OriginalFirstCollarJets`. The normalized regular pressure and the separate inverse-hb/d-sigma split remain the same as in the [accepted quarter-slice integral](CURRENT_ORIGINAL_FIRST_MICRO_M15_2026_10_16.md).

The complete actual selected-N source partition is

```
[sc/2,1/4], [1/4,1/2], [1/2,3/4], [3/4,1].
```

The constructor checks every original axial-frame partition before evaluating. Regular pressure quadrature splits at all three interior quarter boundaries, retaining directed cell/source errors and all original physical factors.

## Uniform correction errors over the complete prefix

An endpoint's local B cap alone cannot bound the earlier volume integral once that integral crosses source cells. The new row owner therefore unions all source cells in the positive-kernel history prefix with all terminal local cells, and uses this union for the uniform finite-N U error. Mean errors still integrate the original partial cell widths; they are not reset at a quarter boundary.

Deduplication removes repeated references to the same object only. Distinct axial-frame cells with the same radial boundaries remain present. Q/Q_Z/Q_ZZ error recovery then includes these complete U/m histories. U_y/Q_y/Z1 errors retain the same original bridge_first slow B_y/N and undivided B_phi/B_Z_phi bounds and physical Pstar factor.

The inherited reduced-S1 boundary, volume and finite RQ_y/RQ_yZ error formulas can therefore be used over the complete first-micro prefix. Nonlinear delta-delta cross terms, actual selected integer, phase origin and physical R^2*hb*ds measure remain present. Symmetric correction sectors are absolute-error enclosures; their proxy zero does not assert a selected signed correction.

## Finite sigma endpoint and singular memory

The upper half-chart uses the original sigma_jets value directly as the singular coefficient, with a zero log atom. At exact s=1, sigma is exactly1. This avoids dividing by (1-s)^2 or forming an infinite odds times a vanishing coefficient. The lower sigma(sc/2) remains its nonzero original log atom and is subtracted once.

The leading singular contribution is still epsilon^2*A0*(sigma(s)-sigma(sc/2)), plus its complete positive d-sigma width error. Even though sigma-prime vanishes at s=1, this accumulated singular pressure moment does not vanish or restart. The output retains the same thirteen-sector partition as the quarter-slice API.

## Focused evidence and continuation

The focused driver reuses the actual Laurent-source and reduced-S1 physical value/Z1 identities, checks exact accepted quarter-history replays, all original prefix/axial-frame cells, serialized increment errors, the exact upper sigma=1 endpoint, empty-prefix zeros, lower singular history and cell refinement. It also evaluates an interior point in the final source cell. The focused run exited 0 in 96.177 seconds; saved-report inspection exited 0: 19 original packets, three exact accepted quarter replays, 56 prefix/axial-frame conditions, 168 serialized increment-error checks, six exact sigma=1 conditions, 42 empty-prefix rows, 12 retained lower histories, six refinement comparisons, 32 product-error sign cases and four rejected invalid queries. The 2,343,031-byte report SHA256 is `09c6509957f7bbbfc84cdc802d9dfbeb0110898de81e6de9fb3e924fa3fd72d0`; it pins 1751 inputs and 35 implementations.

Use `python scripts/run_original_first_micro_full_m15.py --inspect-report` for saved evidence; `--write-report` rebuilds after relevant implementation changes. Existing read-only worker `/root/n1_atlas_map`: **GPT-5.6 Luna / max**; no new agent. It confirmed the leading [0,1] source domain, whole-prefix correction union, endpoint treatment and retained error-only scope.

Next continue second_micro from the unchanged s=1 cumulative inlet, with its actual constant-hb leading controls, original absolute radius and incoming correction histories. Then complete remaining active charts, compact finite-N pressure and independent source-correlated Rp cancellation. Global five-moment numerical joins, common N, analytic preheat/exact heat/energy, admissible stress/flat remainder, genuine temporal n>=2 recursion, oscillatory correction and corrected Cartesian NS remain open.
