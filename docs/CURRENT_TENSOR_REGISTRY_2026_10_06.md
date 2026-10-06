Historical native-registry milestone. Continuation: [CURRENT_PHYSICAL_TENSOR_LOCATOR_2026_10_06.md](CURRENT_PHYSICAL_TENSOR_LOCATOR_2026_10_06.md); actual implicit coordinates and original radius inverse candidates now exist, while exact physical seam/global cover remain open.

# One-graph native full tensor registry (2026-10-06)

Implementation and scoped one-graph native tensor registry receipt: commit [cad58441](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/cad58441fbe661b08695f6bfbdee0a3d6fc93fcf).

All **33 existing full tensor regions**, **32 adjacent traces** and **14 internal support traces** now have one source-bound callable registry. It reuses the checked nonlinear core/axis/angular/pressure/energy/history graph through the complete Gamma exterior. Every native view and its original signed contributions remain available.

The 14 support traces are exactly six patch, four pulse-end and four outer-angular boundaries. The core axis is a separate Cartesian analytic extension of the same region; it adds no region or counted support boundary. The primitive velocity/pressure atlas remains 14 adjacent / 8 internal.

This completes the native registry layer. An automatic physical-coordinate locator and global physical cover certificate are still open. The returned values are directed enclosures of source functions, not converged physical point values. Original admissibility/lift, actual coefficient recovery, completed temporal flatness, prescribed energy and corrected NS remain open.

## Implementation and API

- `experiments/root_st073/lei_ren_part1_paper_compliant_current_tensor_registry.py`.
- Complete registry/source bindings: matching `.json`.
- Independent checker/receipt: matching `_check.py` and `_check.json`.
- Focused controller: `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currenttensorregistry`.
- Warm construction: `CurrentTensorRegistry(background=checked_current_angular_internal_background_tensor)`.

```python
registry.native('pulse_gap_end', Z, phase, log_tau, theta, viscosity)
registry.native('core_positive_radius', Z, rho, log_tau, theta, viscosity)
registry.core_cartesian(Z, X, Y, log_tau, viscosity)
registry.trace('gap_coordinate', Z, log_tau, theta, viscosity)
registry.trace('outer_angular_-3_-1', Z, log_tau, theta, viscosity)
registry.exterior(Z, log_tau, theta, viscosity)
```

Native coordinates must match the named region. `rho=0` routes through the nonsingular axis source. Requests spanning zero should use `core_cartesian`; the original positive-radius cylindrical chart retains its domain guard. Z endpoints represent source limits at physical infinity.

`native` and `core_cartesian` return the complete original view plus `canonical_signed_component_groups`. Most regional views have 71 common groups; postpulse zero meridional remainder rows are added as exact zero source rows without replacing any nonzero contribution. The Cartesian core uses its explicit derivative-index layout, with 210 groups. `trace` keeps `original_complete_trace`, including the source equality theorem and all common bounds; it aliases the one different angular-entry common-row key. Bounds are not substituted for fields.

## Original owners and coordinates

| Native region(s) | Same-graph callable | Coordinate |
| --- | --- | --- |
| core_positive_radius | core.chart('core', ...) / axis.axis / axis.cartesian | positive rho in (0,4]; nonsingular Cartesian disk rho in [0,4] |
| bridge_first / bridge_second / bridge_macro | bridge.chart | phases [0,1], [1,2], macro coverage fraction [0,1] |
| switch_first / switch_second / switch_power | micro.chart / power.chart | phases [0,1], [1,2]; power fraction [0,1] |
| reshape | reshape.chart | phase [0,1] |
| inner_reference / axial_restore / restore_buffer | restore.chart | phases [0,1], [0,1], offset [-7,-6] |
| actual_patch | patch.chart | x=R/Rm in [1,e]; exact named support/Rh limits |
| Rh_reference / O2_slope / O2_axial / O2_buffer | o2.chart | [-5,0], [0,1], phase [0,1], buffer offset [0,11] |
| O3_slope_mu / O3_power | o3.chart / incoming.chart | offset [0,1], phase [0,1] |
| pulse_entrance / pulse_main / pulse_exit | incoming.chart / main.chart | xi [0,.02], [.02,10], [10,11] |
| pulse_gap / pulse_gap_end | gap.chart | xi [11,12]; phase [0,1], exact reciprocal s=-d/mu |
| pulse_end | end.end | s [-4,0] |
| flatten / outer_power | flatten.chart | offset [0,100]; phase [0,1] |
| outer_angular | angular.angular | s [-4,0] |
| steep_entry | entry.entry | offset [0,1] |
| steep_power / steep_exit / waiting | o7.chart | [0,1] |
| heat_collar / heat_exterior | heat.chart / heat.unbounded_exterior | [0,3], >=3 / full unbounded Gamma |

Each registry entry binds the exact original method signature, full method AST and file hash. Per-region `source_coordinate` and `derivative_pullback` describe the coordinate contract. Core uses D_y=rho*d_rho with Euler/Stirling conversion; bridge first/second use hb^-j phase conversion; macro receives no extra inverse-length factor. The O2 axial phase and buffer_offset remain distinct in the original route. Gap-end uses its exact reciprocal source adapter. Other coverage phases/offsets are handled by their original operators, with no derivative rescaling in the registry.

The 32 adjacent records pair the ordered native regions and bind the actual admitted receipt row for each named seam. Exact endpoint coordinates are supplied by the pinned original `interface` recipe. The 14 internal records similarly bind each individual support row and its exact endpoint recipe. Source equality is inherited from these admitted original proofs; a connected list alone does not prove physical-coordinate cover.

## Targeted evidence

- One fresh complete view for every one of the 33 regional routes.
- Separate exact axis, nonzero Cartesian disk, and unbounded Gamma views: 36 views total, 2834 groups / 25513 retained signed contributions.
- All original callable AST/signature bindings independently replayed; shared family/source/datum and exact owner paths enforced.
- Fresh exact reciprocal gap, end/flatten, angular/entry, patch support, pulse support and angular support trace routes.
- Original factor schemas checked; invalid native/time/viscosity/axis-angle inputs and foreign owner/route rejected.
- Source metadata and routing reviewed read-only by GPT-5.6 Luna / max.
- Focused controller uses the real fresh checker receipt after complete source-hash/manifest matching; no historical producer rebuild or redundant second numeric routing run.
- Three changed Python files compiled, staged whitespace and working/index source hashes checked before publication.
- Gates `current_33_region_full_tensor_native_registry_available` and `current_32_adjacent_and_14_internal_tensor_trace_routes_bound` true. Both physical locator/cover gates and all global cone/NS/energy/point/temporal-recursion gates remain false.

## Completed tasks and next construction

- [x] F57C-global-cover1a: 33 native source-owner/domain records from one checked graph.
- [x] F57C-global-cover1b: each of the 32 adjacent and 14 support traces mapped to original endpoint recipe and actual individual receipt evidence; primitive atlas kept separate.
- [x] F57C-global-cover1c-native: all regional native routes, exact axis, reciprocal gap, microscopic bridge/switch and unbounded Gamma APIs available.
- [x] F57C-global-cover1d-native: full views and signed factor sectors retained with common regional groups and explicit Cartesian core indices.
- [ ] F57C-global-cover1e-1: recover and bind the actual implicit similarity/physical lambda relation and constant-viscosity pullback from the existing physical mapper. Implement a directed physical-coordinate locator in log radius; avoid materializing huge/small scale factors or replacing source intervals with their midpoints.
- [ ] F57C-global-cover1e-2: invert each original regional radius formula to its native phase/offset, resolve exact boundaries using the admitted trace recipes, and return a source enclosure with its locator error. Respect the full Cartesian core at rho0 and finite physical |Z|<1.
- [ ] F57C-global-cover1f: independently verify physical radial coverage from core through all 32 boundaries to the unbounded Gamma exterior, with no missing domain or wrong coordinate factor. Keep this a finite requested positive-time cover certificate.
- [ ] F57C-cone1a/b: recover the paper's stress cone/lift conditions for this same signed tensor; compute correlated margins in each required region. Retain failed margins as construction targets.
- [ ] F57D-recursion1a/b: derive and recover the actual n=1 coefficient functions, pressure/stress and per-order moment repair on the same core interval. Remove or absorb the certified nonzero leading origin term E_z=sqrt(nu)*k_axis*tau^((-3+delta)/2); it cannot be flat unchanged.
- [ ] F57D-recursion2a/b: distinct n>=2 operators/nonlinear forcing and at least two nontrivial orders, with streamfunction/vector-potential localization before curl and matched outer histories.
- [ ] F57D-recursion3: finite-order remainder estimates and actual smooth coefficient sum; certify temporal flatness only for the completed corrected background.
- [ ] F57C-points1 / energy1: convergent same-source physical u/v/w/p with error control, and finite energy on the prescribed localization/domain.
- [ ] F57E-oscillation1 / corrected1: both oscillatory families/mean corrections, averaged flux cancellation and independent corrected Cartesian NS residual.
- [ ] F57F-dynamics1: measured contraction, relative axial slenderness, swirl/vorticity growth and material winding, separately from genuine coefficient recursion.

The long-term reconstruction goal remains active. Continue construction using checked inputs and publish achieved scope only.
