# Current downstream physical maps — 2026-10-05

CurrentDownstreamPhysicalAssembly maps the fourteen accepted downstream charts through the original Cartesian and physical-time operators. It uses the current dispatcher's nested core as the parameter source and the exact same Rh-to-Rp history object. The adapter does not initialize the legacy global core owner.

Run on branch codex/st073-transition-next:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentphysicalmaps
```

```python
from lei_ren_part1_paper_compliant_current_downstream_physical_assembly import CurrentDownstreamPhysicalAssembly
field = CurrentDownstreamPhysicalAssembly()
packet = field.evaluate("O3_power", Z=".3", coordinate=1, log_tau="-2", theta=".4")
```

## Current ownership and coordinate coverage

The fourteen source owners are switch_first, switch_second, switch_power, reshape, inner_reference, axial_restore, restore_buffer, actual_patch, Rh_reference, O2_slope, O2_axial, O2_buffer, O3_slope_mu and O3_power. Existing chart phase selectors and domains are preserved.

For each chart the adapter provides all 35 Cartesian spatial multiindices of total order at most four, for ux, uy, uz and pressure, plus their first time derivative at fixed physical position. The source dispatcher, common pressure definition and parameter origins are current. Whole-domain reports use Z in [-1,1], the full angular bound and each original chart domain.

The original radius and normalized_sources callables are inherited unchanged. The original evaluate operator is explicitly called with the injected current dispatcher. Independent existing coordinate fixtures remain current by hashes: 140 Cartesian derivatives, four fixed-position time derivatives and 120 microscopic/macroscopic source-scale rows. Those fixtures certify reusable operators; they do not transfer the historical 33-chart admission to current sources.

## Physical units and source factors

The mapper retains the exact microscopic phase width and formal width powers before source bounds, pressure/swirl amplitude factors, the moving cylindrical basis and the time dependence of the implicit similarity coordinate. Huge positive exponentials stay in signed logarithmic source form.

For actual_patch, the input already contains derivatives of the full physical radial profile divided by the fixed Rm unit. Each derivative row is divided by sqrt(x0) at the evaluation point to convert to the fixed current R unit. This unit is frozen after differentiation. Applying another Leibniz normalization would remove the physical radial prefactor incorrectly.

A new independent finite patch fixture differentiates two full radial functions through total order four, checks 30 mixed rows, restores the physical amplitude and verifies that a constant shape still has a nonzero first radial derivative equal to one half in current units.

## Acceptance and scope

The focused checker verifies all fourteen owner receipts, current source/parameter graph, retained signed source rows, source-supported exact zeros and unfinished-scope flags. The producer records proof with certification false; runtime certification becomes true only after loading current_downstream_physical_assembly_check.json.

This delivers source-field enclosures and signed physical derivative bounds. It does not select nonlinear point coefficients or enclosure midpoints. The core/axis, bridge, pulse and heat charts remain outside this adapter's admitted scope. The external Rp-to-pulse source join remains open.

Production point fields, full shared leading-source/remainder admission, completed global stress tensor and admissibility, independently bounded flat remainder, required-domain energy, true n-dependent coefficient recursion, oscillatory stress cancellation and corrected dynamics remain open.
