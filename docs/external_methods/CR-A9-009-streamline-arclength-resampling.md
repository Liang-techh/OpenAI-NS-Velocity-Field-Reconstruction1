# CR-A9-009 — arc-length streamline resampling

## Current blocker

The velocity-first path already has candidate-agnostic streamline integration and deterministic seed planning in open work. Dense 3D display still has a practical geometry/rendering issue: integrators may emit nonuniform point spacing, while per-segment colored rendering can create tens of thousands of graphics objects for roughly 200 lines. A fixed point count on a smooth arclength parameterization makes cross-time geometry easier to compare and keeps rendering bounded to one object per streamline.

## External source screened

- Repository: `scipy/scipy`
- Commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- Public API: `scipy.interpolate.PchipInterpolator`
- Relevant source: `scipy/interpolate/_cubic.py`
- License: BSD-3-Clause
- Classification: **directly reusable public API**

SciPy describes PCHIP as a C1 piecewise-cubic shape-preserving interpolator that requires a strictly increasing independent coordinate and avoids overshoot for monotone data. This increment uses only that public API; no SciPy interpolation implementation is copied.

## Migration scope and difference

The repository supplies the arclength coordinate itself from each input polyline, collapses only exact consecutive duplicate samples, and uses PCHIP independently on Cartesian x/y/z. Endpoints are restored exactly. The module additionally supports height or speed coloring and records polyline segment-spacing variation; these are repository visualization concerns, not SciPy behavior.

No new dependency is added because SciPy is already required by the project.

## Truth boundary

This changes only the rendered polyline representation. It does not alter or refit `velocity(x,y,z,t)`, pressure, forcing, domain, support, viscosity, optimizer, residual, validation samples, or thresholds. Smooth/resolution-stable streamlines are visualization evidence only and cannot establish PDE validity, public OpenAI correspondence, hidden-field recovery, paper exactness, singularity, or blow-up.
