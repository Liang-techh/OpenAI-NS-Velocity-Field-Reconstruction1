# CR-A9-006 — structured-grid velocity interpolation

## Classification

**Directly reusable public API.** This increment uses SciPy's public `scipy.interpolate.RegularGridInterpolator` API. No SciPy interpolation implementation is copied into this repository.

## External source

- Repository: `scipy/scipy`
- Screened commit: `eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`
- Relevant source: `scipy/interpolate/_rgi.py`
- License: BSD-3-Clause (`LICENSE.txt` at the screened commit)
- API behavior relied on: rectilinear N-D grids, linear interpolation, strictly monotone coordinate axes, and `bounds_error=True` for fail-closed out-of-domain evaluation.

## Why this lowers the current delivery blocker

The project goal is a callable, save/loadable 3D time-dependent velocity field usable by MATLAB/Python visualization. Existing/open work provides the production candidate lane, direct-grid consistency audits, and interchange formats, but a saved structured velocity grid is not itself a callable field. `StructuredVelocityField` closes that delivery gap by turning one frozen sampled grid back into the common `velocity(points,time)->[u,v,w]` interface.

This is intentionally distinct from PR #27: that PR audits trilinear visualization-grid error against direct candidate calls; this module packages and reloads a frozen rectilinear space-time field and delegates interpolation to SciPy's independently maintained public implementation.

## Migration scope and differences

Migration scope is API-only:

- sample a frozen velocity callable on declared `x,y,z,t` axes;
- store `[t,x,y,z,component]` values in a pickle-free NPZ bundle;
- reload the bundle and evaluate it with linear `RegularGridInterpolator`;
- reject extrapolation, malformed axes, bad shapes, non-finite velocity values, or claim-state tampering.

Differences from SciPy: this repository adds only the velocity-field schema, save/load wrapper, candidate identity, and truth-boundary metadata. It does not reimplement interpolation mathematics.

## Regression / interface smoke

Focused calibration uses a fully space-time affine manufactured velocity. Linear interpolation must reproduce it at 257 random off-grid points and non-grid times to floating-point tolerance. A save/load round trip must preserve the callable interface and grid output. Mutation checks cover duplicate coordinates, spatial/time extrapolation, NaN source velocity, and attempted promotion of `pde_validated=true`.

## Truth boundary

Interpolation/save/load success means only that a frozen sampled field can be delivered and visualized reproducibly. It does **not** establish Navier–Stokes residual acceptance, divergence/support validity, visual correspondence to OpenAI public imagery, hidden-field identification, paper exactness, or blow-up. Interpolation error must remain separate from independent PDE validation, and no acceptance threshold is changed by this module.
