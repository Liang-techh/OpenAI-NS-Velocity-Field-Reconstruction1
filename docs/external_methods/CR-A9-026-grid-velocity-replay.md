# CR-A9-026 — truth-bounded regular-grid velocity replay

## Goal

Make an already sampled `(time,x,y,z,component)` Cartesian velocity artifact directly callable again as `velocity(x,y,z,t)->[u,v,w]` for Python visualization and downstream export, without requiring the original analytic candidate implementation.

This is a delivery/visualization adapter only. Interpolation error is not a PDE residual estimate, and interpolated derivatives are not accepted for PDE/divergence promotion.

## External result screened

- source repository: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API: `scipy.interpolate.RegularGridInterpolator`
- source file: `scipy/interpolate/_rgi.py`
- license: BSD-3-Clause terms in SciPy `LICENSE.txt`
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none; this repository already requires `scipy>=1.10,<2`

The SciPy API supplies interpolation on a rectilinear grid and supports fail-closed out-of-bounds evaluation through `bounds_error=True`. This repository fixes the method to `linear`, requires strictly increasing physical axes, forbids extrapolation, preserves `[u,v,w]` component order, adds immutable identity/provenance, rejects an exact all-zero grid, and owns the governed NPZ save/load schema.

## Migration difference / truth boundary

The adapter does not reconstruct hidden OpenAI samples, infer frame times, fit a camera, choose a candidate, alter the sampled velocity, or claim interpolation is an independent PDE checker. `velocity_export_ready=true` means only that the frozen grid can be replayed as a callable visualization field. `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, and `blowup_proved` remain false.

## Direct contribution

A frozen candidate can now follow the delivery chain

`analytic candidate -> sampled (t,x,y,z,[u,v,w]) grid -> governed NPZ -> callable GridVelocityReplay -> Python streamlines/vorticity / further MATLAB-VTK export`.

This reduces dependence on candidate-specific analytic code after the final field is frozen and makes a saved grid usable as a stable visualization interchange artifact.
