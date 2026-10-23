# CR-A9-025 — streamline curvature profile

- source repo: `scipy/scipy`
- screened commit: `f0371a854cc7031df9a1ff920c304979b8c7d93e`
- public API: `scipy.interpolate.CubicSpline`
- source area: `scipy/interpolate/_cubic.py`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none (`scipy>=1.10,<2` is already a project dependency)

The repository uses exact natural cubic interpolation only to obtain first and second derivatives of an already-computed 3-D streamline as a function of cumulative polyline arclength. It adds no smoothing parameter and does not integrate the streamline itself.

Repository-local logic owns input validation, arclength parameterization, curvature `||r' x r''|| / ||r'||^3`, total turning, tortuosity, immutable outputs, provenance, and fail-closed truth-state metadata.

This metric is a visualization-geometry diagnostic. It does not fit camera pose, register to public imagery, infer hidden OpenAI velocities/times, modify `[u,v,w]`, or substitute for PDE/divergence/support/energy validation. Curvature depends on the supplied streamline integration and sampling quality, so comparison should use the same integration and resampling contract across candidates.
