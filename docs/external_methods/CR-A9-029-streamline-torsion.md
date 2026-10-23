# CR-A9-029 — streamline torsion / handedness profile

## External method screened

- Source repository: `scipy/scipy`
- Screened commit: `27923e3c62756666088c4b2384daf9789be6f558`
- Public API: `scipy.interpolate.make_interp_spline`
- Source file: `scipy/interpolate/_bsplines.py`
- License: SciPy BSD-3-Clause terms (`LICENSE.txt`)
- Classification: **direct migration / public API only**
- Copied upstream implementation: none
- Dependency delta: none; this repository already requires `scipy>=1.10,<2`

SciPy supplies an interpolating B-spline and its derivative evaluation. This module uses an exact quintic spline (`k=5`) on the cumulative input-polyline-length parameter so first, second, and third derivatives of an already-integrated 3-D streamline can be evaluated from one interpolant. No smoothing, optimizer, registration, camera model, or velocity-field reconstruction is imported from SciPy.

## Repository-local scope

`measure_streamline_torsion_profile(...)` computes the intrinsic curve torsion

`tau = dot(cross(r', r''), r''') / ||cross(r', r'')||^2`

on one supplied finite 3-D polyline. It reports signed/absolute/RMS/max torsion, integrated signed/absolute torsion, valid-arclength fraction, and a handedness-coherence diagnostic. Only machine-scale near-zero-curvature denominator samples are excluded; their omitted arclength fraction is reported rather than silently imputed. The returned arrays are immutable.

The diagnostic is invariant to rigid translation/rotation and curve reparameterization/reversal, while reflection reverses the signed torsion as expected for intrinsic chirality. It deliberately does not fit away reflections: a mirrored helix is a different handedness observation.

## Difference from existing repository work

Existing streamline work already integrates trajectories, measures axis-specific turns/pitch, checks integration-step convergence, and measures curvature. This increment does not reintegrate streamlines or repeat curvature. It adds the missing intrinsic 3-D twist/chirality channel, which can distinguish two curves with similar curvature or projected pitch but different out-of-plane helical geometry.

## Direct contribution to final `[u,v,w]`

For frozen candidate fields evaluated under the same seed/time/integrator contract, torsion gives a scalar/profile diagnostic for whether the true 3-D streamline geometry has a coherent helical handedness and how strongly it twists. That can help rank which candidate deserves final MATLAB/Python/VTK visualization without changing the field or relying only on camera-dependent appearance.

## Truth boundary and limitations

This is visualization geometry only. It does not infer OpenAI hidden velocity values, hidden frame times, a camera, or a vortex-line identity. It does not alter `velocity`, pressure, forcing, energy normalization, residual definitions, validation samples, or any acceptance threshold. Attractive or coherent torsion is not Navier–Stokes validation and cannot promote `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `paper_exact`, `openai_field_identified`, or `blowup_proved`.

Torsion uses third curve derivatives and is consequently more sensitive than curvature to streamline integration error and polyline sampling. Candidate-to-candidate comparisons must therefore hold seed/time/integration/sampling contracts fixed and should be paired with existing integration-step and frozen-grid resolution audits. No universal torsion or handedness pass threshold is introduced here.
