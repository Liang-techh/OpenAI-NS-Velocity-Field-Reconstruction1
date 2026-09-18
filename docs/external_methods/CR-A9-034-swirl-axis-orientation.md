# CR-A9-034 — swirl-axis orientation external-method screen

## Source and classification

- Source repository: `guilindner/VortexFitting`
- Audited source commit: `083e25850594c8e5efdd2e38a85d6f0c8f8d026d`
- Relevant source: `docs/methodology.rst`, section **Swirling strength criterion**
- License: MIT (`LICENSE` at the audited commit)
- Classification: **suitable for minimal reimplementation / mathematical definition only**
- Copied upstream code: **none**
- New runtime dependency: **none**

The upstream methodology states that a 3-D velocity-gradient tensor with one real eigenvalue and one complex-conjugate pair can be written using the real-eigenvalue eigenvector `nu_r`; the imaginary part `lambda_ci` of the complex pair is the local swirling strength.  CR-A9-034 uses only that public mathematical decomposition: `lambda_ci` weights the local orientation, while the eigenvector belonging to the remaining real eigenvalue supplies an unsigned local swirl axis.

## Why this is useful for the repository goal

`docs/VISUAL_TARGET.md` records the public OpenAI announcement image as a slender, vertically oriented, axis-centred vortex with inward/helical trajectories and axial stretching.  It also records that the paper's physical cylindrical coordinate uses the `z` direction as the axis.  A candidate can therefore be audited for whether its *intrinsic 3-D swirl axes* actually align with a caller-declared physical axis, instead of accepting a camera-dependent 2-D impression of verticality.

This lowers the chance that a visually convenient projection, shear-dominated structure, or arbitrary camera orientation is mistaken for the intended three-dimensional vortex geometry.  It does **not** recover the hidden camera, trajectory seeds, hidden frame time, or OpenAI velocity samples.

## Local implementation scope

`src/openai_ns_reconstruction/constrained_swirl_axis_orientation.py`:

1. Samples only the caller's public `velocity(points,time)->[u,v,w]`.
2. Builds all nine components of `grad(u)` with centered Cartesian second-order finite differences; it does not impose `div(u)=0` while differentiating.
3. Computes the full local eigensystem of `grad(u)`.
4. Uses only points with a numerically resolved complex-conjugate pair.  The imaginary-part guard is `64*eps*max(1,||grad_u||_F)` and is documented as floating-point noise suppression, not an acceptance threshold.
5. Uses the eigenvector belonging to the remaining real eigenvalue as an **unsigned** local axis.  Eigenvector sign is intentionally discarded, so the diagnostic cannot manufacture handedness.
6. Reports swirl-strength-weighted axis alignment, tilt quantiles, an unsigned orientation tensor/principal axis, orientation concentration, velocity RMS, and same-operator divergence RMS.
7. Leaves all scientific truth states false and defines no post-hoc visualization pass threshold.

## Deliberate differences from the upstream package

The upstream package contains vortex detection/localization choices aimed at experimental-flow analysis, including optional normalization and peak localization.  CR-A9-034 does **not** migrate those choices.  In particular it does not:

- normalize `lambda_ci` by a profile RMS;
- choose an empirical vortex threshold;
- detect peaks or fit a vortex model;
- replace a derivative component using an incompressibility identity;
- rotate/register the candidate to improve the reported alignment;
- fit a camera or image projection;
- change the candidate velocity, pressure, forcing, support, energy normalization, or residual threshold.

The output is a morphology diagnostic only.  A high axis-alignment score is neither PDE validation nor proof of correspondence to OpenAI's hidden field.

## Verification added in this PR

Focused analytic regressions cover:

- exact recovery of the axis and `lambda_ci` for rigid-body rotation;
- an oblique rotation axis and invariance of orientation under uniform velocity scaling;
- nonzero simple shear producing no false swirl axis;
- a deliberately wrong declared target axis staying visibly mismatched rather than being registered away;
- fail-closed handling of out-of-window time, out-of-box bounds, even/undersized grids, zero target axis, bad candidate identity/provenance, exact-zero velocity, malformed velocity shape, and non-finite velocity values.

Local focused test command used before upload:

```text
PYTHONPATH=/tmp/a9_034/src python -m pytest -q -W error /tmp/a9_034/tests/test_constrained_swirl_axis_orientation.py
```

Result: `5 passed in 0.10s`.  `py_compile` also completed successfully for the module and test file.
