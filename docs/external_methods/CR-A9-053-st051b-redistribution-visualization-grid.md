# CR-A9-053 — ST051-B frozen-redistribution visualization grid

## Task

One minimal Constrained Agent 9 delivery increment: turn the already-screened ST051-B + frozen inner/mid swirl redistribution child into a deterministic, self-describing 3-D time-varying `[u,v,w]` grid consumable from MATLAB or Python.

This is not a new basis, parameter search, image fit, pressure/forcing fit, or PDE acceptance test. It samples one immutable candidate-side transform so that the next fixed-seed 3-D streamline/vorticity rendering step does not have to reconstruct research branches by hand.

## Frozen internal sources

- repository base: `main@f0193d66c9d92948b4820ebcb70263673995b324`;
- ST051-B reconstruction: PR #460 exact head `4b784f1b8457af2ead49295631d834d4e882000b`;
- reported original ST051-B raw candidate SHA-256: `0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d`;
- frozen redistribution transfer: PR #469 exact head `95f94188270f02511c8bdf660822f77530c9983a`;
- transform identity inherited from PR #458: `alpha_C=2.520520814687742`, inner window `(0.30,1.05)`, outer window `(0.95,1.85)`, gain `.025`;
- PR #471 independently showed the same transform increases actual ST051-B cumulative material winding on the frozen 48-path contract without sacrificing aggregate contraction or axial-pair separation.

Classification for the PR #469 transform: **direct internal method reuse / independent reimplementation**. Migrated scope is only the frozen radial windows, `alpha_C`, gain, azimuthal-only velocity transform, and one positive common normalization restoring the parent reference energy. It is not rebalanced or refitted on this branch.

## External method

The grid writer uses the public `scipy.io.netcdf_file` API to create NetCDF3 64-bit-offset files. Provenance is pinned to `scipy/scipy@b12c772edbc1fe0d3db9481cdfcc2e311569cb25`; license **BSD-3-Clause**. Classification: **direct migration / public API only**. No SciPy implementation is copied and no new dependency is introduced beyond the repository's existing SciPy use.

NetCDF is chosen as a cross-language carrier rather than as a numerical method. Python can read it with SciPy-compatible NetCDF readers; MATLAB can use `ncinfo`/`ncread` directly.

## Frozen export contract

Default export is `33^3` Cartesian points over the registered box `[-2,2]^3` at `t=.25,.50,.75`. Variables are

- `time`, `x`, `y`, `z`;
- `u(time,x,y,z)`, `v(time,x,y,z)`, `w(time,x,y,z)`;
- `speed(time,x,y,z)`.

The file records exact parent/transform heads, component/dimension order, a deterministic hash over coordinates and `[u,v,w]`, and fail-closed truth-state metadata. Reload verifies axes, shapes, finite/nonzero values, speed consistency, metadata and checksum before the artifact is accepted.

The exact workflow also checks the independently recomputed common reference-energy scale against PR #469's `1.0014791925672812` before export. Drift in the pinned ST051 checkout, parent recipe/truth state, transform constants or normalization aborts the build.

## MATLAB / Python consumption

MATLAB example:

```matlab
x = ncread('ST051-B-frozen-redistribution.nc','x');
y = ncread('ST051-B-frozen-redistribution.nc','y');
z = ncread('ST051-B-frozen-redistribution.nc','z');
t = ncread('ST051-B-frozen-redistribution.nc','time');
u = ncread('ST051-B-frozen-redistribution.nc','u');
v = ncread('ST051-B-frozen-redistribution.nc','v');
w = ncread('ST051-B-frozen-redistribution.nc','w');
```

Python example:

```python
from scipy.io import netcdf_file
with netcdf_file('ST051-B-frozen-redistribution.nc', 'r', mmap=False) as f:
    u = f.variables['u'].data.copy()
```

These examples only load the sampled candidate. They are not visual-correspondence or PDE tests.

## Constraint governance / truth boundary

Unchanged registered scientific contract: `nu=.01`, physical domain `R^3`, evaluation box `[-2,2]^3`, smooth zero extension outside `r<2, |z|<2`, time `[.25,.75]`, restricted two-parameter forcing family only, `E(.25)=1`, separate optimization/validation data, divergence max/L2 `1e-5`, momentum max/L2 `1e-3`.

The velocity transform changes the candidate and therefore **does not inherit ST051-B pressure, forcing or held-out momentum evidence**. This increment performs no pressure reconstruction and no held-out PDE residual. It does not use `u->0`, free `f=R(u,p)`, hidden OpenAI data, an image-derived numerical target, or a visual pass threshold.

The NetCDF may truthfully mark `velocity_export_ready=true` for this candidate-specific grid, while all of the following remain hard false: `production_candidate_selected`, `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, `blowup_proved`.

Direct contribution to final `velocity(x,y,z,t)->[u,v,w]`: the current best visualization-side ST051-B swirl child becomes a checksum-verified cross-language 3-D velocity artifact. The next useful Agent-9 step is fixed-seed streamline/vorticity rendering/comparison on this exact child, not another exporter or swirl basis.
