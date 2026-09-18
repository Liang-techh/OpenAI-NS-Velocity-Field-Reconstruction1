# CR-A9-045 — ST048-S governed visualization grid

## Scope

Constrained Agent 9 minimal delivery increment. This round does not add a new basis, optimizer, pressure/forcing direction, PDE metric, pathline metric, vorticity metric, or visual pass score. It turns the frozen ST048-S challenger from PR #390 into one exact-source, cross-language Cartesian `[u,v,w]` grid artifact that MATLAB and Python can consume directly.

## Fresh audit and non-duplication

Before implementation, Agent 9 reread current `README.md`, `AGENTS.md`, `docs/PROJECT_GOAL.md`, `docs/CURRENT_CHECKPOINT.md`, `docs/AGENT_TASKS.md`, audited `configs/constraints*`, current open PRs, recent commits and CI.

Relevant occupied lanes remain separate:

- PR #397 owns the frozen 48-material-path comparison for ST048-S/B;
- PR #398 owns target-free ST048-S/B 3-D vorticity-morphology transfer;
- older PRs own generic streamline/pathline integration, MATLAB export, VTK export, frozen-grid replay and render packing;
- PR #390 owns construction and independent PDE replay of ST048-S/B.

PR #398 routes ST048-S as the cleaner visualization-oriented residual challenger for the next fixed 3-D visualization. The missing seam addressed here is candidate-specific delivery: a deterministic ST048-S velocity grid tied to the exact source head and raw candidate identity, not another generic exporter or diagnostic.

## Frozen source

- source repository: `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`
- source PR: #390
- source head: `97695a86f85ce68fb4ae70c41fc81c904d655183`
- candidate: `ST048-S`
- parent: `ST047-E`
- archived parent raw SHA-256: `dfe6e51af93c9c42f1322798b82f89523e6a193c1b6495894cb60f63b4eb07e0`
- archived ST048-S raw SHA-256: `6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09`

The readable PR #390 recipe reconstructs the mathematical field. Regenerated JSON metadata is not claimed byte-identical to the archived raw candidate.

## Canonical constraint boundary

No registered scientific value changes. The contract remains:

- `nu = 0.01`;
- physical domain `R^3`;
- evaluation box `[-2,2]^3`;
- smooth zero extension outside `r<2`, `|z|<2`;
- time interval `[0.25,0.75]`;
- only the preregistered restricted two-parameter forcing family, with no free/residual-defined force;
- `E(0.25)=1 +/- 0.001`;
- optimization and held-out validation remain separate;
- original divergence and full-momentum gates remain unchanged, including `1e-5` divergence and `1e-3` momentum thresholds.

ST048-S still fails the original full-momentum `1e-3` acceptance gate in PR #390. Grid/export success does not alter that result.

## Minimal delivery contract

The artifact freezes:

- times `0.25, 0.50, 0.75`;
- full registered Cartesian box `[-2,2]^3`;
- exact `33^3` CI artifact resolution;
- dimensions `(time,x,y,z)`;
- separate NetCDF variables `u`, `v`, `w`, `speed`;
- deterministic grid SHA-256 over source identity, axes, times and all `[u,v,w]` bytes;
- fail-closed source head, candidate raw SHA, component order and false truth states.

The CLI can generate another odd resolution in `[5,129]`, but the checked CI artifact is fixed at `33^3`. Changing resolution creates a different grid SHA and is representation/delivery evidence only.

MATLAB can consume the artifact with `ncinfo` / `ncread`; Python can use `scipy.io.netcdf_file` or another compatible NetCDF reader. Derivatives from the sampled/interpolated grid are explicitly not PDE-acceptance evidence.

## External result screened / migrated

- source repo: `scipy/scipy`
- screened commit: `b12c772edbc1fe0d3db9481cdfcc2e311569cb25`
- public API: `scipy.io.netcdf_file`
- upstream source area: `scipy/io/_netcdf.py`
- license: BSD-3-Clause
- classification: **direct migration / public API only**
- copied upstream implementation: none
- dependency delta: none (`scipy>=1.10,<2` is already required)

SciPy supplies only NetCDF serialization/deserialization. Candidate reconstruction, exact-source pinning, grid sampling, checksum semantics, truth-state checks and provenance are repository-local.

## Truth boundary

Hard false: `visualization_ready`, `visual_correspondence_verified`, `pde_validated`, `source_correspondence_verified`, `paper_exact`, `openai_field_identified`, `blowup_proved`.

The export uses no OpenAI numerical velocity, hidden frame time, hidden streamline seeds, camera registration, image fit, visual acceptance threshold, free `f=R(u,p)`, collapsed-velocity route, or threshold relaxation.

## Direct contribution to final `[u,v,w]`

ST048-S is currently the better-supported visualization-oriented member of the ST048 residual challengers, but it lived only behind an open research recipe. This increment makes the exact frozen challenger available as a standard multi-time Cartesian velocity artifact for immediate MATLAB/Python 3-D rendering and for the already-owned streamline/vorticity/VTK comparison lanes, while retaining an explicit identity and truth boundary.

It does not claim ST048-S visually matches the OpenAI public field. The next useful step is to consume this fixed grid or the exact callable candidate in one fixed 3-D render/streamline comparison rather than add another generic diagnostic.
