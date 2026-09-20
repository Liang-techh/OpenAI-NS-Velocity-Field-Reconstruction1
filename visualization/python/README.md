# ST054 Python visualization

This directory provides a small headless Python rendering entry for the same immutable
continuous ST054-Q2 / ST054-M3 snapshots used by the native MATLAB viewer.

From the repository root:

```bash
python -m pip install -e '.[dev]'
python visualization/python/st054_visualize.py \
  --model ST054-Q2 \
  --time 0.50 \
  --output outputs/st054_q2_t050.png \
  --receipt outputs/st054_q2_t050.json
```

The script calls the published checksum-bound `load_st054(...).velocity(x,y,z,t)`
directly. It does not interpolate a cached velocity grid and it does not fit a new
field. The left panel renders instantaneous bidirectional streamlines from a frozen
48-seed protocol; line color is speed. The right panel shows a centered-difference
`omega_z = d_x v - d_y u` slice at `z=0` as a visualization diagnostic.

## Frozen smoke protocol

- models: `ST054-Q2`, `ST054-M3`;
- CI times: `0.25`, `0.50`, `0.75`;
- streamline seeds: radii `.6/.9/1.2`, `z=+/- .3`, eight azimuths = 48 seeds;
- instantaneous streamline geometry: normalized-velocity, bidirectional vectorized RK4;
- step `.025`, 90 steps per direction;
- vorticity slice: `65 x 65` Cartesian samples on `[-1.6,1.6]^2`, `z=0`;
- fixed camera: elevation `18`, azimuth `-58`.

These numbers are rendering/engineering choices only. There is deliberately no image
similarity threshold and no numerical target extracted from an OpenAI image or video.
The Actions smoke only checks that both saved continuous fields can produce finite,
nonempty render data and PNGs while preserving the published checksum and false
scientific truth flags.

## External-method screening

The rendering layer uses Matplotlib through public APIs only; no Matplotlib source is
copied. The screened upstream repository was `matplotlib/matplotlib` at
`a0acad0712a8571509ce37b864e651255ff9a4fe` (Matplotlib License Agreement), classified
as **direct migration / public API only**.

SciPy `solve_ivp` / DOP853 was also screened at `scipy/scipy@2c2998c8d9d5fc63d3e6e3d330db46ee476379b2`
(BSD-3-Clause), but is **not adopted as the streamline integrator** in this increment.
The already-established fixed-step, normalized-velocity batch RK4 protocol is cheaper
to evaluate against this finite-basis continuous field and keeps the rendering smoke
deterministic. SciPy remains an inherited dependency of the published ST054 loader.

The seed/camera/normalized-streamline convention is a **suitable reimplementation of
an internal protocol** previously used by CR-A9-064 (`ee80de11775bd45157d114000378b6ee3422ff46`),
not an external or OpenAI-sourced numerical target.

## Truth boundary

A successful render means only that the saved continuous field is usable from Python
for this visualization path. It does **not** mean the candidate passed the registered
Navier--Stokes residual gates, matches the hidden OpenAI numerical field, has verified
visual correspondence, is paper-exact, or proves blow-up. In receipts,
`visualization_ready`, `visual_correspondence_verified`, `source_correspondence_verified`,
`pde_validated`, `paper_exact`, `openai_field_identified`, and `blowup_proved` remain
false. The repository-wide ST006 PDE baseline and CR001 scientific thresholds are not
changed by this visualization adapter.
