# CR-A9-004 — truth-bounded VTK velocity-series export

## External method screening

- Classification: **suitable to reimplement minimally**.
- Source repository: `nschloe/meshio`.
- Screened commit: `b2ee99842e119901349fdeee06b5bf61e01f450a`.
- License: MIT (`LICENSE.txt` at the screened commit).
- Relevant public behavior: meshio supports legacy VTK output and its VTK writer emits standard VTK dataset headers; its README links the published VTK file-format specification.
- Migration scope: no meshio source code is copied and no new dependency is added. This repository reimplements only the tiny ASCII legacy `RECTILINEAR_GRID` subset needed for `x,y,z` coordinates, a 3-component point-data velocity vector, and scalar speed.
- Why not direct migration: meshio is a broad mesh-format package and the final artifact here is a structured velocity field, so adding the whole dependency would be disproportionate. The minimal writer keeps the visualization path inspectable and deterministic.
- Difference from the source implementation: this writer is intentionally restricted to one structured visualization schema, one vector field, scalar speed, and a truth-bounded JSON manifest. It is not a general mesh reader/writer.

## Direct contribution to the final velocity field

A frozen callable `velocity(points, t) -> [u,v,w]` can be sampled at declared times and exported as a deterministic `.vtk` series that ParaView and other VTK-aware tools can open directly. This removes a visualization handoff blocker without changing the candidate, optimizer, pressure, forcing, or validation rules.

## Truth boundary

`visualization_export_ready=true` means only that the callable was sampled and serialized in the declared VTK schema. It does **not** establish Navier–Stokes residual acceptance, divergence-free structure, support/boundary correctness, resemblance to OpenAI public media, hidden-field identity, paper exactness, or blow-up.

No domain, viscosity, time interval, forcing family, nontriviality normalization, residual norm, validation seed, or acceptance threshold is modified by this increment.
