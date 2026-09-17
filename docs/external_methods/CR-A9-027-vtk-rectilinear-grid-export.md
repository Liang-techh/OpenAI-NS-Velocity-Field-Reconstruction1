# CR-A9-027 — VTK rectilinear velocity-grid export

## Task

Constrained Agent 9 minimal delivery/visualization increment. The goal is to move an already-frozen Cartesian velocity grid closer to direct 3-D visualization without changing the candidate or any scientific acceptance rule.

## External result screened

- source repository: `Kitware/VTK`
- screened commit: `68c140cae7f80a0514dcd121f5c9526dd19e43f7`
- public specification: VTK simple legacy file format, `RECTILINEAR_GRID`, `POINT_DATA`, and `VECTORS`
- upstream license: BSD 3-clause terms in `Copyright.txt`
- classification: **suitable to reimplement / format-contract only**
- copied upstream implementation: none
- dependency delta: none

The VTK legacy specification describes rectilinear geometry with monotonically increasing x/y/z coordinate arrays and point-associated vector attributes. The repository-local implementation writes only that public interchange contract. It does not vendor VTK, call private VTK APIs, or copy a VTK writer implementation.

## Migration scope

`export_velocity_grid_vtk(...)` accepts frozen arrays in repository layout `(time,x,y,z,component)` with component order `[u,v,w]`. It writes one ASCII legacy `RECTILINEAR_GRID` file per time slice, plus a JSON manifest binding candidate SHA, grid SHA, provenance, file hashes, time values, speed diagnostics, component order, and the VTK point-order convention.

VTK point tuples are emitted x-fastest, then y, then z. A scalar `speed` field is included alongside the vector `velocity` field so ParaView/VTK-compatible readers can color the same frozen field without reconstructing magnitude externally.

## Repository-specific differences and truth boundary

This is deliberately narrower than VTK's own writer stack:

- ASCII only; no binary, XML, VTKHDF, compression, parallel I/O, or arbitrary datasets;
- rectilinear Cartesian grids only;
- one file per time slice rather than hidden time interpolation;
- no resampling, extrapolation, camera fitting, streamline integration, vorticity reconstruction, threshold selection, or candidate selection;
- exact-zero time slices fail closed, but this is not a replacement for the registered `E(0.25)=1±0.001` nontriviality gate;
- exported/interpolated finite differences are not valid evidence for the registered PDE gate;
- successful export means only `velocity_export_ready=true` for this interchange artifact.

The manifest keeps `visualization_ready=false`, `visual_correspondence_verified=false`, `pde_validated=false`, `paper_exact=false`, `openai_field_identified=false`, and `blowup_proved=false`.

## Direct contribution to final `[u,v,w]`

The intended path is:

`frozen candidate -> sampled (t,x,y,z,[u,v,w]) grid -> governed VTK files -> ParaView / VTK-capable Python tooling / downstream 3-D visualization`.

This removes the need for downstream visualization software to reimplement Eq45, support tapers, temporal schedules, or candidate-specific analytic code merely to consume the frozen velocity vectors.

## Remaining limitation

The exporter does not decide which current candidate is final, does not establish visual correspondence to OpenAI public imagery, and does not assess spatial-grid resolution sensitivity. Final visualization promotion still requires a frozen candidate identity, common visual frame, public-observable comparison, and resolution-stable 3-D morphology while PDE/support status remains reported independently.
